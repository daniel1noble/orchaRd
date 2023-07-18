#' @title r2_ml
#' @description R2 (R-squared) for mixed (mulitlevel) models, based on Nakagawa & Schielzeth (2013).
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @param data Data frame used to fit the \code{rma.mv} or \code{rma} model object
#' @param boot The number of parametric bootstrap iterations, if desired. Defaults to \code{NULL}. A setting of 1000 is recommended as a minimum number of iterations.
#' @return A data frame containing all model results, including: mean effect size estimate, confidence and prediction intervals, with estimates converted back to r.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @references Nakagawa, S, and Schielzeth, H. 2013. A general and simple method for obtaining R2 from generalized linear mixed‐effects models. *Methods in Ecology and Evolution* 4(2): 133-142.
#' @examples
#' \dontrun{
#' data(lim)
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data=lim)
#' R2 <- r2_ml(lim_MR,data=lim, boot = 10)
#' }
#' @export
r2_ml <- function(model, data, boot = NULL) {

  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}

  if(any(model$tau2 > 0)) { stop("Sorry. At the moment r2_ml cannot take models with heterogeneous variance.")}

  R2 <- R2_calc(model)

  if(!is.null(boot)){

    if(any(class(model) %in% c("robust.rma")) == TRUE){stop("Sorry, bootstrapping currently doesn't work for robust.rma objects. Please use rma.mv instead.")}
    # Simulate the vector of effect sizes
    sim <- metafor::simulate.rma(model, nsim=boot) # Add try catch here? DN

    # Get formula from model object.
    random_formula <- model$random
    mods_formula <- metafor::formula.rma(model, type = "mods") #in case moderators
    vi <- model$vi

    pb <- progress::progress_bar$new(total = boot,
                                     format = "Bootstrapping [:bar] :percent ETA: :eta",
                                     show_after = 0)
    # Parametric bootstrap
    R2 <- sapply(sim, function(ysim) {
      # The model
      tmp <- tryCatch(metafor::rma.mv( ysim, vi,
                     mods = mods_formula,
                     random = random_formula,
                     data = data))
      R2s <- R2_calc(tmp)
      pb$tick()
      Sys.sleep(1 / boot)
      return(R2s)
    })

    # Summarise the bootstrapped distribution.
    R2 <- data.frame(t(apply(R2, 1, stats::quantile, probs=c(0.5, .025, .975))))
    R2 <-  round(R2, digits = 3)
    colnames(R2) = c("Est.", "2.5%", "97.5%")
}

return(R2)

}

#' @title R2_calc
#' @description Calculated R2 (R-squared) for mixed (mulitlevel) models, based on Nakagawa & Schielzeth (2013).
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' data(lim)
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1,
#' random=list(~1|Article, ~1|Datapoint), data=lim)
#' R2 <- R2_calc(lim_MR)
#' }
#' @export

R2_calc <- function(model){
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  # fixed effect variance
  fix <- stats::var(as.numeric(as.vector(model$b) %*% t(as.matrix(model$X))))

  # marginal
  R2m <- fix / (fix + sum(model$sigma2))

  # conditional. Need to remove 'residual' variance; assume this is the sigma level with the largest k. Really the only way we can get that to work.
  R2c <- (fix + sum(model$sigma2) - model$sigma2[which(model$s.nlevels.f == max(model$s.nlevels.f))]) /
    (fix + sum(model$sigma2))

  R2s <- c(R2_marginal = R2m, R2_conditional = R2c)
  return(R2s)
}
