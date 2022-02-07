#' @title r2_ml
#' @description R2 (R-squared) for mixed (mulitlevel) models, based on Nakagawa & Schielzeth (2013)
#' @param model Model object of class 'rma.mv', 'rma'
#' @return A data frame containing all the model results including mean effect size estimate, confidence and prediction intervals with estimates converted back to r
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

  R2 <- R2_calc(model)

  if(!is.null(boot)){
    # Simulate the vector of effect sizes
    sim <- simulate(model, nsim=boot)

    # Get formula from model object.
    random_formula <- model$random
    mods_formula <- formula(model, type = "mods") #in case moderators
    vi <- model$vi

    # Paramatric bootsrap
    R2 <- sapply(sim, function(ysim) {
      # The model
      tmp <- rma.mv( ysim, vi,
                     mods = mods_formula,
                     random = random_formula,
                     data = data)
      R2s <- R2_calc(tmp)
      return(R2s)
    })

    # Summarise the bootstrapped distribution.
    R2 <- data.frame(t(apply(R2, 1, quantile, probs=c(0.5, .025, .975))))
    R2 <-  round(R2, digits = 3)
    colnames(R2) = c("Est.", "2.5%", "97.5%")
}

return(R2)

}


R2_calc <- function(model){
  # fixed effect variance
  fix <- var(as.numeric(as.vector(model$b) %*% t(as.matrix(model$X))))

  # marginal
  R2m <- fix / (fix + sum(model$sigma2))

  # conditional
  R2c <- (fix + sum(model$sigma2) - model$sigma2[length(model$sigma2)]) /
    (fix + sum(model$sigma2))

  R2s <- c(R2_marginal = R2m, R2_conditional = R2c)
  return(R2s)
}
