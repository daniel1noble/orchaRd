#' @title m1_ml
#' @description M1 for mulilevel meta-analytic models, based on Yang et al. (2023). Under multilevel models, we can have multiple M1 - TODO - we need to cite original M1 paper
#' @param model Model object of class \code{rma.mv} or \code{rma}. Currently only model objects using the \code{mods} argument work (e.g., \code{mod = ~1}).
#' @param boot Number of simulations to run to produce 95 percent confidence intervals for M1. Default is \code{NULL}, where only the point estimate is provided.
#' @return A data frame containing all the model results including mean effect size estimate, confidence, and prediction intervals
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \donttest{
#' library(metafor)
#' # NOTE: boot is set LOW here for speed; use boot >= 1000 in practice.
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl,
#'   sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E,
#'   m2i = MeanE, var.names = c("SMD", "vSMD"), data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD,
#'   random = list(~1 | StudyNo, ~1 | EffectID), data = english)
#' m1_ml(english_MA)
#' m1_ml(english_MA, boot = 10)
#' }
#' @references TODO
#' @export

m1_ml <- function(model,
                  boot = NULL) {

  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}

  if(any(model$tau2 > 0)) { stop("Sorry. At the moment m1_ml cannot take models with heterogeneous variance.")}

  M1s <- ml_m1(model)

  # Extract the data from the model object
  data <- model$data

  # Check if missing values exist and use complete case data
  if(any(model$not.na == FALSE)){
    data <- data[model$not.na,]
  }

  if(!is.null(boot)){
    # Simulate the vector of effect sizes
    sim <- metafor::simulate.rma(model, nsim=boot) # Add try catch here? DN

    # Get formula from model object.
    random_formula <- model$random
    mods_formula <- metafor::formula.rma(model, type = "mods") #in case moderators
    vi <- model$vi

    # Parametric bootstrap
    pb <- progress::progress_bar$new(total = boot,
                                     format = "Bootstrapping [:bar] :percent ETA: :eta",
                                     show_after = 0)
    if(is.null(mods_formula)){
      M1_each <- sapply(sim, function(ysim) {
        tmp <- tryCatch(metafor::rma.mv(ysim, vi,
                                        random = random_formula, data = data))
        pb$tick()
        Sys.sleep(1/boot)
        M1 <- ml_m1(tmp)
      })
    } else{
      M1_each <- sapply(sim, function(ysim) {

        # The model
        tmp <- tryCatch(metafor::rma.mv( ysim, vi,
                                         mods = mods_formula,
                                         random = random_formula,
                                         data = data))
        pb$tick()
        Sys.sleep(1 / boot)
        M1 <- ml_m1(tmp)
        return(M1) })
    }
    # Summarise the bootstrapped distribution.
    M1s_each_95 <- data.frame(t(apply(M1_each, 1, stats::quantile, probs=c(0.5, .025, .975))))
    M1s <-  round(M1s_each_95, digits = 3)
    colnames(M1s) = c("Est.", "2.5%", "97.5%")
  }

  return(M1s)
}


#' @title ml_m1
#' @description Calculated M1 for mulilevel meta-analytic models
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @return A named numeric vector of M1 heterogeneity values (a proportion-type index on the standard-deviation scale, bounded between 0 and 1). The first element, \code{M1_Total}, is the overall value across all random effects; each remaining element (named \code{M1_<level>} after a random-effect level of the model) is the value for that level.
#' @export

ml_m1 <- function(model){

  # total m1
  M1_total <- sqrt(sum(model$sigma2)) / (abs(model$beta[[1]]) + sqrt(sum(model$sigma2)))
  # m1 at different levels
  M1_each <-  sqrt(model$sigma2) / (abs(model$beta[[1]]) + sqrt(sum(model$sigma2)))
  names(M1_each) <- paste0("M1_", model$s.names)
  names(M1_total) <- "M1_Total"

  M1s <- c(M1_total, M1_each)

  return(M1s)
}

