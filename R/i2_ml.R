#' @title i2_ml
#' @description I2 (I-squared) for mulilevel meta-analytic models, based on Nakagawa & Santos (2012). Under multilevel models, we can have multiple I2 (see also Senior et al. 2016). Alternatively, the method proposed by Wolfgang Viechtbauer can also be used.
#' @param model Model object of class \code{rma.mv} or \code{rma}. Currently only model objects using the \code{mods} argument work (e.g., \code{mod = ~1}).
#' @param method Method used to calculate I2. Two options exist: a ratio-based calculation proposed by Nakagawa & Santos (\code{"ratio"}), or Wolfgang Viechtbauer's matrix method (\code{"matrix"}).
#' @param boot Number of simulations to run to produce 95 percent confidence intervals for I2. Default is \code{NULL}, where only the point estimate is provided.
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
#' i2_ml(english_MA)
#' i2_ml(english_MA, method = "matrix")
#' i2_ml(english_MA, boot = 10)
#' }
#' @references Senior, A. M., Grueber, C. E., Kamiya, T., Lagisz, M., O’Dwyer, K., Santos, E. S. A. & Nakagawa S. 2016. Heterogeneity in ecological and evolutionary meta-analyses: its magnitudes and implications. Ecology 97(12): 3293-3299.
#'  Nakagawa, S, and Santos, E.S.A. 2012. Methodological issues and advances in biological meta-analysis.Evolutionary Ecology 26(5): 1253-1274.
#' @export

i2_ml <- function(model,
                  method = c("ratio", "matrix"),
                  boot = NULL) {

  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}

  if (!is.null(boot) && identical(model$backend, "glmmTMB")) {
    stop("Bootstrap refits are not supported for objects converted with glmmTMB_to_rma().",
         call. = FALSE)
  }

  if(any(model$tau2 > 0)) { stop("Sorry. At the moment i2_ml cannot take models with heterogeneous variance.")}

  ## evaluate choices
  method <- match.arg(method)

  if (method == "matrix") {
    # Wolfgang Viechtbauer's method
    I2s <- matrix_i2(model)
  } else {
    # Nakagawa & Santos (2012)
    I2s <- ratio_i2(model)
  }

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
      I2_each <- sapply(sim, function(ysim) {
        tmp <- tryCatch(metafor::rma.mv(ysim, vi,
                                        random = random_formula, data = data))
        pb$tick()
        Sys.sleep(1/boot)
        if (method == "matrix") {
          I2 <- matrix_i2(tmp)
        }
        else {
          I2 <- ratio_i2(tmp)
        }
        return(I2)
      })
    } else{
     I2_each <- sapply(sim, function(ysim) {

              # The model
             tmp <- tryCatch(metafor::rma.mv( ysim, vi,
                             mods = mods_formula,
                           random = random_formula,
                             data = data))
             pb$tick()
             Sys.sleep(1 / boot)

            if(method == "matrix"){
              I2 <- matrix_i2(tmp)
            } else {
              I2 <- ratio_i2(tmp)
            }

             return(I2) })
    }
      # Summarise the bootstrapped distribution.
       I2s_each_95 <- data.frame(t(apply(I2_each, 1, stats::quantile, probs=c(0.5, .025, .975))))
               I2s <-  round(I2s_each_95, digits = 3)
      colnames(I2s) = c("Est.", "2.5%", "97.5%")
  }

  return(I2s)
}

#' @title matrix_i2
#' @description Calculated I2 (I-squared) for mulilevel meta-analytic models, based on a matrix method proposed by Wolfgang Viechtbauer.
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @return A named numeric vector of I2 (I-squared) values expressed as percentages. The first element, \code{I2_Total}, is the total heterogeneity; each remaining element (named \code{I2_<level>} after a random-effect level of the model) is the heterogeneity attributable to that level.
#' @examples \donttest{
#' library(metafor)
#' data(english)
#' english <- escalc(
#'   measure = "SMD", n1i = NStartControl,
#'   sd1i = SD_C, m1i = MeanC,
#'   n2i = NStartExpt, sd2i = SD_E, m2i = MeanE,
#'   var.names = c("SMD", "vSMD"), data = english)
#' english_MA <- rma.mv(
#'   yi = SMD, V = vSMD,
#'   random = list(~1 | StudyNo, ~1 | EffectID),
#'   data = english)
#' matrix_i2(english_MA)
#' }
#' @export
matrix_i2 <- function(model){

  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}

    W <- solve(model$V)
    X <- stats::model.matrix(model)
    P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
    I2_total <- 100* (sum(model$sigma2) / (sum(model$sigma2) + (model$k - model$p) / sum(diag(P))))
    I2_each <- 100* (model$sigma2 / (sum(model$sigma2) + (model$k - model$p) / sum(diag(P))))
    names(I2_each) <- paste0("I2_", model$s.names)
    names(I2_total) <- "I2_Total"
    I2s <- c(I2_total, I2_each)
    return(I2s)
}


#' @title ratio_i2
#' @description I2 (I-squared) for mulilevel meta-analytic models based on Nakagawa & Santos (2012). Under multilevel models, we can have a multiple I2 (see also Senior et al. 2016).
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @return A named numeric vector of I2 (I-squared) values expressed as percentages. The first element, \code{I2_Total}, is the total heterogeneity; each remaining element (named \code{I2_<level>} after a random-effect level of the model) is the heterogeneity attributable to that level.
#' @examples \donttest{
#' library(metafor)
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl,
#' sd1i = SD_C, m1i = MeanC, n2i = NStartExpt,
#' sd2i = SD_E, m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD,
#' random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' ratio_i2(english_MA)
#' }
#' @export
ratio_i2 <- function(model){

  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}

  # sigma2_v = typical sampling error variance
  sigma2_v <- sum(1 / model$vi) * (model$k - 1) /
              (sum(1 / model$vi)^2 - sum((1 / model$vi)^2))

  # s^2_t = total variance
  I2_total <- 100 * (sum(model$sigma2) / (sum(model$sigma2) + sigma2_v))
   I2_each <- 100 * (model$sigma2 / (sum(model$sigma2) + sigma2_v))
  names(I2_each) <- paste0("I2_", model$s.names)
  names(I2_total) <- "I2_Total"

  I2s <- c(I2_total, I2_each)
  return(I2s)
}
