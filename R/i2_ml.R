#' @title i2_ml
#' @description I2 (I-squared) for mulilevel meta-analytic models, based on Nakagawa & Santos (2012). Under multilevel models, we can have a multiple I2 (see also Senior et al. 2016). Alternatively, the method proposed by Wolfgang Viechtbauer (http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate?s[]=multilevel) can also be used.
#' @param model Model object of class 'rma.mv', 'rma'
#' @param method Method used to calculate I2. Two options exist, Nakagawa & Santos ('ns') or Wolfgang Viechtbauer's method ("wv").
#' @param boot Number of simulations to run to produce 95% CI's for I2. Default is NULL and only point estimate is provided.
#' @return A data frame containing all the model results including mean effect size estimate, confidence and prediction intervals with estimates converted back to r
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl, sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E, m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD, random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' I2 <- i2_ml(english_MA)
#' }
#' @references Senior, A. M., Grueber, C. E., Kamiya, T., Lagisz, M., O’Dwyer, K., Santos, E. S. A. & Nakagawa S. 2016. Heterogeneity in ecological and evolutionary meta-analyses: its magnitudes and implications. *Ecology* 97(12): 3293-3299.
#' @references Nakagawa, S, and Santos, E.S.A. 2012. Methodological issues and advances in biological meta-analysis. *Evolutionary Ecology* 26(5): 1253-1274.
#' @export

i2_ml <- function(model, method = c("ns", "wv"), boot = NULL) {

  ## evaluate choices
  method <- match.arg(method)

  # Wolfgang Viechtbauer's method
  if (method == "wv") {
    W <- solve(model$V)
    X <- model.matrix(model)
    P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
    I2_total <- sum(model$sigma2) / (sum(model$sigma2) + (model$k - model$p) / sum(diag(P)))
    I2_each <- model$sigma2 / (sum(model$sigma2) + (model$k - model$p) / sum(diag(P)))
    names(I2_each) <- paste0("I2_", model$s.names)
    # or Nakagawa & Santos (2012); they usually produce identical values
  } else {
    # sigma2_v = typical sampling error variance
    sigma2_v <- sum(1 / model$vi) * (model$k - 1) / (sum(1 / model$vi)^2 - sum((1 / model$vi)^2))
    I2_total <- sum(model$sigma2) / (sum(model$sigma2) + sigma2_v) # s^2_t = total variance
    I2_each <- model$sigma2 / (sum(model$sigma2) + sigma2_v)
    names(I2_each) <- paste0("I2_", model$s.names)
  }

  # putting all together
  I2s <- c(I2_total = I2_total, I2_each)

  if(!is.null(boot)){
    # Simulate the vector of effect sizes
    sim <- simulate(model, nsim=boot)

    # Get formula from model object. This is needed for the function to work. Slightly tricky to generalise but doable with careful checks
    random_formula <- paste0("~ 1 | ", sapply(model$mf.r, function(x) names(x)), collapse = " , ")
      mods_formula <- formula(model, type = "mods") #in case moderators


     I2_total <- sapply(sim, function(ysim) { # Need to get this working with formula of model
      tmp <- rma.mv(ysim, vSMD, random = list(random_formula), data=english)
      100 * sum(tmp$sigma2) / (sum(tmp$sigma2) + tmp$vi)
    })

    apply(I2_total, 1, quantile, probs=c(.025, .975))
  }

  return(I2s)
}


