#' @title cvh1_ml
#' @description CVH1 for mulilevel meta-analytic models, based on Yang et al. (2023). Under multilevel models, we can have multiple CVH1. TODO - we need to cite original CVH1 paper
#' @param model Model object of class \code{rma.mv} or \code{rma}. Currently only model objects using the \code{mods} argument work (e.g., \code{mod = ~1}).
#' @param boot Number of simulations to run to produce 95 percent confidence intervals for I2. Default is \code{NULL}, where only the point estimate is provided.
#' @return A data frame containing all the model results including mean effect size estimate, confidence, and prediction intervals
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' # IMPORTANT NOTE ** boot = 10 is set LOW deliberately to make the models run fast. You should always run for at least boot = 1000
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl,
#' sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E,
#' m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD,
#' random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' CVH1_eng_1 <- cvh1_ml(english_MA, boot = 10)
#' CVH1_eng_2 <- cvh1_ml(english_MA)
#'
#' ## Fish example
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
#' random = list(~1 | group_ID, ~1 | es_ID),
#' mods = ~ experimental_design + trait.type + deg_dif + treat_end_days,
#' method = "REML", test = "t", data = warm_dat,
#' control=list(optimizer="optim", optmethod="Nelder-Mead"))
#' CVH1_fish_1 <- cvh1_ml(model, boot = 10)
#' CVH1_fish_2 <- cvh1_ml(model)
#'
#' # Lim example
#' data(lim)
#' # Add in the sampling variance
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' # Lets fit a meta-regression - I will do Article non-independence.
#' The phylogenetic model found phylogenetic effects, however, instead we could fit Phylum as a fixed effect and explore them with an Orchard Plot
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data=lim)
#' CVH1_lim_1 <- cvh1_ml(lim_MR, boot = 10)
#' CVH1_lim_2 <- cvh1_ml(lim_MR)
#' }
#' @references TODO
#' @export

cvh1_ml <- function(model,
                  boot = NULL) {

  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}

  if(any(model$tau2 > 0)) { stop("Sorry. At the moment cvh1_ml cannot take models with heterogeneous variance.")}

    CVH1s <- ml_cvh1(model)

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
      CVH1_each <- sapply(sim, function(ysim) {
        tmp <- tryCatch(metafor::rma.mv(ysim, vi,
                                        random = random_formula, data = data))
        pb$tick()
        Sys.sleep(1/boot)
        CVH1 <- ml_cvh1(tmp)
      })
    } else{
      CVH1_each <- sapply(sim, function(ysim) {

        # The model
        tmp <- tryCatch(metafor::rma.mv( ysim, vi,
                                         mods = mods_formula,
                                         random = random_formula,
                                         data = data))
        pb$tick()
        Sys.sleep(1 / boot)
        CVH1 <- ml_cvh1(tmp)
        return(CVH1) })
    }
    # Summarise the bootstrapped distribution.
    CVH1s_each_95 <- data.frame(t(apply(CVH1_each, 1, stats::quantile, probs=c(0.5, .025, .975))))
    CVH1s <-  round(CVH1s_each_95, digits = 3)
    colnames(CVH1s) = c("Est.", "2.5%", "97.5%")
  }

  return(CVH1s)
}


#' @title ml_cvh1
#' @description Calculated CVH1 for mulilevel meta-analytic models
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @export

ml_cvh1 <- function(model){

 # total cvh1
 CVH1_total <- sqrt(sum(model$sigma2)) / abs(model$beta[[1]])
 # cvh1 at different levels
 CVH1_each <-  sqrt(model$sigma2) / abs(model$beta[[1]])
 names(CVH1_each) <- paste0("CVH1_", model$s.names)
 names(CVH1_total) <- "CVH1_Total"

 CVH1s <- c(CVH1_total, CVH1_each)

 return(CVH1s)
}

