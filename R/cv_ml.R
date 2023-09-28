#' @title cv_ml
#' @description CV (I-squared) for mulilevel meta-analytic models, based on Yang et al. (2023). Under multilevel models, we can have multiple CV. TODO - we need to cite original CV paper
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
#' CV_eng_1 <- cv_ml(english_MA, boot = 10)
#' CV_eng_2 <- cv_ml(english_MA)
#'
#' ## Fish example
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
#' random = list(~1 | group_ID, ~1 | es_ID),
#' mods = ~ experimental_design + trait.type + deg_dif + treat_end_days,
#' method = "REML", test = "t", data = warm_dat,
#' control=list(optimizer="optim", optmethod="Nelder-Mead"))
#' CV_fish_1 <- cv_ml(model, boot = 10)
#' CV_fish_2 <- cv_ml(model)
#'
#' # Lim example
#' data(lim)
#' # Add in the sampling variance
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' # Lets fit a meta-regression - I will do Article non-independence.
#' The phylogenetic model found phylogenetic effects, however, instead we could fit Phylum as a fixed effect and explore them with an Orchard Plot
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data=lim)
#' CV_lim_1 <- cv_ml(lim_MR, boot = 10)
#' CV_lim_2 <- cv_ml(lim_MR)
#' }
#' @references TODO
#' @export

cv_ml <- function(model,
                  boot = NULL) {

  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}

  if(any(model$tau2 > 0)) { stop("Sorry. At the moment cv_ml cannot take models with heterogeneous variance.")}

    CVs <- ml_cv(model)

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
      CV_each <- sapply(sim, function(ysim) {
        tmp <- tryCatch(metafor::rma.mv(ysim, vi,
                                        random = random_formula, data = data))
        pb$tick()
        Sys.sleep(1/boot)
        CV <- ml_cv(tmp)
      })
    } else{
      CV_each <- sapply(sim, function(ysim) {

        # The model
        tmp <- tryCatch(metafor::rma.mv( ysim, vi,
                                         mods = mods_formula,
                                         random = random_formula,
                                         data = data))
        pb$tick()
        Sys.sleep(1 / boot)
        CV <- ml_cv(tmp)
        return(CV) })
    }
    # Summarise the bootstrapped distribution.
    CVs_each_95 <- data.frame(t(apply(CV_each, 1, stats::quantile, probs=c(0.5, .025, .975))))
    CVs <-  round(CVs_each_95, digits = 3)
    colnames(CVs) = c("Est.", "2.5%", "97.5%")
  }

  return(CVs)
}


#' @title ml_cv
#' @description Calculated CV for mulilevel meta-analytic models
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @export

ml_cv <- function(model){

 # total cv
 CV_total <- sqrt(sum(model$sigma2)) / abs(model$beta[[1]])
 # cv at different levels
 CV_each <-  sqrt(model$sigma2) / abs(model$beta[[1]])
 names(CV_each) <- paste0("CV_", model$s.names)
 names(CV_total) <- "CV_Total"

 CVs <- c(CV_total, CV_each)

 return(CVs)
}

