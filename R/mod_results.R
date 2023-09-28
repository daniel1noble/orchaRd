############# Key functions #############

#' @title mod_results
#' @description Using a \pkg{metafor} model object of class \code{rma} or \code{rma.mv}, this function creates a table of model results containing the mean effect size estimates for all levels of a given categorical moderator, and their corresponding confidence and prediction intervals. The function is capable of calculating marginal means from meta-regression models, including those with multiple moderator variables of mixed types (i.e. continuous and categorical variables).
#' @param model \code{rma.mv} model object
#' @param mod Moderator variable of interest that one wants marginal means for. Defaults to the intercept, i.e. \code{"1"}.
#' @param group The grouping variable that one wishes to plot beside total effect sizes, k. This could be study, species, or any grouping variable one wishes to present sample sizes for.
#' @param by Character vector indicating the name that predictions should be conditioned on for the levels of the moderator.
#' @param at List of levels one wishes to predict at for the corresponding variables in \code{by}. Used when one wants marginalised means. This argument can also be used to suppress levels of the moderator when argument \code{subset = TRUE}. Provide a list as follows: \code{list(mod = c("level1", "level2"))}.
#' @param weights How to marginalize categorical variables. The default is \code{weights = "prop"}, which weights moderator level means based on their proportional representation in the data. For example, if "sex" is a moderator, and males have a larger sample size than females, then this will produce a weighted average, where males are weighted more towards the mean than females. This may not always be ideal. In the case of sex, for example, males and females are roughly equally prevalent in a population. As such, you can give the moderator levels equal weight using \code{weights = "equal"}.
#' @param subset Used when one wishes to only plot a subset of levels within the main moderator of interest defined by \code{mod}. Default is \code{FALSE}, but use \code{TRUE} if you wish to subset levels of a moderator plotted (defined by \code{mod}) for plotting. Levels one wishes to plot are specified as a list, with the level names as a character string in the \code{at} argument. For subsetting to work, the \code{at} argument also needs to be specified so that the \code{mod_results} function knows what levels one wishes to plot.
#' @param N The name of the column in the data specifying the sample size so that each effect size estimate is scaled to the sample size, N. Defaults to \code{NULL}, so that precision is used for scaling each raw effect size estimate instead of sample size.
#' @param upper Logical, defaults to \code{TRUE}, indicating that the first letter of the character string for the moderator variable should be capitalized.
#' @param ... Additional arguments passed to \code{emmeans::emmeans()}.
#' @return A data frame containing all the model results including mean effect size estimate, confidence and prediction intervals
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples \dontrun{
#' # Simple eklof data
#' data(eklof)
#' eklof<-metafor::escalc(measure="ROM", n1i=N_control, sd1i=SD_control,
#' m1i=mean_control, n2i=N_treatment, sd2i=SD_treatment, m2i=mean_treatment, data = eklof)
#' # Add the unit level predictor
#' eklof$Datapoint<-as.factor(seq(1, dim(eklof)[1], 1))
#' # fit a MLMR - accouting for some non-independence
#' eklof_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~ Grazer.type, random=list(~1|ExptID,
#' ~1|Datapoint), data = eklof)
#' results <- mod_results(eklof_MR, mod = "Grazer.type", group = "ExptID")
#'
#' # Fish example demonstrating marginalised means
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
#' random = list(~1 | group_ID, ~1 | es_ID),
#' mods = ~ experimental_design + trait.type + deg_dif + treat_end_days,
#' method = "REML", test = "t",
#' control=list(optimizer="optim", optmethod="Nelder-Mead"), data = warm_dat)
#'   overall <- mod_results(model, group = "group_ID")
#' across_trait <- mod_results(model, group = "group_ID", mod = "trait.type")
#' across_trait_by_degree_diff <- mod_results(model, group = "group_ID",
#' mod = "trait.type", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif")
#' across_trait_by_degree_diff_at_treat_end_days10 <- mod_results(model, group = "group_ID",
#' mod = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = 10),
#' by = "deg_dif",data = warm_dat)
#' across_trait_by_degree_diff_at_treat_end_days10And50 <- mod_results(model, group = "group_ID",
#' mod = "trait.type", at = list(deg_dif = c(5, 10, 15),
#'  treat_end_days = c(10, 50)), by = "deg_dif")
#' across_trait_by_treat_end_days10And50 <- mod_results(model, group = "group_ID",
#' mod = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)),
#' by = "treat_end_days")
#' across_trait_by_treat_end_days10And50_ordinaryMM <- mod_results(model, group = "group_ID",
#' mod = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)),
#' by = "treat_end_days", weights = "prop")
#'
#' # Fish data example with a heteroscedastic error
#' model_het <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list(~1 | group_ID, ~1 + trait.type| es_ID), mods = ~ trait.type + deg_dif, method = "REML", test = "t", rho = 0, struc = "HCS", control=list(optimizer="optim", optmethod="Nelder-Mead"), data = warm_dat)
#' HetModel <- mod_results(model_het, group = "group_ID", mod = "trait.type", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif", weights = "prop")
#' orchard_plot(HetModel, xlab = "lnRR")
#' }
#' @export
#'
#'
# We will need to make sure people use "1" or"moderator_names"

mod_results <- function(model, mod = "1", group,  N = NULL,  weights = "prop", by = NULL, at = NULL, subset = FALSE, upper = TRUE, ...){

  if(any(grepl("-1|0", as.character(model$formula.mods)))){
    warning("It is recommended that you fit the model with an intercept. Unanticipated errors can occur otherwise.")
  }

  if(any(model$struct %in% c("GEN", "HCS"))){
    warning("We noticed you're fitting an ~inner|outer rma model ('random slope'). There are circumstances where the prediction intervals for such models are calculated incorrectly. Please check your results carefully.")
  }

  if(missing(model)){
    stop("Please specify the 'model' argument by providing rma.mv or rma model object. See ?mod_results")
  }

  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class rma.mv, rma, or robust.rma")}

  if(missing(group)){
    stop("Please specify the 'group' argument by providing the name of the grouping variable. See ?mod_results")
  }


  if(is.null(stats::formula(model))){ ##**NOTE** Not sure we need this bit of code anymore. Left here for now
    #model <- stats::update(model, "~1")
    model$formula.mods <- ~ 1
    #dat_tmp <- model$data$`1` <- "Intrcpt"
    #model$data <- dat_tmp
  }

   if(model$test == "t"){
    df_mod = as.numeric(model$ddf[[1]])
  } else{
    df_mod = 1.0e6 # almost identical to z value
  }

# Extract the data from the model object
data <- model$data 

# Check if missing values exist and use complete case data
if(any(model$not.na == FALSE)){
	data <- data[model$not.na,]
}

  if(is.character(data[[mod]]) | is.factor(data[[mod]]) | is.null(data[[mod]])) {
    grid <- emmeans::qdrg(formula = stats::formula(model), at = at, data = data, coef = model$b,
                          vcov = stats::vcov(model), df = model$k-1) ## NOTE: Added data argument emmeans >vers 1.7.4. Object is unstable so feeding in the relevant arguments from model object directly. Note, we should think about df!
    mm <- emmeans::emmeans(grid, specs = mod, df = df_mod, by = by, weights = weights, ...)

    # getting prediction intervals
    mm_pi <- pred_interval_esmeans(model, mm, mod = mod)

    if(is.null(by)){
      mod_table <- data.frame(name = firstup(as.character(mm_pi[,1]), upper = upper),
                              estimate = mm_pi[,"emmean"],
                              lowerCL = mm_pi[,"lower.CL"],
                              upperCL = mm_pi[,"upper.CL"],
                              lowerPR = mm_pi[,"lower.PI"],
                              upperPR = mm_pi[,"upper.PI"])

    } else{
      mod_table <- data.frame(name = firstup(as.character(mm_pi[,1]), upper = upper),
                              condition = mm_pi[,2], estimate = mm_pi[,"emmean"],
                              lowerCL = mm_pi[,"lower.CL"],
                              upperCL = mm_pi[,"upper.CL"],
                              lowerPR = mm_pi[,"lower.PI"],
                              upperPR = mm_pi[,"upper.PI"])
    }

    # Extract data
    data2 <- get_data_raw(model, mod, group, N, at = at, subset)

    mod_table$name <- factor(mod_table$name,
                             levels = mod_table$name,
                             labels = mod_table$name)

  } else{
    at2 <- list(mod = seq(min(data[,mod], na.rm = TRUE), max(data[,mod], na.rm = TRUE), length.out = 100))
    names(at2) <- mod
    grid <- emmeans::qdrg(formula =  stats::formula(model), data = data, coef = model$b,
                          vcov = stats::vcov(model), df = model$k-1, at = c(at2, at))  # getting 100 points. Fixing this to make it more general
    mm <- emmeans::emmeans(grid, specs = mod, by = c(mod, by), weights = weights, df = df_mod)

    # getting prediction intervals
    mm_pi <- pred_interval_esmeans(model, mm, mod = mod)

    if(is.null(by)){
      mod_table <- data.frame(moderator = mm_pi[,1],
                              estimate = mm_pi[,"emmean"],
                              lowerCL = mm_pi[,"lower.CL"],
                              upperCL = mm_pi[,"upper.CL"],
                              lowerPR = mm_pi[,"lower.PI"],
                              upperPR = mm_pi[,"upper.PI"])
    } else{
      mod_table <- data.frame(moderator = mm_pi[,1],
                              condition = mm_pi[,2],
                              estimate = mm_pi[,"emmean"],
                              lowerCL = mm_pi[,"lower.CL"],
                              upperCL = mm_pi[,"upper.CL"],
                              lowerPR = mm_pi[,"lower.PI"],
                              upperPR = mm_pi[,"upper.PI"])
    }

    # extract data
    data2 <- get_data_raw_cont(model, mod, group, N, by = by)

  }


  output <- list(mod_table = mod_table,
                 data = data2)

  class(output) <- c("orchard", "data.frame")

  return(output)
}




############# Key Sub-functions #############

#' @title pred_interval_esmeans
#' @description Function to get prediction intervals (credibility intervals) from \code{esmeans} objects (\pkg{metafor}).
#' @param model \code{rma.mv} object.
#' @param mm result from \code{emmeans::emmeans} object.
#' @param mod Moderator of interest.
#' @param ... other arguments passed to function.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @export


pred_interval_esmeans <- function(model, mm, mod, ...){

        tmp <- summary(mm)
        tmp <- tmp[ , ]
  test.stat <- stats::qt(0.975, tmp$df[[1]])

  if(length(model$tau2) <= 1 | length(model$gamma2) <= 1){ # Note this should fix #46 but code is repetitive and needs to be cleaned up. Other issue is how this plays with different rma. objects. uni models will treat slots for gamma NULL and we need to deal with this. 
                 sigmas <- sum(model$sigma2)
                 taus   <- model$tau2
                 gamma2 <- ifelse(is.null(model$gamma2), 0, model$gamma2)
                     PI <- test.stat * base::sqrt(tmp$SE^2 + sigmas + taus + gamma2)
        } else {
                 sigmas <- sum(model$sigma2)
                 taus   <- model$tau2
                 gammas <- model$gamma2
                  w_tau <- model$g.levels.k
                w_gamma <- model$g.levels.k

            if(mod == "1"){
                tau <- weighted_var(taus, weights = w_tau)
              gamma <- weighted_var(gamma, weights = w_gamma)
                     PI <- test.stat * sqrt(tmp$SE^2 + sigmas + tau + gamma)

            } else {
               PI <- test.stat * sqrt(tmp$SE^2 + sigmas + taus + gammas)
            }
        }

  tmp$lower.PI <- tmp$emmean - PI
  tmp$upper.PI <- tmp$emmean + PI

  # renaming "overall" to ""
  if(tmp[1,1] == "overall"){tmp[,1] <- "intrcpt"}

return(tmp)
}

#' @title get_data_raw
#' @description Collects and builds the data used to fit the \code{rma.mv} or \code{rma} model in \pkg{metafor}.
#' @param model \code{rma.mv} object.
#' @param mod the moderator variable.
#' @param group The grouping variable that one wishes to plot beside total effect sizes, k. This could be study, species, or whatever other grouping variable one wishes to present sample sizes.
#' @param N The name of the column in the data specifying the sample size, N. Defaults to \code{NULL}, so precision is plotted instead of sample size.
#' @param at List of moderators. If \code{at} is equal to \code{mod} then levels specified within \code{at} will be used to subset levels when \code{subset = TRUE}. Otherwise, it will marginalise over the moderators at the specified levels.
#' @param subset Whether or not to subset levels within the \code{mod} argument. Defaults to \code{FALSE}.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a data frame
#' @export
#' @examples \dontrun{
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list(~1 | group_ID, ~1 | es_ID), mods = ~ experimental_design + trait.type + deg_dif + treat_end_days, method = "REML", test = "t", data = warm_dat, control=list(optimizer="optim", optmethod="Nelder-Mead"))
#'  test <- get_data_raw(model, mod = "trait.type", group = "group_ID", at = list(trait.type = c("physiology", "morphology")))
#'  test2 <- get_data_raw(model, mod = "1", group = "group_ID")
#'
#'  data(english)
#'  # We need to calculate the effect sizes, in this case d
#'  english <- escalc(measure = "SMD", n1i = NStartControl, sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E, m2i = MeanE, var.names=c("SMD","vSMD"))
#'  model <- rma.mv(yi = SMD, V = vSMD, random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#'  test3 <-  get_data_raw(model, mod = "1", group = "StudyNo")}

get_data_raw <- function(model, mod, group, N = NULL, at = NULL, subset = TRUE){
  if(missing(group)){
    stop("Please specify the 'group' argument by providing the name of the grouping variable. See ?mod_results")
  }
 
# Extract the data from the model object
  data <- model$data 

# Check if missing values exist and use complete case data
  if(any(model$not.na == FALSE)){
    data <- data[model$not.na,]
  }

  if(!is.null(at) & subset){
    # Find the at slot in list that pertains to the moderator and extract levels
    at_mod <- at[[mod]]
    position2 <- which(data[,mod] %in% at_mod)
    # Subset the data to only the levels in the moderator
    data <- data[position2,]
    yi <- model$yi[position2]
    vi <- model$vi[position2]
    type <- attr(model$yi, "measure")
  } else {
    # Extract effect sizes
    yi <- model$yi
    vi <- model$vi
    type <- attr(model$yi, "measure")
  }
  if(mod == "1"){
    moderator <- "Intrcpt"
  }else{
    # Get moderator
    moderator <- as.character(data[[mod]]) # Could default to base instead of tidy
    moderator <- firstup(moderator)
  }
  # Extract study grouping variable to calculate the
  stdy <- data[[group]] # Could default to base instead of tidy
  data_reorg <- data.frame(yi, vi, moderator, stdy, type)
  #names(data_reorg)[4] <- "stdy" # sometimes stdy gets replaced by group's names
  row.names(data_reorg) <- 1:nrow(data_reorg)

  if(is.null(N) == FALSE){
    data_reorg$N <- data[ ,N]
  }

  return(data_reorg)
}

#' @title get_data_raw_cont
#' @description Collects and builds the data used to fit the \code{rma.mv} or \code{rma} model in \pkg{metafor} when a continuous variable is fit within a model object.
#' @param model \code{rma.mv} object.
#' @param mod the moderator variable.
#' @param group The grouping variable that one wishes to plot beside total effect sizes, k. This could be study, species or whatever other grouping variable one wishes to present sample sizes.
#' @param N  The name of the column in the data specifying the sample size, N. Defaults to \code{NULL} so that precision is plotted instead of sample size.
#' @param by Character name(s) of the 'condition' variables to use for grouping into separate tables.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a data frame
#' @export

#TODO what if there is no "by"

get_data_raw_cont <- function(model, mod, group, N = NULL, by){
  if(missing(group)){
    stop("Please specify the 'group' argument by providing the name of the grouping variable. See ?mod_results")
  }

  # Extract the data from the model object
  data <- model$data 

# Check if missing values exist and use complete case data
  if(any(model$not.na == FALSE)){
    data <- data[model$not.na,]
  }

  # Extract effect sizes
  yi <- model$yi
  vi <- model$vi
  type <- attr(model$yi, "measure")
  # Get moderator
  moderator <- data[[mod]] # Could default to base instead of tidy
  #names(moderator) <  "moderator"
  if(is.null(by)){
  condition <- data[ , by]
  }else{
    condition <- data[[by]]
  }
  #names(condition) <  "condition"
  # Extract study grouping variable to calculate the
  stdy <- data[[group]] # Could default to base instead of tidy
  data_reorg <- data.frame(yi, vi, moderator, condition, stdy, type)
  # if(!is.na(names(data_reorg)[names(data_reorg) == by]) == TRUE) {  ## FAILING HERE
  #   names(data_reorg)[names(data_reorg) == by] <- "condition"
  # }
  #names(data_reorg)[5] <- "stdy" # sometimes stdy gets replaced by group's names
  row.names(data_reorg) <- 1:nrow(data_reorg)

  if(is.null(N) == FALSE){
    data_reorg$N <- data[ ,N]
  }

  return(data_reorg)
}

############# Helper-functions #############

#' @title firstup
#' @description Uppercase moderator names
#' @param x a character string
#' @param upper logical indicating if the first letter of the character string should be capitalized. Defaults to TRUE.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a character string with all combinations of the moderator level names with upper case first letters
#' @export

firstup <- function(x, upper = TRUE) {
        if(upper){
        substr(x, 1, 1) <- toupper(substr(x, 1, 1))
        x
        } else{ x }
      }


#' @title print.orchard
#' @description Print method for class 'orchard'
#' @param x an R object of class orchard
#' @param ... Other arguments passed to print
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a data frame
#' @export
#'

print.orchard <- function(x, ...){
    return(print(x$mod_table))
}

#' @title weighted_var
#' @description Calculate weighted variance
#' @param x A vector of tau2s to be averaged
#' @param weights Weights, or sample sizes, used to average the variance
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a vector with a single weighted variance
#' @export
#'

weighted_var <- function(x, weights){
    weight_var <- sum(x * weights) / sum(weights)
    return(weight_var)
}


#' @title num_studies
#' @description Computes how many studies are in each level of categorical moderators of a \code{rma.mv} model object.
#' @param mod Character string describing the moderator of interest.
#' @param data Raw data from object of class "orchard"
#' @param group A character string specifying the column name of the study ID grouping variable.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a table with the number of studies in each level of all parameters within a \code{rma.mv} or \code{rma} object.
#' @export
#' @examples
#' \dontrun{data(fish)
#'warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list( ~1 | es_ID,~1 | group_ID), mods = ~experimental_design-1, method = "REML", test = "t", data = warm_dat, control=list(optimizer="optim", optmethod="Nelder-Mead"))
#' num_studies(model$data, experimental_design, group_ID)
#' }

num_studies <- function(data, mod, group){

  # Summarize the number of studies within each level of moderator
   table <- data        %>%
            dplyr::group_by({{mod}}) %>%
            dplyr::summarise(stdy = length(unique({{group}})))

   table <- table[!is.na(table$moderator),]
  # Rename, and return
    colnames(table) <- c("Parameter", "Num_Studies")
      return(data.frame(table))

}
