############# Key functions #############

#' @title mod_results
#' @description Using a metafor model object of class rma or rma.mv it creates a table of model results containing the mean effect size estimates for all levels of a given categorical moderator, their corresponding confidence intervals and prediction intervals. Function can calculate marginal means from meta-regression models with single or multiple moderator variables that are both continuous or categorical.
#' @param model rma.mv object
#' @param mod Moderator variable of interest that one wants marginal means for. Defaults to intercept "1".
#' @param group The grouping variable that one wishes to plot beside total effect sizes, k. This could be study, species or whatever other grouping variable one wishes to present sample sizes.
#' @param by The 'condition' variable that one wishes to have the mean for the moderator vary.
#' @param at The 'condition' that one wishes to calculate the means at, but is not presented in output
#' @param data The data frame used to fit the rma.mv model object
#' @param weights how to marginalize categorical variables. The default is weights = "prop", which wights means for moderator levels based on their proportional representation in the data. For example, if "sex" is a moderator, and males have a larger sample size than females, then this will produce a weighted average, where males are weighted more towards the mean than females. This may not always be ideal. IN the case if sex, for example, males and females are roughly equally prevalent in a population. As such, you can give the moderator levels equal weight using weights = "equal".
#' @param ... Additonal arguments passed to emmeans::emmeans()
#' @return A data frame containing all the model results including mean effect size estimate, confidence and prediction intervals
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples \dontrun{
#' # Simple eklof data
#' data(eklof)
#' eklof<-metafor::escalc(measure="ROM", n1i=N_control, sd1i=SD_control,
#' m1i=mean_control, n2i=N_treatment, sd2i=SD_treatment, m2i=mean_treatment,
#' data=eklof)
#' # Add the unit level predictor
#' eklof$Datapoint<-as.factor(seq(1, dim(eklof)[1], 1))
#' # fit a MLMR - accouting for some non-independence
#' eklof_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~ Grazer.type, random=list(~1|ExptID,
#' ~1|Datapoint), data=eklof)
#' results <- mod_results(eklof_MR, mod = "Grazer.type", group = "ExptID", data=eklof)
#'
#' # Fish example demonstrating marginalised means
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list(~1 | group_ID, ~1 | es_ID), mods = ~ experimental_design + trait.type + deg_dif + treat_end_days, method = "REML", test = "t", data = warm_dat, control=list(optimizer="optim", optmethod="Nelder-Mead"))
#'   overall <- mod_results(model, group = "group_ID", data = warm_dat)
#' across_trait <- mod_results(model, group = "group_ID", mod = "trait.type", data = warm_dat)
#' across_trait_by_degree_diff <- mod_results(model, group = "group_ID", mod = "trait.type", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif", data = warm_dat)
#' across_trait_by_degree_diff_at_treat_end_days10 <- mod_results(model, group = "group_ID", mod = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = 10), by = "deg_dif",data = warm_dat)
#' across_trait_by_degree_diff_at_treat_end_days10And50 <- mod_results(model, group = "group_ID", mod = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), by = "deg_dif", data = warm_dat)
#' across_trait_by_treat_end_days10And50 <- mod_results(model, group = "group_ID", mod = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), by = "treat_end_days", data = warm_dat)
#' across_trait_by_treat_end_days10And50_ordinaryMM <- mod_results(model, group = "group_ID", mod = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), by = "treat_end_days", weights = "prop", data = warm_dat)
#'
#' # Fish data example with a heteroscedastic error
#' model_het <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list(~1 | group_ID, ~1 + trait.type| es_ID), mods = ~ trait.type + deg_dif, method = "REML", test = "t", rho = 0, struc = "HCS", data = warm_dat, control=list(optimizer="optim", optmethod="Nelder-Mead"))
#' HetModel <- mod_results(model_het, group = "group_ID", mod = "trait.type", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif", weights = "prop", data = warm_dat)
#' orchard_plot(HetModel, xlab = "lnRR", data = warm_dat)
#' }
#' @export
#'
#'
# We will need to make sure people use "1" or"moderator_names"

marginal_means <- function(model, mod = "1", group, data, weights = "prop", by = NULL, at = NULL, ...){

  if(all(class(model) %in% c("rma.mv", "rma")) == FALSE) {stop("Sorry, you need to fit a metafor model of class rma.mv or rma")}

  # Extract data
   data2 <- get_data_raw(model, mod, group, data)

      model$data <- data
     grid <- emmeans::qdrg(object = model, at = at)

    if(model$test == "t"){
      df_mod = as.numeric(model$ddf[[1]])
    } else{
      df_mod = 100000000 # almost indentifcal to z value
    }

       mm <- emmeans::emmeans(grid, specs = mod, df = df_mod, by = by, weights = weights, ...)
    mm_pi <- pred_interval_esmeans(model, mm, mod = mod)


    if(is.null(by)){
      mod_table <- data.frame(name = firstup(as.character(mm_pi[,1])),
                              estimate = mm_pi[,"emmean"],
                              lowerCL = mm_pi[,"lower.CL"],
                              upperCL = mm_pi[,"upper.CL"],
                              lowerPR = mm_pi[,"lower.PI"],
                              upperPR = mm_pi[,"upper.PI"])

    } else{
      mod_table <- data.frame(name = firstup(as.character(mm_pi[,1])),
                              condition = mm_pi[,2], estimate = mm_pi[,"emmean"],
                              lowerCL = mm_pi[,"lower.CL"],
                              upperCL = mm_pi[,"upper.CL"],
                              lowerPR = mm_pi[,"lower.PI"],
                              upperPR = mm_pi[,"upper.PI"])

    }

    mod_table$name <- factor(mod_table$name,
                             levels = mod_table$name,
                             labels = mod_table$name)

    output <- list(mod_table = mod_table,
                        data = data2)

    class(output) <- c("orchard", "data.frame")

  return(output)
}


############# Key Sub-functions #############

#' @title pred_interval_esmeans
#' @description Function to get prediction intervals (credibility intervals) from esmeans objects (metafor)
#' @param model rma.mv object
#' @param mm result from emmeans::emmeans object'
#' @param mod Moderator of interest
#' @param ... other arguments passed to function
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @export

# TODO I think we can make it general and get PI with hetero from predict function
# TODO we want to add gamma too
# TODO - if we model Hetero for a particular categorical variable then, we do not want to
# TODO - this is actually challenging to make it general (they can have tau and gamma but the moderator needs to be the same)
# TODO - warning for just one kind of categorical variables for taus and gammas

pred_interval_esmeans <- function(model, mm, mod, ...){

        tmp <- summary(mm)
  test.stat <- stats::qt(0.975, tmp$df[[1]])

  if(length(model$tau2) <= 1){
                 sigmas <- sum(model$sigma2)
                     PI <- test.stat * base::sqrt(tmp$SE^2 + sigmas)
        } else {
            sigmas <- sum(model$sigma2)
            taus   <- model$tau2
                 w <- model$g.levels.k

            if(mod == "1"){
              tau <- weighted_var(taus, weights = w)
                     PI <- test.stat * sqrt(tmp$SE^2 + sigmas + tau)

            } else {
               PI <- test.stat * sqrt(tmp$SE^2 + sigmas + taus)
            }
        }

  tmp$lower.PI <- tmp$emmean - PI
  tmp$upper.PI <- tmp$emmean + PI

  # renaming "overall" to ""
  if(tmp[1,1] == "overall"){tmp[,1] <- "intrcpt"}

return(tmp)
}

#' @title get_data_raw
#' @description Collects and builds the data used to fit the rma.mv or rma model in metafor
#' @param model rma.mv object
#' @param mod the moderator variable
#' @param group The grouping variable that one wishes to plot beside total effect sizes, k. This could be study, species or whatever other grouping variable one wishes to present sample sizes.
#' @param data The data frame used to fit the rma.mv model object
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a data frame
#' @export
#' @examples \dontrun{
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list(~1 | group_ID, ~1 | es_ID), mods = ~ experimental_design + trait.type + deg_dif + treat_end_days, method = "REML", test = "t", data = warm_dat, control=list(optimizer="optim", optmethod="Nelder-Mead"))
#'  test <- get_data_raw(model, mod = "trait.type", group = "group_ID", data = warm_dat)
#'  test2 <- get_data_raw(model, mod = "1", group = "group_ID", data = warm_dat)
#'
#'  data(english)
#'  # We need to calculate the effect sizes, in this case d
#'  english <- escalc(measure = "SMD", n1i = NStartControl, sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E, m2i = MeanE, var.names=c("SMD","vSMD"), data = english)
#'  model <- rma.mv(yi = SMD, V = vSMD, random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#'  test3 <-  get_data_raw(english_MA, mod = "1", group = "StudyNo", data = english)}

get_data_raw <- function(model, mod, group, data){

  if(missing(group)){
    stop("Please specify the 'group' argument by providing the name of the grouping variable. See ?marginal_means.")
  }

  if(missing(data)){
    stop("Please specify the 'data' argument by providing the data used to fit the model. See ?marginal_means")
  }

  # Extract data
  # Check first if missing data exists
  if(length(attr(model$X, "dimnames")[[1]]) > 0){
    # full model delete missing values so need to adjust
     position <- as.numeric(attr(model$X, "dimnames")[[1]])
         data <- data[position, ] }

    if(mod == "1"){
      moderator <- "Intrcpt"
    }else{
      # Get moderator
       moderator <- as.character(data[,mod]) # Could default to base instead of tidy
       moderator <- firstup(moderator)
    }

    # Extract study grouping variable to calculate the
      stdy <- data[,group] # Could default to base instead of tidy

    # Extract effect sizes
        yi <- model$yi
        vi <- model$vi
      type <- attr(model$yi, "measure")

  data_reorg <- data.frame(yi, vi, moderator, stdy, type)
  row.names(data_reorg) <- 1:nrow(data_reorg)
  return(data_reorg)
}

############# Helper-functions #############

#' @title firstup
#' @description Uppercase moderator names
#' @param x a character string
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a character string with all combinations of the moderator level names with upper case first letters
#' @export

firstup <- function(x) {
        substr(x, 1, 1) <- toupper(substr(x, 1, 1))
        x
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
#' @description Computes how many studies are in each level of categorical moderators of a rma.mv model object.
#' @param mod Character string describing the moderator of interest.
#' @param data Raw data from object of class "orchard"
#' @param group A character string specifying the column name of the study ID grouping variable.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a table with the number of studies in each level of all parameters within a rma.mv or rma object.
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

  # Rename, and return
    colnames(table) <- c("Parameter", "Num_Studies")
      return(data.frame(table))

}
