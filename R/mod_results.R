#' @title get_est
#' @description Function gets estimates from rma objects (metafor)
#' @param model rma.mv object 
#' @param mod the name of a moderator. If meta-analysis (i.e. no moderator, se mod = "Int")
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au 
#' @export

get_est <- function (model, mod) {
  name <- firstup(as.character(stringr::str_replace(row.names(model$beta), {{mod}}, "")))

  estimate <- as.numeric(model$beta)
  lowerCL <- model$ci.lb
  upperCL <- model$ci.ub 
  
  table <- tibble::tibble(name = factor(name, levels = name, labels = name), estimate = estimate, lowerCL = lowerCL, upperCL = upperCL)

  return(table)
}

#' @title get_pred
#' @description Function to get prediction intervals (crediblity intervals) from rma objects (metafor)
#' @param model rma.mv object 
#' @param mod the name of a moderator 
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @export

get_pred <- function (model, mod) {

  name <- firstup(as.character(stringr::str_replace(row.names(model$beta), {{mod}}, "")))
  len <- length(name)
  
  if(len != 1){
  newdata <- matrix(NA, ncol = len, nrow = len)
  
  for(i in 1:len) {
    # getting the position of unique case from X (design matrix)
    pos <- which(model$X[,i] == 1)[[1]]
    # I think this is the other way around but it is diag(len) so fine
    newdata[, i] <- model$X[pos,]
    }
  pred <- metafor::predict.rma(model, newmods = newdata)
  }
  else {
    pred <- metafor::predict.rma(model)
    }
  lowerPR <- pred$cr.lb
  upperPR <- pred$cr.ub 
  
  table <- tibble::tibble(name = factor(name, levels = name, labels = name), lowerPR = lowerPR, upperPR = upperPR)
  return(table)
}

#' @title pred_interval_esmeans
#' @description Function to get prediction intervals (crediblity intervals) from esmeans objects (metafor)
#' @param model rma.mv object 
#' @param esmeans result from emmeans::emmeans object
#' @param ... other arguments passed to function
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @export
pred_interval_esmeans <- function(model, mm, ...){

        tmp <- summary(mm)
  test.stat <- qt(0.975, tmp$df)
  
  if(length(model$tau2) <= 1){         
                 sigmas <- sum(model$sigma2)
                     PI <- test.stat * sqrt(tmp$SE^2 + sigmas)
        } else {
            sigmas <- sum(model$sigma2)
            taus   <- model$tau2
                 w <- model$g.levels.k
            
            if(pred == "1"){
              tau <- weighted_var(taus, weights = w)
                     PI <- test.stat * sqrt(tmp$SE^2 + sigmas + tau)

            } else { 
               PI <- test.stat * sqrt(tmp$SE^2 + sigmas + taus)
            }
        }
  
  tmp$lower.PI <- tmp$emmean - PI
  tmp$upper.PI <- tmp$emmean + PI

return(tmp)
}

#' @title marginalised_means
#' @description Function to to get marginalised means from met-regression models with single or multiple moderator variables that are both continuous or categorical.
#' @param model rma.mv object 
#' @param data data frame used to fit rma.mv model
#' @param pred predictor variable of interest that one wants marginalised means for. 
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @example \dontrun{
#'warm_dat <- data(fish)
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list(~1 | group_ID, ~1 | es_ID), mods = ~ experimental_design + trait.type + deg_dif + treat_end_days, method = "REML", test = "t", data = warm_dat,                               control=list(optimizer="optim", optmethod="Nelder-Mead")) 
#'   overall <- marginalised_means(model, data = warm_dat)
#' across_trait <- marginalised_means(model, data = warm_dat, pred = "trait.type")
#' across_trait_by_degree_diff <- marginalised_means(model, data = warm_dat, pred = "trait.type", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif")
#' across_trait_by_degree_diff_at_treat_end_days10 <- marginalised_means(model, data = warm_dat, pred = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = 10), by = "deg_dif")
#' across_trait_by_degree_diff_at_treat_end_days10And50 <- marginalised_means(model, data = warm_dat, pred = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), by = "deg_dif")
#' across_trait_by_treat_end_days10And50 <- marginalised_means(model, data = warm_dat, pred = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), by = "treat_end_days")
#' across_trait_by_treat_end_days10And50_ordinaryMM <- marginalised_means(model, data = warm_dat, pred = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), by = "treat_end_days", weights = "prop")
#' }
#' @export
#' 
marginalised_means <- function(model, data, pred = "1", by = NULL, at = NULL, ...){
     model$data <- data

     grid <- emmeans::qdrg(object = model, at = at)
       mm <- emmeans::emmeans(grid, specs = pred, df = model$df, by = by, ...)
    mm_pi <- pred_interval_esmeans(model, mm, pred = pred)


    if(is.null(by)){
      mod_table <- tibble::tibble(name = mm_pi[,1], estimate = mm_pi[,"emmean"], lowerCL = mm_pi[,"lower.CL"], upperCL = mm_pi[,"upper.CL"], lowerPI = mm_pi[,"lower.PI"], upperPI = mm_pi[,"upper.PI"])
    
    } else{
      mod_table <- tibble::tibble(name = mm_pi[,1], mod = mm_pi[,2], estimate = mm_pi[,"emmean"], lowerCL = mm_pi[,"lower.CL"], upperCL = mm_pi[,"upper.CL"], lowerPI = mm_pi[,"lower.PI"], upperPI = mm_pi[,"upper.PI"])
      
    }

    output <- list(mod_table = mod_table, 
                data = data)

    class(output) <- "orchard"

  return(output)
}


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

#' @title get_data
#' @description Collects and builds the data used to fit the rma.mv or rma model in metafor
#' @param model rma.mv object
#' @param mod the moderator variable
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a data frame
#' @export
#' 
get_data <- function(model, mod){
     X <- as.data.frame(model$X)
 names <- vapply(stringr::str_split(colnames(X), {{mod}}), function(x) paste(unique(x), collapse = ""), character(1L))

  moderator <- matrix(ncol = 1, nrow = dim(X)[1])

  for(i in 1:ncol(X)){
      moderator <- ifelse(X[,i] == 1, names[i], moderator)
  }
    moderator <- firstup(moderator)
    yi <- model$yi
    vi <- model$vi
  type <- attr(model$yi, "measure")

data <- data.frame(yi, vi, moderator, type)
return(data)

}

#' @title mod_results
#' @description Using a metafor model object of class rma or rma.mv it creates a table of model results containing the mean effect size estimates for all levels of a given categorical moderator, their corresponding confidence intervals and prediction intervals
#' @param model rma.mv object 
#' @param mod the name of a moderator; put "Int" if the intercept model (meta-analysis) or no moderators. 
#' @return A data frame containing all the model results including mean effect size estimate, confidence and prediction intervals
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{data(eklof) 
#' eklof<-metafor::escalc(measure="ROM", n1i=N_control, sd1i=SD_control,
#' m1i=mean_control, n2i=N_treatment, sd2i=SD_treatment, m2i=mean_treatment,
#' data=eklof)
#' # Add the unit level predictor
#' eklof$Datapoint<-as.factor(seq(1, dim(eklof)[1], 1))
#' # fit a MLMR - accouting for some non-independence
#' eklof_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~ Grazer.type-1, random=list(~1|ExptID,
#' ~1|Datapoint), data=eklof)
#' results <- mod_results(eklof_MR, mod = "Grazer.type")
#' }
#' @export

mod_results <- function(model, mod) { 

	if(all(class(model) %in% c("rma.mv", "rma")) == FALSE) {stop("Sorry, you need to fit a metafor model of class rma.mv or rma")}

  data <- get_data(model, mod)

	# Get confidence intervals
	CI <- get_est(model, mod)

	# Get prediction intervals
	PI <- get_pred(model, mod)

	model_results <- list(mod_table = cbind(CI, PI[,-1]), data = data)

	class(model_results) <- "orchard"

	return(model_results)

}
# TODO - I think we can improve `mod` bit?

#' @title print.orchard
#' @description Print method for class 'orchard'
#' @param object x an R object of class orchard
#' @param ... Other arguments passed to print
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a data frame
#' @export
#' 
print.orchard <- function(object, ...){
    return(object$mod_table)
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
