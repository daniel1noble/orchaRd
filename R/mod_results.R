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
