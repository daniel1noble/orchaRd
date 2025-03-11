# Auxiliary functions 

# The basic data structure is a data frame, and arguments are usually
# columns from this data frame. 

#' Validate Group Argument
#'
#' Checks if a group variable meets the requirements for analysis:
#' - It is not missing or NULL
#' - It exists in the data frame
#' - It is a categorical variable
#' - It has no missing values
#'
#' @param data data.frame. A data frame containing the group variable
#' @param column_name character. The name of the column to check
#'
#' @return logical. TRUE if the group variable is valid, otherwise stops with an error message
#' @keywords internal
.is_group_valid <- function(data, column_name) {
  if (missing(column_name) || is.null(column_name)) {
    stop("Please specify the 'group' argument by providing the name of the grouping variable.",
         call. = FALSE)
  }

  if (!.column_exists(data, column_name)) {
    stop(sprintf("Column '%s' is not in '%s'", column_name, deparse(substitute(data))),
         call. = FALSE)
  } else {
    col <- data[[column_name]]
  }
  
  if (!.column_is_categorical(col)) {
    stop(sprintf("Column '%s' isn't categorical. Group variable must be categorical.",
                 column_name),
         call. = FALSE)
  } 

  if (anyNA(col)) {
    stop(sprintf("Column '%s' has %d NAs. Group variable can't have NAs",
                 column_name, sum(is.na(col))),
         call. = FALSE)
  }

  return(TRUE)
}

#' Check if Column Exists
#'
#' Verifies if a specified column exists in the data frame
#'
#' @param data data.frame. A data frame to check
#' @param column_name character. The name of the column to look for
#'
#' @return logical. TRUE if the column exists, FALSE otherwise
#' @keywords internal
.column_exists <- function(data, column_name) {
    if (column_name %in% names(data)) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

#' Check if Column is Categorical
#'
#' Determines if a column can be used as a categorical variable.
#' Valid categorical variables include:
#' - Character vectors
#' - Factors
#' - Integer values (common for ID numbers used as factors)
#'
#' @param col vector. A vector to check (character, factor, or numeric)
#'
#' @return logical. TRUE if the column can be used as a category, FALSE otherwise
#' @keywords internal
.column_is_categorical <- function(col) {
    # Check if a column can be used as a category.
    # It accepts characters, factors and, because the use of ID numbers 
    # is common, numbers used as 'numeric factors'.
    if (is.character(col) || is.factor(col) || (is.numeric(col) && all(col %% 1 == 0))) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}


# ---------------------------------------------------------------------
# Validate 'model' argument

#' Validate Model Argument
#'
#' Checks if a model argument is valid for meta-analysis:
#' - It is not missing or NULL
#' - It is a metafor package object (rma.mv, rma, etc.)
#'
#' @param model object. A model object from the metafor package (rma.mv, rma, rma.uni, or robust.rma)
#'
#' @return logical. TRUE if the model is valid, otherwise stops with an error message
#' @keywords internal
.is_model_valid <- function(model) {
  if (missing(model) || is.null(model)) {
    stop("Incorrect argument 'model'. Please specify the 'model' argument by providing rma.mv or rma model object.",
         call. = FALSE)
  } 

  if (!.is_metafor_object(model)) {
    stop("Incorrect argument 'model'. Please specify the 'model' argument by providing rma.mv or rma model object.",
         call. = FALSE)
  }

  return(TRUE)
}

#' Check if Object is from metafor Package
#'
#' Verifies if an object is of class rma.mv, rma, rma.uni, or robust.rma
#'
#' @param obj object. An R object to check for metafor class inheritance
#'
#' @return logical. TRUE if the object is from metafor package, FALSE otherwise
#' @keywords internal
.is_metafor_object <- function(obj) {
  if (!inherits(obj, c("rma.mv", "rma", "rma.uni", "robust.rma"))) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}
