# Auxiliary functions 

# The basic data structure is a data frame, and arguments are usually
# columns from this data frame. 

# ---------------------------------------------------------------------
# Validate 'group' argument
#
# A 'group' argument is required in most functions. 
# It is valid if:
# - It is not missing or NULL
# - It is a column in the data frame
# - It is a categorical variable
# - It has no missing values

.is_group_valid <- function(data, column_name) {
  if (missing(column_name) || is.null(column_name)) {
    stop("Please specify the 'group' argument by providing the name of the grouping variable.",
         call. = FALSE)
  }

  if (!.column_exists(data, column_name)) {
    stop(sprintf("Column '%s' is not in '%s'", column_name, deparse(substitute(data))),
         call. = FALSE)
  }
  
  col <- data[[column_name]]
  
  if (!.column_is_categorical(col)) {
    stop(sprintf("Column '%s' isn't categorical. Group variable must be categorical.",
                 column_name),
         call. = FALSE)
  } else if (anyNA(col)) {
    stop(sprintf("Column '%s' has %d NAs. Group variable can't have NAs",
                 column_name, sum(is.na(col))),
         call. = FALSE)
  }

  return(TRUE)
}

.column_exists <- function(data, column_name) {
    if (column_name %in% names(data)) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

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

.is_metafor_object <- function(obj) {
  if (!inherits(obj, c("rma.mv", "rma", "rma.uni", "robust.rma"))) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

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


