# Auxiliary functions 

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

.is_metafor_object <- function(obj) {
    if (!inherits(obj, c("rma.mv", "rma", "rma.uni"))) {
        stop(sprintf("'%s' is not a metafor object.", deparse(substitute(obj))),
             call. = FALSE)

    }
}

# A group is valid if:
#   - Is a column from data.
#   - Is a categorical column
#   - Has not NAs
.is_group_valid <- function(data, column_name) {
    if (!.column_exists(data, column_name)) {
        stop(sprintf("Column '%s' is not in '%s'",
                     column_name, deparse(substitute(data))),
             call. = FALSE)
    }
    
    col <- data[[column_name]]
    
    if (!.column_is_categorical(col)) {
        stop(sprintf("Column '%s' isn't categorical. Group variable must be categorical.", column_name),
             call. = FALSE)
    } else if (anyNA(col)) {
        stop(sprintf("Column '%s' has %d NAs. Group variable can't have NAs",
                     column_name, sum(is.na(col))),
             call. = FALSE)
    }
}
