#' Leave-One-Out Analysis for Meta-Analytic Models
#'
#' Performs a leave-one-out (LOO) analysis on a meta-analytic model from
#' the **metafor** package by iteratively removing each level of a grouping
#' variable and refitting the model.
#'
#' @param model A meta-analytic model fitted using **metafor**.
#' @param group A factor or categorical variable specifying the leave-one-out groups.
#'
#' @return A data frame with:
#'   - `left_out`: The removed group.
#'   - `b`: The estimated coefficient.
#'   - `ci_lb`: The lower confidence bound.
#'   - `ci_ub`: The upper confidence bound.
#'
#' @details Uses `run_leave1out()` for model refitting and `models_to_dataframe()`
#'   to structure the results.
#'
#' @seealso [metafor::rma()]
#'
#' @examples
#' \dontrun{
#' res <- metafor::rma.mv(lnrr, lnrr_vi, random = ~ 1 | paper_ID, data = fish)
#' loo_results <- leave_one_out(res, group = "paper_ID")
#' }
#'
#' @export

leave_one_out <- function(model, group, vcalc = NULL) {
  # Validate inputs
  if (!inherits(model, c("rma.mv", "rma", "rma.uni"))) {
    stop("Model must be a metafor rma object", call. = FALSE)
  }
  if (!(group %in% names(model$data))) {
    stop(sprintf("Group variable '%s' not found in model data", group), call. = FALSE) 
  }

  #
  # TODO: Should validate vcalc here!!
  #
  if (!is.null(vcalc)) {
    if (!is.list(vcalc)) {
      stop("vcalc must be a list with the arguments for the 'vcalc' function: e.g., vcalc = list(vi = lnrr_vi, cluster = 'paper_ID', obs = 'es_ID', rho = 0.5)",
           call. = FALSE)
    }
  }

  models_outputs <- run_leave1out(model, group, vcalc)
  loo_dataframe <- models_to_dataframe(models_outputs)
  return(loo_dataframe)
}


#' Fit Meta-Analytic Models While Omitting Each Group
#'
#' Iteratively refits a meta-analytic model, leaving out one level of a specified 
#' grouping variable in each iteration.
#'
#' @param model A fitted model object containing a \code{data} element, which must be 
#'   a data frame with all model variables.
#' @param group A character string specifying the column in \code{model$data} that 
#'   defines the groups to be omitted one at a time.
#'
#' @details The function removes each unique group from \code{model$data} one at a time, 
#'   refitting the model using \code{update()}. If an update fails, a warning is issued, 
#'   and \code{NULL} is returned for that group.
#'
#' @return A named list of models, each fitted after omitting one group. Names correspond 
#'   to the omitted group IDs.
#'
#' @keywords internal

run_leave1out <- function(model, group, vcalc = NULL) {
  group_ids <- unique(model$data[[group]])
  if (length(group_ids) < 2) {
    stop("Need at least 2 groups for leave-one-out analysis", call. = FALSE)
  }

  tmp_model <- model

  models_outputs <- lapply(group_ids, function(id_left_out) {
    new_data <- subset(model$data, model$data[[group]] != id_left_out)

    tmp_res <- tryCatch({
      # If vcalc is provided, calculate the VCV matrix
      if (!is.null(vcalc)) {
        #
        #
        # TODO: This part is a mess, should be cleaned up
        # WARNING: 
        # There is no validation of the vcalc arguments, so careful with that axe eugene!
        #
        tmp_VCV <- vcalc(vi      = new_data[[vcv_args$vi]],
                         cluster = new_data[[vcv_args$cluster]],
                         obs     = new_data[[vcv_args$obs]],
                         data    = new_data,
                         rho     = vcv_args$rho)
        update(tmp_model, data = new_data, V = tmp_VCV)

      } else {
        # If no variance-covariance matrix is needed, just update the model
        update(tmp_model, data = new_data)
      }
    }, error = function(e) {
      warning(sprintf("Error fitting model when leaving out '%s': %s",
                      id_left_out,
                      e$message))
      return(NULL)
    })

    tmp_res
  })

  names(models_outputs) <- group_ids
  models_outputs
}


#' Convert List Of Model Outputs To Data Frame
#'
#' Extracts the first coefficient and its confidence interval
#' from a list of model outputs and converts them into a structured data frame.
#'
#' @param models_outputs A named list of models, where each model is expected to 
#'   have elements `b`, `ci.lb`, and `ci.ub` representing the coefficient and
#'   its confidence interval.
#'
#' @return A data frame with the following columns:
#'   - `left_out`: Names of the models (from `models_outputs`).
#'   - `b`: The extracted coefficient.
#'   - `ci_lb`: The lower bound of the confidence interval.
#'   - `ci_ub`: The upper bound of the confidence interval.
#'
#' @details The function assumes that each model in `models_outputs` has
#'   a structure with elements `b`, `ci.lb`, and `ci.ub`. The function extracts
#'   the first coefficient and its corresponding confidence interval and
#'   returns them in a tidy format.
#' @keywords internal

models_to_dataframe <- function(models_outputs) {

  # NOTE: Must transpose to get the right shape
  coef_matrix <- t(
    sapply(models_outputs, function(model) {
      c(b     = model$b[1],
        ci_lb = model$ci.lb[1],
        ci_ub = model$ci.ub[1])
    })
  )

  # Create dataframe with proper column names
  outputs_dataframe <- data.frame(
    left_out = names(models_outputs),
    coef_matrix,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  return(outputs_dataframe)
}
