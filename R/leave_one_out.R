#' Leave-One-Out Analysis for Meta-Analytic Models
#'
#' Performs a leave-one-out (LOO) analysis on a meta-analytic model from
#' the **metafor** package by iteratively removing each level of a grouping
#' variable and refitting the model.
#'
#' @param model A meta-analytic model fitted using **metafor**.
#' @param group A factor or categorical variable specifying the leave-one-out groups.
#'
#' @return Same as `mod_results()`, but with the estimates from each model ran
#'   in the leave-one-out analysis, and the effect sizes from each model.
#'
#' @examples
#' \dontrun{
#' res <- metafor::rma.mv(lnrr, lnrr_vi, random = ~ 1 | paper_ID, data = fish)
#' loo_results <- leave_one_out(res, group = "paper_ID")
#' }
#'
#' @export

leave_one_out <- function(model, group, vcalc_args = NULL) {
  # Check model is a metafor object
  if (!inherits(model, c("rma.mv", "rma", "rma.uni"))) {
    stop("Model must be a metafor rma object", call. = FALSE)
  }
  # Check if group is in model data
  if (!(group %in% names(model$data))) {
    stop(sprintf("Group variable '%s' not found in model data", group), call. = FALSE) 
  }
  # Check if we have at least 2 groups
  if (length(unique(model$data[[group]])) < 2) {
  stop("Need at least 2 groups for leave-one-out analysis", call. = FALSE)
  }
  # Check if vcalc is provided. If so, validate the arguments
  if (!is.null(vcalc_args)) {
    validate_vcalc_args(model$data, vcalc_args)
  } 

  # Run leave-one-out analysis
  models_outputs <- run_leave1out(model, group, vcalc_args)
  # Extract estimates
  estimates      <- get_estimates(models_outputs, group)
  # Extract effect sizes from each run
  effect_sizes   <- get_effectsizes(models_outputs, group)

  # Immitates the output of mod_results.
  #   - mod_table: Here are are the estimates from each model ran
  #   - data:  The effect sizes from each model
  output <- list(mod_table = estimates, data = effect_sizes)
  class(output) <- c("orchard", "data.frame")

  return(output)
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

run_leave1out <- function(model, group, vcalc_args = NULL) {
  tmp_model <- model
  group_ids <- unique(model$data[[group]])

  models_outputs <- lapply(group_ids, function(id_left_out) {
    new_data <- subset(model$data, model$data[[group]] != id_left_out)

    tmp_res <- tryCatch({

      # If vcalc is provided, calculate the VCV matrix
      if (!is.null(vcalc_args)) {
        tmp_VCV <- vcalc(vi      = new_data[[vcalc_args$vi]],
                         cluster = new_data[[vcalc_args$cluster]],
                         obs     = new_data[[vcalc_args$obs]],
                         data    = new_data,
                         rho     = vcalc_args$rho)
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

#' Orchard Plot for Leave-One-Out Analysis
#'
#' Performs a leave-one-out analysis on a meta-analytic model and produces an orchard plot.
#'
#' @param model A meta-analytic model from the metafor package.
#' @param mod Character string specifying the model version (default "1").
#' @param group Character string naming the column in model$data to omit iteratively.
#' @param ylab Optional label for the y-axis.
#' @param vcalc_args Optional list of parameters for computing the variance-covariance matrix.
#' @param alpha Numeric value for plot element transparency (default 0.1).
#' @param angle Numeric value for x-axis label angle (default 0).
#' @param g Logical; whether to apply grouping in the plot (default FALSE).
#' @param ... Additional arguments passed to \code{orchard_plot}.
#'
#' @return A ggplot2 object displaying the leave-one-out analysis with reference lines
#'   for the original model's confidence limits.
#'
#' @examples
#' \dontrun{
#'   res <- metafor::rma.mv(lnrr, lnrr_vi, random = ~ 1 | paper_ID, data = fish)
#'   orchard_leave1out(res, group = "paper_ID", ylab = "Study left out")
#' }
#'
#' @export
orchard_leave1out <- function(model,
                              mod = "1", 
                              group,
                              ylab = NULL,
                              vcalc_args = NULL,
                              alpha = 0.1, 
                              angle = 0,
                              g = FALSE,
                              ...) {

  # Extract original model estimates (point estimate and confidence limits)
  orig_table <- mod_results(model, group = group)$mod_table
  orig_results <- orig_table[, c("estimate", "lowerCL", "upperCL"), drop = FALSE]
  
  # Run leave-one-out analysis. Each iteration omits one element from 'group'
  leave1out_results <- leave_one_out(model = model, group = group, vcalc_args = vcalc_args)

  # Plot each result and include reference lines for the original model's confidence limits.
  p <- orchard_plot(leave1out_results,
                    alpha = alpha,
                    angle = angle,
                    g = g,
                    ...)

  p <- p +
    ggplot2::xlab(ylab) +    # xlab uses ylab? Yes. It is confusing, but because of the coord_flip it is like this.
    ggplot2::geom_hline(yintercept = orig_results$lowerCL, linetype = 2, color = "red") +
    ggplot2::geom_hline(yintercept = orig_results$upperCL, linetype = 2, color = "red")

  return(p)
}


validate_vcalc_args <- function(model_data, vcalc_args) {
  if (!is.list(vcalc_args)) {
    stop("vcalc must be a list with the arguments for the 'vcalc' function: e.g., vcalc_args = list(vi = 'lnrr_vi', cluster = 'paper_ID', obs = 'es_ID', rho = 0.5)",
         call. = FALSE)
  }

  # Check if required arguments for vcalc are present
  if (!all(c("vi", "cluster", "obs", "rho") %in% names(vcalc_args))) {
    stop("vcalc_args must contain at list the following elements: 'vi', 'cluster', 'obs', 'rho'", call. = FALSE)
  }

  # TODO: This is not a compelte check, but it's a start
  # Check if the vcalc arguments are present in the model data
  if (is.null(model_data[[vcalc_args$vi]]) || is.null(model_data[[vcalc_args$cluster]]) || is.null(model_data[[vcalc_args$obs]])) {
    stop("One or more of the vcalc arguments are not found in the model data", call. = FALSE)
  }

  return(vcalc_args)
}


#' Get Leave-One-Out Model Estimates
#'
#' Extracts and combines the model estimates from each leave-one-out iteration.
#'
#' @param outputs A named list of model objects from leave-one-out analysis.
#' @param group A string specifying the grouping variable used in the analysis.
#'
#' @return A data frame of model estimates with an added column indicating the omitted group.
#'
#' @keywords internal
get_estimates <- function(outputs, group) {
   # Call `mod_results` for each model ran in the leave-one-out,
   # transform its output to a dataframe, and then rbind()  
   # to create a long data frame with the estimates of all the models.
    estimates <- do.call(rbind, lapply(names(outputs), function(name) {
        res <- mod_results(outputs[[name]], group = group)
        df <- res$mod_table
        df$name <- name
        df
    }))

    row.names(estimates) <- NULL
    return(estimates)
}


#' Get Leave-One-Out Effect Sizes
#'
#' Extracts and aggregates effect size data from each leave-one-out iteration.
#'
#' @param outputs A named list of model objects from leave-one-out analysis.
#' @param group A string specifying the grouping variable used in the analysis.
#'
#' @return A data frame of effect sizes with a column indicating the omitted group.
#'
#' @keywords internal
get_effectsizes <- function(outputs, group) {
    effect_sizes <- do.call(rbind, lapply(names(outputs), function(name) {
        res <- mod_results(outputs[[name]], group = group)
        df <- res$data
        df$moderator <- name  
        df
    }))

    row.names(effect_sizes) <- NULL
    return(effect_sizes)
}
