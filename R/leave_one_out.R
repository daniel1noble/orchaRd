# The function rerun the model leaving out one study at a time and returns a list with the results.

# Ideas:
#   - What happens if the rma model used `subset`?
#   - It can use the `slab` term from metafor to set the paper names. 
#   - It must handle VCV and phylogenetic models!:
#     - If the user opts for VCV, there should be an option to add the `obs` argument used by metafor::vcalc
#   - Phylogenetic models are a problem. Don't know what to do there.

run_leave1out <- function(model, group) {
  # Rerun the model leaving out one element of `group` at a time.
  # Returns a list with the results with the element left out in each
  # run as names

  group_ids <- unique(model$data[[group]])
  tmp_model <- model

  models_outputs <- lapply(group_ids, function(id_left_out) {

    new_data <- subset(model$data, model$data[[group]] != id_left_out)

    tmp_res <- tryCatch({
      update(tmp_model, data = new_data)
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


models_to_dataframe <- function(models_outputs, coef_name) {
  # Extract from each model the `beta` with its confidence interval.
  # Returns a dataframe with columns: left_out, beta,  ci_lb, ci_up.

  outputs_dataframe <- data.frame(
    left_out = names(models_outputs),
    b        = vapply(models_outputs, function(x) { x$b[1]     }, numeric(1)),
    ci_ub    = vapply(models_outputs, function(x) { x$ci.ub[1] }, numeric(1)),
    ci_lb    = vapply(models_outputs, function(x) { x$ci.lb[1] }, numeric(1))
  )

  rownames(outputs_dataframe) <- NULL   # Remove rownames
  outputs_dataframe
}



leave_one_out <- function(model,
                          group,
                          coef_name) {

  # Rerun the model leaving-one-out
  models_outputs <- run_leave1out(model, group)

  models_to_dataframe(models_outputs)
}
