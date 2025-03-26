#' Leave-One-Out Analysis for Meta-Analytic Models
#'
#' Performs a leave-one-out analysis on a meta-analytic model from
#' the **metafor** package by iteratively removing each level of a grouping
#' variable and refitting the model.
#'
#' @param model A meta-analytic model fitted using **metafor**.
#' @param group A factor or categorical variable specifying the leave-one-out groups.
#' @param vcalc_args Optional list of arguments for the variance-covariance calculation using 
#'   metafor's vcalc function. Must include 'vi', 'cluster', 'obs', and 'rho' elements.
#' @param robust_args Optional list of arguments for robust variance estimation using
#'   metafor's robust function. Must include a 'cluster' element.
#'
#' @return Same as `mod_results()`, but with the estimates from each model ran
#'   in the leave-one-out analysis, and the effect sizes from each model.
#'
#' @author Facundo Decunta - fdecunta@agro.uba.ar
#'
#' @examples
#' \dontrun{
#' res <- metafor::rma.mv(lnrr, lnrr_vi, random = ~ 1 | paper_ID, data = fish)
#' loo_results <- leave_one_out(res, group = "paper_ID")
#' 
#' # With variance-covariance calculation
#' loo_results <- leave_one_out(res, group = "paper_ID", 
#'                              vcalc_args = list(vi = "lnrr_vi", 
#'                                               cluster = "paper_ID", 
#'                                               obs = "es_ID", 
#'                                               rho = 0.5))
#'                                               
#' # With robust variance estimation
#' loo_results <- leave_one_out(res, group = "paper_ID", 
#'                              robust_args = list(cluster = "paper_ID"))
#' }
#'
#' @export

leave_one_out <- function(model, group, vcalc_args = NULL, robust_args = NULL, phylo_args = NULL) {
  # Check model is a metafor object
  .is_model_valid(model)
  # Check if group is in model data
  .is_group_valid(model$data, group)

  # Check if we have at least 2 groups
  if (length(unique(model$data[[group]])) < 2) {
    stop("Need at least 2 groups for leave-one-out analysis", call. = FALSE)
  }

  if (!is.null(vcalc_args)) {
    .validate_vcalc_args(model$data, vcalc_args)
  } 

  if (!is.null(robust_args)) {
    .validate_robust_args(model$data, robust_args)
  }

  if (!is.null(phylo_args)) {
    .validate_phylo_args(model, phylo_args) 
  }

  # Run leave-one-out analysis
  models_outputs <- .run_leave1out(model, group, vcalc_args, robust_args, phylo_args)
  estimates      <- .get_estimates(models_outputs, group)
  effect_sizes   <- .get_effectsizes(models_outputs, group)

  # Immitates the output of mod_results().
  #   - mod_table: In this case, the estimates from each model ran
  #   - data:  Effect sizes and sampling variances from each model
  output <- list(mod_table = estimates, data = effect_sizes)
  class(output) <- c("orchard", "data.frame")

  return(output)
}


#' Fit Multiple Meta-Analytic Models For Leave-One-Out Analysis
#'
#' Iteratively refits a meta-analytic model, leaving out one level of a specified 
#' grouping variable in each iteration.
#'
#' @param model A fitted model object containing a \code{data} element, which must be 
#'   a data frame with all model variables.
#' @param group A character string specifying the column in \code{model$data} that 
#'   defines the groups to be omitted one at a time.
#' @param vcalc_args Optional list of arguments for the variance-covariance calculation using 
#'   metafor's vcalc function. Must include 'vi', 'cluster', 'obs', and 'rho' elements.
#' @param robust_args Optional list of arguments for robust variance estimation using
#'   metafor's robust function. Must include a 'cluster' element.
#'
#' @details The function removes each unique group from \code{model$data} one at a time, 
#'   refitting the model using \code{update()}. If an update fails, a warning is issued, 
#'   and \code{NULL} is returned for that group.
#'
#' @return A named list of models, each fitted after omitting one group. Names correspond 
#'   to the omitted group IDs.
#'
#' @author Facundo Decunta - fdecunta@agro.uba.ar
#'
#' @keywords internal

.run_leave1out <- function(model, group, vcalc_args = NULL, robust_args = NULL, phylo_args = NULL) {
  group_ids <- unique(model$data[[group]])

  models_outputs <- lapply(group_ids, function(id_left_out) {
    # Create a new call to fit the model. Modify the data to leave out the group
    # and change de VCV and phylo matrix if needed. Then evaluate the new call.

    tmp_model_call <- model$call
    tmp_model_call$data <- subset(model$data, model$data[[group]] != id_left_out)

    # If vcalc_args are provided, create a temporary VCV matrix
    if (!is.null(vcalc_args)) {
      tmp_model_call$V <- .create_tmp_vcv(tmp_model_call$data, vcalc_args)
    }

    # If the model uses phylogenetic matrix, recalculate it using the original tree
    # The model object contains the correlation matrices in 'R'. This is a list
    # where the names are the random effects and the elements are the correlation matrices.
    # So, first create the new matrix, then use it as the matrix linked to phylo_args$species_colname.
    if (!is.null(phylo_args)) {
      tmp_phylo_matrix <- .create_tmp_phylo_matrix(tmp_model_call$data, phylo_args) 
      tmp_model_call$R[[phylo_args$species_colname]] <- tmp_phylo_matrix
    }

    # Evaluate the new call. If something happens, return NULL.
    # In some cases the fixed or random effects are not represented
    # when one group is left out and the model fails to fit.
    tmp_res <- tryCatch({
      eval(tmp_model_call)
    }, error = function(e) {
      warning(sprintf("Error fitting model when leaving out '%s': %s", 
                      id_left_out, e$message))
      return(NULL)
    })

    if(!is.null(robust_args)) {
      # TODO: This way of calling robust must be simpler

      # cluster_var has to be a vector, not a string. robust_args$cluster is a string.
      cluster_var <- tmp_model_call$data[[robust_args$cluster]]

      if(!is.null(robust_args$clubSandwich)) {
        clubSandwich_arg <- robust_args$clubSandwich
      } else {
        clubSandwich_arg <- FALSE
      }

      tmp_res <- metafor::robust(tmp_res, cluster = cluster_var, clubSandwich = clubSandwich_arg)
    }

    # Return the model output so it is saved in 'models_outputs' list
    tmp_res
  })

  names(models_outputs) <- group_ids

  return(models_outputs)
}

#' Create Temporary Variance-Covariance Matrix
#'
#' Creates a variance-covariance matrix for a subset of data using metafor's vcalc function.
#'
#' @param data A data frame containing the variables specified in vcalc_args.
#' @param vcalc_args A list of arguments for metafor::vcalc function, including:
#'   \itemize{
#'     \item vi: Name of the variance column
#'     \item cluster: Name of the clustering variable column
#'     \item obs: Name of the observation ID column
#'     \item rho: Correlation coefficient
#'   }
#'
#' @return A variance-covariance matrix for use in meta-analytic models.
#'
#' @author Facundo Decunta - fdecunta@agro.uba.ar
#'
#' @keywords internal
.create_tmp_vcv <- function(data, vcalc_args) {
  tryCatch({
    metafor::vcalc(vi      = data[[vcalc_args$vi]],
                   cluster = data[[vcalc_args$cluster]],
                   obs     = data[[vcalc_args$obs]],
                   data    = data,
                   rho     = vcalc_args$rho)
  }, error = function(e) {
    stop(sprintf("Error creating VCV: %s", e$message))
  })
}


#' Validate Variance-Covariance Calculation Arguments
#'
#' Ensures that the arguments provided for variance-covariance calculation are
#' valid and refer to existing variables in the model data.
#'
#' @param model_data A data frame containing the variables used in the model.
#' @param vcalc_args A list of arguments for the metafor::vcalc function.
#'
#' @return The validated vcalc_args list.
#'
#' @author Facundo Decunta - fdecunta@agro.uba.ar
#'
#' @keywords internal
.validate_vcalc_args <- function(model_data, vcalc_args) {
  if (!is.list(vcalc_args)) {
    stop("vcalc must be a list with the arguments for the 'vcalc' function: e.g., vcalc_args = list(vi = 'lnrr_vi', cluster = 'paper_ID', obs = 'es_ID', rho = 0.5)",
         call. = FALSE)
  }

  # Check if required arguments for vcalc are present
  if (!all(c("vi", "cluster", "obs", "rho") %in% names(vcalc_args))) {
    stop("vcalc_args must contain at least the following elements: 'vi', 'cluster', 'obs', 'rho'", call. = FALSE)
  }

  # TODO: This is not a compelte check, but it's a start
  # Check if the vcalc arguments are present in the model data
  if (is.null(model_data[[vcalc_args$vi]]) || is.null(model_data[[vcalc_args$cluster]]) || is.null(model_data[[vcalc_args$obs]])) {
    stop("One or more of the vcalc arguments are not found in the model data", call. = FALSE)
  }

  return(vcalc_args)
}

#' Validate robust_args
#'
#' Validates that robust_args contains the required parameters and that they
#' reference valid columns in the data.
#'
#' @param model_data A data frame containing the variables used in the model.
#' @param robust_args A list of arguments for the metafor::robust function.
#'
#' @return The validated robust_args list.
#'
#' @author Facundo Decunta - fdecunta@agro.uba.ar
#'
#' @keywords internal
.validate_robust_args <- function(model_data, robust_args) {
  if (!is.list(robust_args)) {
    stop("robust_args must be a list with the arguments for the 'robust' function: e.g., robust_args = list(cluster = 'paper_ID')",
         call. = FALSE)
  }

  if (!("cluster" %in% names(robust_args))) {
    stop("robust_args must contain at least the following elements: 'cluster'", call. = FALSE)
  }

  # Check if the cluster variable is present in the model data
  if (is.null(model_data[[robust_args$cluster]])) {
    stop("The cluster variable specified in robust_args is not found in the model data", call. = FALSE)
  }

  return(robust_args)
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
#' @author Facundo Decunta - fdecunta@agro.uba.ar
#'
#' @keywords internal
.get_estimates <- function(outputs, group) {
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
#' @author Facundo Decunta - fdecunta@agro.uba.ar
#'
#' @keywords internal
.get_effectsizes <- function(outputs, group) {
    effect_sizes <- do.call(rbind, lapply(names(outputs), function(name) {
        res <- mod_results(outputs[[name]], group = group)
        df <- res$data
        df$moderator <- name  
        df
    }))

    row.names(effect_sizes) <- NULL
    return(effect_sizes)
}

#' Prune Tree For Leave-One-Out
#' 
#' @param tree A phylogenetic tree object.
#' @param species_names A vector of species names.
#' 
#' @author Daniel Noble  - daniel.noble@anu.edu.au
#' @author Facundo Decunta - fdecunta@agro.uba.ar
#'
.prune_tree <- function(tree, species_names) {
  tree_species <- tree$tip.label
  data_species <- unique(species_names)

  species_to_prune <- setdiff(tree_species, data_species)

  if (length(species_to_prune) > 0) {
    tree <- ape::drop.tip(tree, species_to_prune)
  }

  return(tree)
}

#' Create Temporary Phylogenetic Matrix For Leave-One-Out
#' 
#' @param data A data frame containing the variables specified in phylo_args.
#' @param phylo_args A list of arguments for the phylogenetic matrix calculation, including:
#'  \itemize{
#'    \item tree: A phylogenetic tree object
#'    \item species_colname: Name of the column in the data that is linked to the phylo matrix (species names)
#'  }
#'
#' @return A phylogenetic matrix for use in meta-analytic models.
#'
#' @author Facundo Decunta -
#' @keywords internal
.create_tmp_phylo_matrix <- function(data, phylo_args) {
  orig_tree <- phylo_args$tree
  species_colname <- phylo_args$species_colname

  # Remove species that are left out
  pruned_tree <- .prune_tree(orig_tree, data[[species_colname]])
  pruned_tree <- ape::compute.brlen(pruned_tree)

  # Compute the phylo matrix 
  tmp_phylo_matrix <- ape::vcv.phylo(pruned_tree, cor = TRUE)

  return(tmp_phylo_matrix)
}

#' Validate Phylogenetic Arguments
# phylo_args must be a list with the phylogenetic tree and
# the name of the column in the data tha is linked to the phylo matrix (species names).
.validate_phylo_args <- function(model, phylo_args) {
  if (!is.list(phylo_args)) {
    stop("phylo_args must be a list with the arguments for the phylogenetic matrix calculation: e.g., phylo_args = list(tree = tree, species_colname = 'species')",
         call. = FALSE)
  }

  if (!all(c("tree", "species_colname") %in% names(phylo_args))) {
    stop("phylo_args must contain at least the following elements: 'tree', 'species_colname'", call. = FALSE)
  }

  if (class(phylo_args$tree) != "phylo") {
    stop("The 'tree' argument in phylo_args must be a phylogenetic tree object", call. = FALSE)
  }

  # Check if the species_colname is the name of a random factor linked to a matrix in the model
  linked_random_factors <- names(model$Rfix[model$Rfix == TRUE])
  if (!(phylo_args$species_colname %in% linked_random_factors)) {
    stop("The 'species_colname' argument in phylo_args must be the name of the random factor linked to the phylo matrix in the model",
         call. = FALSE)
  }

  return(phylo_args)
}


