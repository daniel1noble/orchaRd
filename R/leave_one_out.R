#' Leave-One-Out Analysis for Meta-Analytic Models
#'
#' Performs a leave-one-out analysis on a meta-analytic model from
#' the **metafor** package by iteratively removing each level of a grouping
#' variable and refitting the model.
#'
#' @param model A meta-analytic model fitted using **metafor** package functions such as rma.mv().
#' @param group A character string specifying the column name in model$data that contains the 
#'        grouping variable for which levels will be iteratively removed.
#' @param vcalc_args Optional list of arguments for variance-covariance calculation using 
#'        metafor's vcalc function. Must include:
#'        \itemize{
#'          \item \code{vi}: Name of the variance column
#'          \item \code{cluster}: Name of the clustering variable column
#'          \item \code{obs}: Name of the observation ID column
#'          \item \code{rho}: Correlation coefficient between effect sizes
#'        }
#' @param robust_args Optional list of arguments for robust variance estimation using
#'        metafor's robust function. Must include:
#'        \itemize{
#'          \item \code{cluster}: Name of the clustering variable column
#'          \item \code{clubSandwich}: Logical indicating whether to use clubSandwich method (optional)
#'        }
#' @param phylo_args Optional list of arguments for phylogenetic matrix calculation using
#'        ape's vcv function. Must include:
#'        \itemize{
#'          \item \code{tree}: A phylogenetic tree object of class "phylo"
#'          \item \code{species_colname}: Name of the column in model data linked to the tree tips
#'        }
#'
#' @return An object of class "orchard" containing:
#'   \itemize{
#'     \item \code{mod_table}: A data frame with model estimates from each leave-one-out iteration,
#'           with an additional column indicating which group was omitted.
#'     \item \code{data}: A data frame with effect sizes from each iteration, with an additional
#'           column indicating which group was omitted.
#'     \item \code{orig_mod_results}: The results from the original model without any omissions,
#'           as returned by mod_results().
#'   }
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

  # Get the original model results
  orig_mod_results <- mod_results(model, group = group)

  # Immitates the output of mod_results().
  #   - mod_table: In this case, the estimates from each model ran
  #   - data:  Effect sizes and sampling variances from each model
  # Plus one more element:
  #   - orig_mod_results: mod_results() output from the original model
  output <- list(mod_table = estimates,
                 data = effect_sizes,
                 orig_mod_results = orig_mod_results)
  class(output) <- c("orchard", "data.frame")

  return(output)
}


#' Fit Multiple Meta-Analytic Models For Leave-One-Out Analysis
#'
#' Iteratively refits a meta-analytic model, leaving out one level of a specified 
#' grouping variable in each iteration. This internal function handles the actual model
#' refitting process for the leave-one-out analysis.
#'
#' @param model A fitted metafor model object containing a \code{data} element and \code{call} object.
#' @param group A character string specifying the column in \code{model$data} that 
#'        defines the groups to be omitted one at a time.
#' @param vcalc_args Optional list of arguments for the variance-covariance calculation using 
#'        metafor's vcalc function. See \code{leave_one_out()} for details.
#' @param robust_args Optional list of arguments for robust variance estimation using
#'        metafor's robust function. See \code{leave_one_out()} for details.
#' @param phylo_args Optional list of arguments for phylogenetic matrix calculation.
#'        See \code{leave_one_out()} for details.
#'
#' @details The function creates a subset of the original data for each unique value in the
#'          grouping variable, removing that group and refitting the model using the same
#'          specification as the original model. If variance-covariance matrices or 
#'          phylogenetic corrections were used in the original model, these are 
#'          recalculated for each subset. If a model fails to fit, NULL is returned for that group.
#'
#' @return A named list of models, each fitted after omitting one group. Names correspond 
#'         to the omitted group IDs. Any failed model fits will be NULL.
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
    tmp_model_call$data <- model$data[model$data[[group]] != id_left_out, ]

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
#'     \item rho: Correlation coefficient between effect sizes
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
#' valid and refer to existing variables in the model data. Performs checks on
#' the structure and content of the vcalc_args list.
#'
#' @param model_data A data frame containing the variables used in the model.
#' @param vcalc_args A list of arguments for the metafor::vcalc function.
#'
#' @return The validated vcalc_args list if all checks pass.
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

  # Check if the vcalc arguments are present in the model data
  if (is.null(model_data[[vcalc_args$vi]]) || is.null(model_data[[vcalc_args$cluster]]) || is.null(model_data[[vcalc_args$obs]])) {
    stop("One or more of the vcalc arguments are not found in the model data", call. = FALSE)
  }

  return(vcalc_args)
}

#' Validate Robust Variance Estimation Arguments
#'
#' Validates that robust_args contains the required parameters and that they
#' reference valid columns in the data. Ensures that the cluster variable exists
#' in the provided model data.
#'
#' @param model_data A data frame containing the variables used in the model.
#' @param robust_args A list of arguments for the metafor::robust function.
#'
#' @return The validated robust_args list if all checks pass.
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
#' Extracts and combines the model estimates from each leave-one-out iteration into a 
#' single data frame. Applies mod_results() to each model and combines the resulting
#' mod_table values.
#'
#' @param outputs A named list of model objects from leave-one-out analysis where each
#'        name corresponds to the omitted group ID.
#' @param group A character string specifying the grouping variable used in the analysis.
#'
#' @return A data frame of model estimates with an added column 'name' indicating the omitted group.
#'         Contains the same structure as the mod_table from mod_results() for each model iteration.
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
#' Extracts and aggregates effect size data from each leave-one-out iteration into a
#' single data frame. Applies mod_results() to each model and combines the resulting
#' data values.
#'
#' @param outputs A named list of model objects from leave-one-out analysis where each
#'        name corresponds to the omitted group ID.
#' @param group A character string specifying the grouping variable used in the analysis.
#'
#' @return A data frame of effect sizes with a column 'moderator' indicating the omitted group.
#'         Contains the same structure as the data element from mod_results() for each model iteration.
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

#' Prune Phylogenetic Tree For Leave-One-Out Analysis
#' 
#' Removes species from a phylogenetic tree that are not present in the current data subset
#' during leave-one-out analysis.
#'
#' @param tree A phylogenetic tree object of class "phylo".
#' @param species_names A vector of species names that should remain in the pruned tree.
#'        These are the species present in the current data subset.
#' 
#' @return A pruned phylogenetic tree containing only the tips for species in species_names.
#'
#' @author Daniel Noble  - daniel.noble@anu.edu.au
#' @author Facundo Decunta - fdecunta@agro.uba.ar
#'
#' @keywords internal
.prune_tree <- function(tree, species_names) {
  tree_species <- tree$tip.label
  data_species <- unique(species_names)

  species_to_prune <- setdiff(tree_species, data_species)

  if (length(species_to_prune) > 0) {
    tree <- ape::drop.tip(tree, species_to_prune)
  }

  return(tree)
}

#' Create Temporary Phylogenetic Matrix For Leave-One-Out Analysis
#' 
#' Creates a correlation matrix from a phylogenetic tree for the species remaining
#' in the data subset during leave-one-out analysis.
#'
#' @param data A data frame containing the variables specified in phylo_args,
#'        representing the current data subset after removing one group.
#' @param phylo_args A list of arguments for the phylogenetic matrix calculation, including:
#'   \itemize{
#'     \item tree: A phylogenetic tree object of class "phylo"
#'     \item species_colname: Name of the column in the data that contains species names
#'   }
#'
#' @return A phylogenetic correlation matrix for the species in the data subset.
#'
#' @author Facundo Decunta - fdecunta@agro.uba.ar
#'
#' @keywords internal
.create_tmp_phylo_matrix <- function(data, phylo_args) {
  orig_tree <- phylo_args$tree
  species_colname <- phylo_args$species_colname

  # Remove species that are left out
  pruned_tree <- .prune_tree(orig_tree, data[[species_colname]])
  pruned_tree <- ape::compute.brlen(pruned_tree)

  # Compute the phylo matrix 
  tmp_phylo_matrix <- ape::vcv(pruned_tree, corr = TRUE)

  return(tmp_phylo_matrix)
}

#' Validate Phylogenetic Arguments
#'
#' Checks the validity of the arguments provided for phylogenetic matrix calculations.
#' Ensures the tree is a proper phylogenetic object and that the species column name
#' is correctly linked to a random effect in the model.
#'
#' @param model A metafor model object with random effects.
#' @param phylo_args A list containing the arguments for the phylogenetic matrix calculation:
#'   \itemize{
#'     \item tree: A phylogenetic tree object of class "phylo"
#'     \item species_colname: Name of the column in the data corresponding to species names,
#'           which should be a random factor linked to a matrix in the model
#'   }
#'
#' @return The validated phylo_args list if all checks pass.
#'
#' @keywords internal
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
