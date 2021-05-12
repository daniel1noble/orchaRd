# test

remotes::install_github("rvlenth/emmeans")
library(orchaRd)
library(metafor)
library(emmeans)


#' @title marginalised_means
#' @description Function to to get marginalised means from met-regression models with single or multiple moderator variables that are both continuous or categorical.
#' @param model rma.mv object
#' @param data data frame used to fit rma.mv model
#' @param pred predictor variable of interest that one wants marginalised means for.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @example \dontrun{
warm_dat <- fish
model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list(~1 | group_ID, ~1 | es_ID), mods = ~ experimental_design + trait.type + deg_dif + treat_end_days, method = "REML", test = "t", data = warm_dat,                               control=list(optimizer="optim", optmethod="Nelder-Mead"))
  overall <- marginalised_means(model, data = warm_dat)
across_trait <- marginalised_means(model, data = warm_dat, pred = "trait.type")
across_trait_by_degree_diff <- marginalised_means(model, data = warm_dat, pred = "trait.type", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif")
across_trait_by_degree_diff_at_treat_end_days10 <- marginalised_means(model, data = warm_dat, pred = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = 10), by = "deg_dif")
across_trait_by_degree_diff_at_treat_end_days10And50 <- marginalised_means(model, data = warm_dat, pred = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), by = "deg_dif")
across_trait_by_treat_end_days10And50 <- marginalised_means(model, data = warm_dat, pred = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), by = "treat_end_days")
across_trait_by_treat_end_days10And50_ordinaryMM <- marginalised_means(model, data = warm_dat, pred = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), by = "treat_end_days", weights = "prop")

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

# current get_data

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


# TODO
# - we can do - get_data2
# - we need to also change -



