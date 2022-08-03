# test
rm(list=ls())


install.packages("devtools")
install.packages("tidyverse")
install.packages("metafor")
install.packages("patchwork")
install.packages("R.rsp")

remotes::install_github("daniel1noble/orchaRd", force = TRUE)
remotes::install_github("rvlenth/emmeans", dependencies = TRUE, build_opts = "")

library(orchaRd)
library(metafor)
library(emmeans)
library(tidyverse)

######################
# creating bubble plot
######################

# Issues
# TODO - if interaction combinations are missing - qdrg
# TODO - qdrg does not work with poly


#model$data

#grid <- qdrg(object = model,  at = list("deg_dif" = seq(1,15, length.out = 100)))
#test <-emmeans(grid, specs = "deg_dif", by =  c("deg_dif", "trait.type"))

data(lim)
lim[, "year"] <- as.numeric(lim$year)
lim$vi<- 1/(lim$N - 3)
model<-rma.mv(yi=yi, V=vi, mods= ~Environment*year, random=list(~1|Article,~1|Datapoint), data=lim)
#grid <- qdrg(object = model,  at = list("year" = seq(1970, 2015, length.out = 100)))
#test <-emmeans(grid, specs = "year", by =  c("year", "Environment"))

test <- mod_results(model, mod = "year", group = "Article", data = lim, weights = "prop", by = "Environment")
bubble_plot(test, mod = "year", legend.pos = "top.left", g = T, data = lim)

test2 <- mod_results(model, mod = "year", group = "Article", data = lim, weights = "prop")
bubble_plot(test2, mod = "year", legend.pos = "top.left", g = T, data = lim)


# Data
data(fish)
warm_dat <- fish


model2 <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list(~1 | group_ID, ~1 | es_ID), mods = ~ experimental_design + trait.type+deg_dif*treat_end_days, method = "REML", test = "t", data = warm_dat, control=list(optimizer="optim", optmethod="Nelder-Mead"))


mod2_results <- mod_results(model2, mod = "deg_dif", group = "group_ID", data = warm_dat)
bubble_plot(mod2_results, mod = "deg_dif",condition.nrow = 3, legend.pos = "bottom.left") +
  ylim(-1.5,1.8)

### poly

# read - this https://github.com/rvlenth/emmeans/issues/43
# TODO

data(lim)
lim[, "year"] <- as.numeric(lim$year)
lim$vi<- 1/(lim$N - 3)

lim <- lim[complete.cases(lim), ]

model<-rma.mv(yi=yi, V=vi, mods= ~poly(year, degree = 2) , random=list(~1|Article,~1|Datapoint), data=lim)
summary(model)

model <- lm(yi ~ poly(year, degree = 2) + Environment, data = lim)

# model<-rma.mv(yi=yi, V=vi, mods= ~ 1 + year, random=list(~1|Article,~1|Datapoint), data=lim)
# summary(model)

grid <- qdrg(object = model, data = lim)
test <-emmeans(model, specs = "year", at = list("year" = seq(min(lim$year) , max(lim$year), length.out = 100)))
plot(test)



######################
# testing orchard_plot
######################
# Data
data(fish)
warm_dat <- fish

# The Model
model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list(~1 | group_ID, ~1 | es_ID), mods = ~ experimental_design + trait.type+deg_dif + treat_end_days, method = "REML", test = "t", data = warm_dat,                               control=list(optimizer="optim", optmethod="Nelder-Mead"))

model0 <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list(~1 | group_ID, ~1 | es_ID), method = "REML", test = "t", data = warm_dat,                               control=list(optimizer="optim", optmethod="Nelder-Mead"))


orchard_plot(model, mod = "1", group = "group_ID", data = warm_dat, xlab = "lnRR")
orchard_plot(model, mod = "experimental_design", group = "group_ID", data = warm_dat, xlab = "lnRR", colour = "mod")
#test <- get_data_raw(model, mod = "1", group = "group_ID", data = warm_dat)
#test <- mod_results(model, mod = "1", group = "group_ID", data = warm_dat)

#get_data_raw(model, mod = "experimental_design", group = "group_ID", data = warm_dat)

#orchard_plot(test, mod = "1", group = "group_ID", data = warm_dat, xlab = "lnRR")

# +
#   scale_fill_manual(values="grey") +
#   scale_colour_manual(values="grey") +
#   scale_x_discrete(labels = c('Intercept'))

# This works fine for me,DN
orchard_plot(model0, xlab = "lnRR", trunk.size = 1, branch.size = 2, twig.size = 0.5,
             angle = 45, group = "group_ID", data = warm_dat, legend.pos = "none")

# +
#   scale_fill_manual(values="grey") +
#   scale_colour_manual(values="grey") +
#   scale_x_discrete(labels = c('Intercept'))



####

mod_results2 <- function(model, mod = "1", group, data, weights = "prop", by = NULL, at = NULL, subset = FALSE, ...){

  if(missing(model)){
    stop("Please specify the 'model' argument by providing rma.mv or rma model object. See ?mod_results")
  }

  if(all(class(model) %in% c("rma.mv", "rma")) == FALSE) {stop("Sorry, you need to fit a metafor model of class rma.mv or rma")}

  if(missing(group)){
    stop("Please specify the 'group' argument by providing the name of the grouping variable. See ?mod_results")
  }

  if(missing(data)){
    stop("Please specify the 'data' argument by providing the data used to fit the model. See ?mod_results")
  }

  model$data <- data

  # categorical and continuous
  test <- emmeans::qdrg(object = model)

  if(is.null(formula(model))){
    model <- stats::update(model, "~1")
  }

  if(model$test == "t"){
    df_mod = as.numeric(model$ddf[[1]])
  } else{
    df_mod = 1.0e6 # almost identical to z value
  }

  if(any(names(test@model.info$xlev) == mod)) {
    grid <- emmeans::qdrg(object = model, at = at)
    mm <- emmeans::emmeans(grid, specs = mod, df = df_mod, by = by, weights = weights, ...)

    # getting prediction intervals
    mm_pi <- pred_interval_esmeans(model, mm, mod = mod)


  if(is.null(by)){
    mod_table <- data.frame(name = firstup(as.character(mm_pi[,1])),
                            estimate = mm_pi[,"emmean"],
                            lowerCL = mm_pi[,"lower.CL"],
                            upperCL = mm_pi[,"upper.CL"],
                            lowerPR = mm_pi[,"lower.PI"],
                            upperPR = mm_pi[,"upper.PI"])

  } else{
    mod_table <- data.frame(name = firstup(as.character(mm_pi[,1])),
                            condition = mm_pi[,2], estimate = mm_pi[,"emmean"],
                            lowerCL = mm_pi[,"lower.CL"],
                            upperCL = mm_pi[,"upper.CL"],
                            lowerPR = mm_pi[,"lower.PI"],
                            upperPR = mm_pi[,"upper.PI"])
  }

    # Extract data
    data2 <- get_data_raw(model, mod, group, data, at = at, subset)

    mod_table$name <- factor(mod_table$name,
                             levels = mod_table$name,
                             labels = mod_table$name)

  } else{
    at <- list(mod = seq(min(data[,mod]), max(data[,mod]), length.out = 100))
    names(at) <- mod
    grid <- emmeans::qdrg(object = model, at = at)  # getting 100 points
    mm <- emmeans::emmeans(grid, specs = mod, by = mod, weights = weights, df = df_mod)

    # getting prediction intervals
    mm_pi <- pred_interval_esmeans(model, mm, mod = mod)
    mod_table <- data.frame(moderator = mm_pi[,1],
                            estimate = mm_pi[,"emmean"],
                            lowerCL = mm_pi[,"lower.CL"],
                            upperCL = mm_pi[,"upper.CL"],
                            lowerPR = mm_pi[,"lower.PI"],
                            upperPR = mm_pi[,"upper.PI"])
    # extract data
    data2 <- get_data_raw2(model, mod, group, data)

  }


  output <- list(mod_table = mod_table,
                 data = data2)

  class(output) <- c("orchard", "data.frame")

  return(output)
}


get_data_raw2 <- function(model, mod, group, data){

  if(missing(group)){
    stop("Please specify the 'group' argument by providing the name of the grouping variable. See ?mod_results")
  }

  if(missing(data)){
    stop("Please specify the 'data' argument by providing the data used to fit the model. See ?mod_results")
  }

  # Extract data
  # Check first if missing data exists
  if(length(attr(model$X, "dimnames")[[1]]) > 0){
    # full model delete missing values so need to adjust
    position <- as.numeric(attr(model$X, "dimnames")[[1]])
    data <- data[position, ] }

    # Extract effect sizes
    yi <- model$yi
    vi <- model$vi
    type <- attr(model$yi, "measure")

    # Get moderator
    moderator <- data[,mod] # Could default to base instead of tidy
    #moderator <- firstup(moderator)

  # Extract study grouping variable to calculate the
  stdy <- data[,group] # Could default to base instead of tidy

  data_reorg <- data.frame(yi, vi, moderator, stdy, type)
  row.names(data_reorg) <- 1:nrow(data_reorg)
  return(data_reorg)
}





#########
# Example 1: Overall marginal means for each level of experimental design
  # Two step with marginal means
  overall <- marginal_means(model, mod = "experimental_design", group = "group_ID",
                            data = warm_dat)
  orchard_plot(overall, xlab = "lnRR", trunk.size = 2, branch.size = 2, twig.size = 0.5,
               angle = 45)

  # Directly with model
  orchard_plot(model, xlab = "lnRR", data = warm_dat, mod = "experimental_design", group = "group_ID", trunk.size = 2, branch.size = 2, twig.size = 0.5, angle = 45, marginal = TRUE)

# Example 2: Overall marginal mean across all designs
  # Two step with marginal means
  overall1.1 <- marginal_means(model, group = "group_ID", data = warm_dat)
  orchard_plot(overall1.1, xlab = "lnRR", trunk.size = 2, branch.size = 1.2, twig.size = 2)

  # Directly with model
  orchard_plot(model, group = "group_ID", data = warm_dat, xlab = "lnRR", trunk.size = 2, branch.size = 1.2, twig.size = 2, marginal = TRUE)

# Example 3: Marginalised overall mean for each temp
  # Two step with marginal means
  overall2 <- marginal_means(model, group = "group_ID", mod = "1", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif",  weights = "prop", data = warm_dat)
  orchard_plot(overall2, xlab = "lnRR", group = "group_ID", condition.lab = "Temparature")

  # Directly with model
  orchard_plot(model, xlab = "lnRR", group = "group_ID", condition.lab = "Temparature", mod = "1", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif",  weights = "prop", data = warm_dat, marginal = TRUE)


# Example 4: Marginalised means within each trait category as opposed to design
# Two step with marginal means
across_trait <- marginal_means(model,  mod = "trait.type", group = "group_ID", data = warm_dat)
orchard_plot(across_trait, xlab = "lnRR")

# With model directly
orchard_plot(model, xlab = "lnRR", mod = "trait.type", group = "group_ID", data = warm_dat, marginal = TRUE)


# Example 5: Marginalised means for each trait category by different temperature differences
  # Two step with marginal means
  across_trait_by_degree_diff <- marginal_means(model, mod = "trait.type", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif", group = "group_ID", data = warm_dat)
  orchard_plot(across_trait_by_degree_diff,  xlab = "lnRR",  data = warm_dat)

  # With model directly
  orchard_plot(model, mod = "trait.type", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif", group = "group_ID", data = warm_dat, xlab = "lnRR", marginal = TRUE)



  # Example 6: Marginalised means for each trait category by different temperature differences at treatment end days held at 10
  # Two step with marginal means
across_trait_by_degree_diff_at_treat_end_days10 <- marginal_means(model, data = warm_dat, mod = "trait.type", at = list(deg_dif = c(5, 10, 15),  treat_end_days = 10), group = "group_ID", by = "deg_dif")
orchard_plot(across_trait_by_degree_diff_at_treat_end_days10, xlab = "lnRR")

# With model directly
orchard_plot(model, data = warm_dat, mod = "trait.type", at = list(deg_dif = c(5, 10, 15),  treat_end_days = 10), group = "group_ID", by = "deg_dif", xlab = "lnRR", marginal = TRUE)


# Example 7: Marginalised means for each trait category by different temperature differences averaging only treatment end days 10 & 50
# Two step with marginal means
across_trait_by_degree_diff_at_treat_end_days10And50 <- marginal_means(model, data = warm_dat, at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), mod = "trait.type", by = "deg_dif", group = "group_ID")
orchard_plot(across_trait_by_degree_diff_at_treat_end_days10And50, xlab = "lnRR")

# With model directly
orchard_plot(model, data = warm_dat, at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), mod = "trait.type", by = "deg_dif", group = "group_ID", xlab = "lnRR", marginal = TRUE)

# Example 8: Marginalised means for each trait category by different temperature differences averaging only treatment end days 10 & 50 with equal weights
# Two step with marginal means
across_trait_by_treat_end_days10And50 <- marginal_means(model, data = warm_dat, mod = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), by = "treat_end_days", weights = "equal", group = "group_ID")
orchard_plot(across_trait_by_treat_end_days10And50, xlab = "lnRR")

# With model directly
orchard_plot(model, data = warm_dat, mod = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), by = "treat_end_days", weights = "equal", group = "group_ID", xlab = "lnRR", marginal = TRUE)

# Example 9: Marginalised means for each trait category assuming heteroscedastic error within each
# Model
model_het <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
                             random = list(~1 | group_ID, ~1 + trait.type| es_ID),
                             mods = ~ trait.type + deg_dif, method = "REML",
                             test = "t", rho = 0, struc = "HCS", data = warm_dat,
                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

# Two step process
HetModel <- marginal_means(model_het,
                           group = "group_ID",
                           mod = "trait.type",
                           at = list(deg_dif = c(5, 10, 15)),
                           by = "deg_dif", weights = "prop", data = warm_dat)
orchard_plot(HetModel, xlab = "lnRR")

# With model directly
orchard_plot(model_het, group = "group_ID",
             mod = "trait.type",
             at = list(deg_dif = c(5, 10, 15)),
             by = "deg_dif", weights = "prop", data = warm_dat,
             xlab = "lnRR", marginal = TRUE)

test
