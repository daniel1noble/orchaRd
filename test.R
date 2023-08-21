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
# emprep for metafor
######################

data(fish)
warm_dat <- fish
model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
random = list(~1 | group_ID, ~1 | es_ID), mods = ~ experimental_design + trait.type + deg_dif + treat_end_days,
method = "REML", test = "t", control=list(optimizer="optim", optmethod="Nelder-Mead"), data = warm_dat)
  

mod_prep1 <- emmprep(model)
emmeans(mod_prep1, specs="1", type="response", weights="proportional")
mod_results(model, group = "group_ID")

emmeans(mod_prep1, specs="trait.type", type="response", weights="proportional")
str(mod_results(model, group = "group_ID", mod = "trait.type"))

model_het <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list(~1 | group_ID, ~1 + trait.type| es_ID), mods = ~ trait.type + deg_dif, method = "REML", test = "t", dfs = "contain", rho = 0, struc = "HCS", control=list(optimizer="optim", optmethod="Nelder-Mead"), data = warm_dat)

mod_prep2 <- emmprep(model_het, at = list(deg_dif = c(5, 10, 15)), by = "deg_dif")
grid_qdrg <- emmeans::qdrg(formula = stats::formula(model_het), at = list(deg_dif = c(5, 10, 15)), by = "deg_dif", data = model_het$data, coef = model_het$b, vcov = stats::vcov(model_het), df = model_het$k-1) 

emmeans(mod_prep2, specs = "trait.type", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif", type="response", weights="proportional")
emmeans(grid_qdrg, specs = "trait.type", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif", type="response", weights="proportional")
mod_results(model_het, group = "group_ID", mod = "trait.type", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif", weights = "prop")

grid <- metafor::emmprep(model_het, at = list(deg_dif = c(5, 10, 15)), by = "deg_dif")
      mm <- emmeans::emmeans(grid, specs =  "trait.type", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif")

    # getting prediction intervals
    mm_pi <- pred_interval_esmeans(model_het, mm, mod =  "trait.type")


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
bubble_plot(test, mod = "year", legend.pos = "top.left", group = "Article", g = T, data = lim)

test2 <- mod_results(model, mod = "year", group = "Article", data = lim, weights = "prop")
bubble_plot(test2, mod = "year", legend.pos = "top.left", group = "Article", g = T, data = lim)

######################
## Bubble plot and hetero
data(lim)
lim[, "year"] <- as.numeric(lim$year)
lim$vi<- 1/(lim$N - 3)
model<-metafor::rma.mv(yi=yi, V=vi, mods= ~Amniotes*year,
                       random=list(~1|Article,~1+Amniotes|Datapoint), rho = 0, str="HCS", data=na.omit(lim))

model<-metafor::rma.mv(yi~Amniotes*year, V=vi,
                       random=list(~1|Article,~1+Amniotes|Datapoint), rho = 0, str="HCS", data=na.omit(lim))

lim_bubble <- orchaRd::mod_results(model, mod = "year", group = "Article",
                                   data = lim, weights = "prop", by = "Amniotes")

orchaRd::bubble_plot(lim_bubble, data = lim, group = "Article", mod = "year", xlab = "Year", legend.pos = "top.left")

## Testing if categorical interactions are a problem. Note here that we get Warning message:
# Redundant predictors dropped from the model. And then this fails. So, I think emmeans just can't have a situation where levels of one variable are dropped because they don't fit in the model. Should be easy to fix this by creating a interaction between two variables and fitting the model differently

modelLim1<-metafor::rma.mv(yi~RU*Order, V=vi,
                       random=list(~1|Article,~1+Amniotes|Datapoint), rho = 0, str="HCS", data=na.omit(lim))

orchard_plot(modelLim1, mod = "year", group = "Article", data = na.omit(lim), xlab = "Zr") # FAILS

lim$new_fac <- with(lim, as.character(interaction(Order, RU))) # Try as character to make sure that when something is dropped the level doesn't remain if a factor.
modelLim2<-metafor::rma.mv(yi~new_fac, V=vi,
                           random=list(~1|Article,~1+new_fac|Datapoint), rho = 0, str="HCS", data=na.omit(lim))

orchard_plot(modelLim2, mod = "new_fac", group = "Article", data = na.omit(lim), xlab = "Zr", cb = FALSE) # WORKS


# Data
data(fish)
warm_dat <- fish


model2 <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list(~1 | group_ID, ~1 | es_ID), mods = ~ experimental_design + trait.type+deg_dif*treat_end_days, method = "REML", test = "t", data = warm_dat, control=list(optimizer="optim", optmethod="Nelder-Mead"))


mod2_results <- mod_results(model2, mod = "deg_dif", group = "group_ID", data = warm_dat)
bubble_plot(mod2_results, mod = "deg_dif",condition.nrow = 3, group = "group_ID", legend.pos = "bottom.left", data = warm_dat) + ylim(-1.5,1.8)

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

  # Directly with model
  orchard_plot(model, xlab = "lnRR", data = warm_dat, mod = "experimental_design", group = "group_ID", trunk.size = 2, branch.size = 2, twig.size = 0.5, angle = 45)

# Example 2: Overall marginal mean across all designs
  # Two step with marginal means

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

filtered <- warm_dat %>% filter(trait.type %in% c("morphology", "physiology"))
model_het_gamma <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
                             random = list(~1 + trait.type | group_ID, ~1 + trait.type| es_ID),
                             mods = ~ trait.type + deg_dif, method = "REML",
                             test = "t", rho = 0, phi = 0, struc = "HCS", data = filtered,
                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

model_nohet <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
                             random = list(~1 | group_ID, ~1 | es_ID),
                             mods = ~ trait.type + deg_dif, method = "REML",
                             test = "t", data = warm_dat,
                             control=list(optimizer="optim", optmethod="Nelder-Mead"))
warm_dat %>% group_by(trait.type) %>% summarise(n = length(unique(group_ID)))

### CHECK R2 and I2 with het models.
orchaRd::r2_ml(model_nohet, data = warm_dat)
orchaRd::r2_ml(model_het, data = warm_dat) # Works, but not correct
orchaRd::r2_ml(model_het, data = warm_dat, boot = 100) # Works, but not correct

orchaRd::i2_ml(model_nohet, data = warm_dat)
orchaRd::i2_ml(model_het, data = warm_dat) # Works, but not correct
orchaRd::i2_ml(model_het, data = warm_dat, boot = 100) # Works, but not correct


R2_calc <- function(model){
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  model = model_het
  if(any(model$tau2 > 0)) {
    # extract number of levels
      g_sigma <- model$s.nlevels
      g_tau   <- model$g.nlevels[2]
    g_gamma   <- model$h.nlevels[2]

    comp_group <- c(g_sigma, g_tau, g_gamma) #use this to get max g level

   # Sample sizes for each level
      k_tau <- model$g.levels.k # extract sample size for each level
    k_gamma <- model$h.levels.k # extract sample size for each level

    # Extract variances for each level.
        tau2 <- model$tau2
      gamma2 <- model$gamma2

    # Calculated the weighted variance, weighted on sample size
      tau_var <- orchaRd::weighted_var(tau2, weights = k_tau)
    gamma_var <- orchaRd::weighted_var(gamma2, weights = k_gamma)

    # Composite variance
    vars <- c(model$sigma2, tau_var, gamma_var)

    # fixed effect variance
    fix <- stats::var(as.numeric(as.vector(model$b) %*% t(as.matrix(model$X))))

    # marginal
    tot_res_var <- sum(vars)

    R2m <- fix / (fix + tot_res_var)

    # conditional
    R2c <- (fix + (tot_res_var - tau_var)) /
      (fix + tot_res_var)

  } else{

    # fixed effect variance
    fix <- stats::var(as.numeric(as.vector(model$b) %*% t(as.matrix(model$X))))

    # marginal
    R2m <- fix / (fix + sum(model$sigma2))

    # conditional
    R2c <- (fix + sum(model$sigma2) - model$sigma2[which(model$s.nlevels.f == max(model$s.nlevels.f))]) /
      (fix + sum(model$sigma2))
  }

  R2s <- c(R2_marginal = R2m, R2_conditional = R2c)
  return(R2s)
}




###
g.levels.k <- model_het$g.levels.k
      tau2 <- model_het$tau2

orchaRd::weighted_var(tau2, weights = g.levels.k)


metafor::formula.rma(model_het, type = "mods")
model_het$random
##################################################
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
