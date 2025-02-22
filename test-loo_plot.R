library(devtools)
library(metafor)
library(tidyverse)


load_all()

# Let's see if by changhin the model data, the VCV should be computed again or not

VCV <- vcalc(lnrr_vi,
             cluster = fish$paper_ID,
             obs     = fish$es_ID,
             data    = fish,
             rho     = 0.5)

res <- rma.mv(lnrr, VCV,
              random = ~ 1 | paper_ID,
              data = fish)

res2 <- rma.mv(lnrr, lnrr_vi,
              random = ~ 1 | paper_ID,
              data = fish)



res

res2

# Check how to output with the data used so it can be used for orchard plot?

outputs <- run_leave1out(res2, group = "paper_ID")

estimates <- get_estimates(outputs, group = "paper_ID")
effect_sizes <- get_effectsizes(outputs, group = "paper_ID")

library(ggplot2)

ggplot() +
  geom_point(data = effect_sizes, aes(x = left_out,
                                      y = yi),
#                                      size = 1/sqrt(vi)),
             alpha = 0.05) +
  geom_pointrange(data = estimates, aes(x = left_out,
                                        y = estimate,
                                        ymin = lowerCL,
                                        ymax = upperCL),
                  color = "darkorange") +
  coord_flip() +
  ylim(-0.2, 0.2)  +
  theme_minimal()



# --- eklof dataset

data(eklof)

# Calculate the effect size
eklof <- escalc(measure = "ROM", n1i = N_control, sd1i = SD_control, m1i = mean_control,
    n2i = N_treatment, sd2i = SD_treatment, m2i = mean_treatment, var.names = c("lnRR",
        "vlnRR"), data = eklof)

# Add the observation level factor
eklof$Datapoint <- as.factor(seq(1, dim(eklof)[1], 1))

# Also, we can get the sample size, which we can use for weighting if we would
# like
eklof$N <- rowSums(eklof[, c("N_control", "N_treatment")])

# fit a meta-regression with the intercept (contrast)
eklof_MR0 <- rma.mv(yi = lnRR, V = vlnRR, mods = ~Grazer.type, random = list(~1 |
    ExptID, ~1 | Datapoint), data = eklof)

summary(eklof_MR0)


results <- mod_results(eklof_MR0, group = "ExptID")

str(results)

# What about mods?

eklof_MR0 <- rma.mv(yi = lnRR,
                    V = vlnRR,
                    mods = ~Grazer.type,
                    random = list(~1 | ExptID,
                                  ~1 | Datapoint),
                    data = eklof)

summary(eklof_MR0)

results <- mod_results(eklof_MR0, mod = "Grazer.type", group = "ExptID")

results

results$data

str(results)

# Imitate the outpute of mod_results but with leave1out

outputs_eklof <- run_leave1out(eklof_MR0, group = "ExptID")
estimates_eklof <- get_estimates(outputs_eklof, group = "ExptID")
effect_sizes_eklof <- get_effectsizes(outputs_eklof, group = "ExptID")

estimates_eklof$name <- estimates_eklof$left_out
effect_sizes_eklof$moderator <- effect_sizes_eklof$left_out

# This immitates the output of mod_results.
# Here the estimates are the estimates of the model with the left out study
# while the data are the effect sizes of each model. 
leave1out_eklof <- list(mod_table = estimates,
                        data      = effect_sizes)
class(leave1out_eklof) <- c("orchard", "data.frame")

str(leave1out_eklof)

loo <- leave_one_out(eklof_MR0, group = "First.author")

loo


loo$data

orchard_plot(loo, xlab = "bla", alpha = 0.2)

