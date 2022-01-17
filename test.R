# test
install.packages("devtools")
install.packages("tidyverse")
install.packages("metafor")
install.packages("patchwork")
install.packages("R.rsp")

devtools::install_github("daniel1noble/orchaRd", force = TRUE, build_vignettes = TRUE)
remotes::install_github("rvlenth/emmeans", dependencies = TRUE, build_opts = "") 

library(orchaRd)
library(metafor)
library(emmeans)
library(tidyverse)

data(fish)
warm_dat <- fish

# model
model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list(~1 | group_ID, ~1 | es_ID), mods = ~ experimental_design + trait.type + deg_dif + treat_end_days, method = "REML", test = "t", data = warm_dat,                               control=list(optimizer="optim", optmethod="Nelder-Mead"))

# marginal overall
overall <- marginal_means(model, mod = "experimental_design", group = "group_ID")
orchard_plot(overall, xlab = "lnRR", trunk.size = 2, branch.size = 2, twig.size = 0.5, angle = 45)
overall1.1 <- marginal_means(model, group = "group_ID")
orchard_plot(overall1.1, xlab = "lnRR", trunk.size = 2, branch.size = 1.2, twig.size = 2)
overall2 <- marginal_means(model, data = warm_dat, mod = "1", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif",  weights = "prop")
orchard_plot(overall2, xlab = "lnRR", condition.lab = "Temparature")

# marginalized stuff
across_trait <- marginal_means(model, data = warm_dat, mod = "trait.type")
orchard_plot(across_trait, xlab = "lnRR")

across_trait_by_degree_diff <- marginal_means(model, data = warm_dat, mod = "trait.type", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif")
orchard_plot(across_trait_by_degree_diff, xlab = "lnRR")

across_trait_by_degree_diff_at_treat_end_days10 <- marginal_means(model, data = warm_dat, mod = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = 10), by = "deg_dif")
orchard_plot(across_trait_by_degree_diff_at_treat_end_days10, xlab = "lnRR")

across_trait_by_degree_diff_at_treat_end_days10And50 <- marginal_means(model, data = warm_dat, mod = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), by = "deg_dif")
orchard_plot(across_trait_by_degree_diff_at_treat_end_days10And50, xlab = "lnRR")

across_trait_by_treat_end_days10And50 <- marginal_means(model, data = warm_dat, mod = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), by = "treat_end_days", weights = "equal")
orchard_plot(across_trait_by_treat_end_days10And50, xlab = "lnRR")
#
across_trait_by_treat_end_days10And50_ordinaryMM <- marginal_means(model, data = warm_dat, mod = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), by = "treat_end_days", weights = "prop")
orchard_plot(across_trait_by_treat_end_days10And50_ordinaryMM, xlab = "lnRR")


# current get_data

data <-across_trait_by_degree_diff_at_treat_end_days10And50$data
mod_table <- across_trait_by_degree_diff_at_treat_end_days10And50$mod_table
label <- "lnRR (effect size)"
angle <- 90
alpha <- 0.5
scale <- 1/sqrt(data$vi)
legend <- "Precision (1/SE)"
condition_no <- 3
group_no <- 4
condition.lab <- "Tempature "
mod_table$K <- as.vector(by(data, data[,"moderator"], function(x) length(x[,"yi"])))

plot <- ggplot2::ggplot() +
  # pieces of fruit (bee-swarm and bubbles)
  ggbeeswarm::geom_quasirandom(data = data, aes(y = yi, x = moderator, size = scale, colour = moderator), alpha=alpha) +

  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "black", alpha = alpha) +
  # creating dots + CI and PI
  ggplot2::geom_linerange(data = mod_table, aes(x = name, ymin = lowerCL, ymax = upperCL), size = 1.2, position = position_dodge2(width = 0.3)) +
  ggplot2::geom_pointrange(data = mod_table, aes(y = estimate, x = name, ymin = lowerPR, ymax = upperPR,  shape = as.factor(condition), fill = name), size = 0.5, position = position_dodge2(width = 0.3)) +
  scale_shape_manual(values =  20 + (1:condition_no)) +
  #ggplot2::geom_errorbar(aes(ymin = lowerPR, ymax = upperPR), position = position_dodge2(), show.legend = FALSE, size = 0.5, alpha = 0.6, width = 0) +
  coord_flip() +
  # 95 %CI: branches
  #ggplot2::geom_errorbarh(aes(xmin = lowerCL, xmax = upperCL),  height = 0, show.legend = FALSE, size = 1.2, position = position_dodge2()) +
  # putting labels
  #ggplot2::annotate('text', x = (max(data$yi) + (max(data$yi)*0.10)), y = (seq(1, group_no, 1)+0.3),
  #                 label= paste("italic(k)==", mod_table$K), parse = TRUE, hjust = "right", size = 3.5) +
  ggplot2::theme_bw() +
  ggplot2::guides(fill = "none", colour = "none") +
  ggplot2::theme(legend.position= c(0, 1), legend.justification = c(0, 1)) +
  ggplot2::theme(legend.title = element_text(size = 9)) +
  ggplot2::theme(legend.direction="horizontal") +
  ggplot2::theme(legend.background = element_blank()) +
  ggplot2::labs(y = label, x = "", size = legend) +
  ggplot2::labs(shape = condition.lab) +
  ggplot2::theme(axis.text.y = element_text(size = 10, colour ="black",
                                            hjust = 0.5,
                                            angle = 0))

plot <- plot +
  ggplot2::annotate('text', y = (max(data$yi) + (max(data$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
                    label= paste("italic(k)==", mod_table$K[1:group_no]), parse = TRUE, hjust = "right", size = 3.5)





# current get_data

data <-overall2$data
mod_table <- overall2$mod_table
label <- "lnRR (effect size)"
angle <- 90
alpha <- 0.5
scale <- 1/sqrt(data$vi)
legend <- "Precision (1/SE)"
condition_no <- 3
group_no <- 4
condition.lab <- "Temparature"
mod_table$K <- as.vector(by(data, data[,"moderator"], function(x) length(x[,"yi"])))

plot <- ggplot2::ggplot() +
  # pieces of fruit (bee-swarm and bubbles)
  ggbeeswarm::geom_quasirandom(data = data, aes(y = yi, x = moderator, size = scale, colour = moderator), alpha=alpha) +

  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "black", alpha = alpha) +
  # creating dots + CI and PI
  ggplot2::geom_linerange(data = mod_table, aes(x = name, ymin = lowerCL, ymax = upperCL), size = 1.2, position = position_dodge2(width = 0.3)) +
  ggplot2::geom_pointrange(data = mod_table, aes(y = estimate, x = name, ymin = lowerPR, ymax = upperPR,  shape = as.factor(condition), fill = name), size = 0.5, position = position_dodge2(width = 0.3)) +
  scale_shape_manual(values =  20 + (1:condition_no)) +
  #ggplot2::geom_errorbar(aes(ymin = lowerPR, ymax = upperPR), position = position_dodge2(), show.legend = FALSE, size = 0.5, alpha = 0.6, width = 0) +
  coord_flip() +
  # 95 %CI: branches
  #ggplot2::geom_errorbarh(aes(xmin = lowerCL, xmax = upperCL),  height = 0, show.legend = FALSE, size = 1.2, position = position_dodge2()) +
  # putting labels
  #ggplot2::annotate('text', x = (max(data$yi) + (max(data$yi)*0.10)), y = (seq(1, group_no, 1)+0.3),
  #                 label= paste("italic(k)==", mod_table$K), parse = TRUE, hjust = "right", size = 3.5) +
  ggplot2::theme_bw() +
  ggplot2::guides(fill = "none", colour = "none") +
  ggplot2::theme(legend.position= c(0, 1), legend.justification = c(0, 1)) +
  ggplot2::theme(legend.title = element_text(size = 9)) +
  ggplot2::theme(legend.direction="horizontal") +
  ggplot2::theme(legend.background = element_blank()) +
  ggplot2::labs(y = label, x = "", size = legend) +
  ggplot2::labs(shape = condition.lab) +
  ggplot2::theme(axis.text.y = element_text(size = 10, colour ="black",
                                            hjust = 0.5,
                                            angle = 0))
# we need to do something about k[1]
plot <- plot +
  ggplot2::annotate('text', y = (max(data$yi) + (max(data$yi)*0.10)), x = (1+0.3), # here
                    label= paste("italic(k)==", mod_table$K[1]), parse = TRUE, hjust = "right", size = 3.5)


##############





