# test

remotes::install_github("rvlenth/emmeans")
library(orchaRd)
library(metafor)
library(emmeans)

data(fish)
warm_dat <- fish

# model
model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list(~1 | group_ID, ~1 | es_ID), mods = ~ experimental_design + trait.type + deg_dif + treat_end_days, method = "REML", test = "t", data = warm_dat,                               control=list(optimizer="optim", optmethod="Nelder-Mead"))

# margainl ovarll
overall <- marginalised_means(model, data = warm_dat)
overall2 <- marginalised_means(model, data = warm_dat, mod = "1", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif",  weights = "prop")


# marginalsed stuff
across_trait <- marginalised_means(model, data = warm_dat, mod = "trait.type")
across_trait_by_degree_diff <- marginalised_means(model, data = warm_dat, mod = "trait.type", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif")
across_trait_by_degree_diff_at_treat_end_days10 <- marginalised_means(model, data = warm_dat, mod = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = 10), by = "deg_dif")
across_trait_by_degree_diff_at_treat_end_days10And50 <- marginalised_means(model, data = warm_dat, mod = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), by = "deg_dif")


across_trait_by_treat_end_days10And50 <- marginalised_means(model, data = warm_dat, mod = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), by = "treat_end_days")

#
across_trait_by_treat_end_days10And50_ordinaryMM <- marginalised_means(model, data = warm_dat, mod = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), by = "treat_end_days", weights = "prop")

# current get_data

data <-overall2$data
mod_table <- overall2$mod_table
label <- "Test"
angle <- 90
alpha <- 0.5

plot <- ggplot2::ggplot() +
  # pieces of fruit (bee-swarm and bubbles)
  ggbeeswarm::geom_quasirandom(data = data, aes(y = yi, x = moderator, size = (1/sqrt(data[,"vi"])), colour = moderator), alpha=alpha) +

  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "black", alpha = alpha) +
  # creating dots + CI and PI
  ggplot2::geom_linerange(data = mod_table, aes(x = name, ymin = lowerCL, ymax = upperCL), size = 1.2, position = position_dodge2(width = 0.3)) +
  ggplot2::geom_pointrange(data = mod_table, aes(y = estimate, x = name, fill = name, ymin = lowerPR, ymax = upperPR), size = 0.6, shape = 21, position = position_dodge2(width = 0.3)) +
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
  ggplot2::labs(y = "lnRR", x = "", size = "N") +
  ggplot2::theme(axis.text.y = element_text(size = 10, colour ="black",
                                            hjust = 0.5,
                                            angle = 0))


