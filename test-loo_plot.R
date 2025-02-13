library(devtools)
library(metafor)


load_all()

res <- rma.mv(lnrr, lnrr_vi,
              random = ~ 1 | paper_ID,
              data = fish)


res


res$ci.lb[1]

res$ci.ub[1]


#res_loo <- leave_one_out(res, group = "paper_ID")

res_loo %>%
  ggplot2::ggplot(ggplot2::aes(x = left_out,
                               y = b,
                               ymin = ci_lb,
                               ymax = ci_ub)) +
  # Add a rectangle to highlight the overall effect and its confidence interval
  ggplot2::annotate("rect", 
                  xmin = -Inf, xmax = Inf,  # Covers all paper IDs (after coord_flip)
                  ymin = res$ci.lb[1], 
                  ymax = res$ci.ub[1], 
                  alpha = 0.5, fill = "gray")
  # Add solid line for the estimated overall effect
  ggplot2::geom_hline(yintercept = res$b[1], linetype = "solid") +
  # Add dotted lines for the estimated confidence interval
  ggplot2::geom_hline(yintercept = res$ci.lb[1], linetype = "dotted") +
  ggplot2::geom_hline(yintercept = res$ci.ub[1], linetype = "dotted") +
  # Add shaded area for the estimated confidence interval
  ggplot2::geom_pointrange(alpha = 0.7, color = "darkorange") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                 panel.grid  = ggplot2::element_blank()) +
  ggplot2::coord_flip() +


load_all()


loo_plot(res, res_loo)
