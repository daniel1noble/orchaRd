#' @title caterpillars
#' @description Using a \pkg{metafor} model object of class \code{rma} or \code{rma.mv} or a results table of class \code{orchard}, this function produces a caterpillar plot from mean effect size estimates for all levels of a given categorical moderator, their corresponding confidence intervals, and prediction intervals.
#' @param object Model object of class \code{rma.mv}, \code{rma} or \code{orchard} table of model results
#' @param mod The name of a moderator variable. Otherwise, "1" for an intercept-only model.
#' @param group The grouping variable that one wishes to plot beside total effect sizes, k. This could be study, species, or whatever other grouping variable one wishes to present sample sizes for.
#' @param xlab The effect size measure label.
#' @param overall Logical, indicating whether to re-label "Intrcpt" (the default label from \code{rma} or \code{rma.mv} intercept only models or meta-analyses) to "Overall". Defaults to \code{TRUE}.
#' @param transfm If set to \code{"tanh"}, a tanh transformation will be applied to effect sizes, converting Zr to a correlation or pulling in extreme values for other effect sizes (lnRR, lnCVR, SMD).  \code{"invlogit"} can be used to convert lnRR to the inverse logit scale. \code{"percentr"} can convert to the percentage change scale when using response ratios and \code{"percent"} can convert to the percentage change scale of an log transformed effect size. Defaults to \code{"none"}.
#' @param n_transfm The vector of sample sizes for each effect size estimate. This is used when \code{transfm = "inv_ft"}. Defaults to NULL.
#' @param colerrorbar Colour of the error bar in the caterpillars plot. Defaults to hex code - "#00CD00". 
#' @param colpoint Point estimate colour in the caterpillars plot. Defaults to hex code - "#FFD700".
#' @param colpoly Polygon colour in the caterpillars plot. Defaults to "red".
#' @param k If \code{TRUE}, it displays k (number of effect sizes) on the plot.
#' @param g If \code{TRUE}, it displays g (number of grouping levels for each level of the moderator) on the plot.
#' @param by Used when one wants marginalised means. Character vector indicating the name that predictions should be conditioned on for the levels of the moderator.
#' @param at Used when one wants marginalised means. List of levels one wishes to predict at for the corresponding variables in \code{by}, but is not presented in output. Provide a list as follows: \code{list(mod = c("level1", "level2"))}.
#' @param weights Used when one wants marginalised means. How to marginalize categorical variables. The default is \code{weights = "prop"}, which weights means for moderator levels based on their proportional representation in the data. For example, if "sex" is a moderator, and males have a larger sample size than females, then this will produce a weighted average, where males are weighted more towards the mean than females. This may not always be ideal. In the case of sex, for example, males and females are roughly equally prevalent in a population. As such, you can give the moderator levels equal weight using \code{weights = "equal"}.
#' @return Caterpillars plot
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' data(eklof)
#' eklof<-metafor::escalc(measure="ROM", n1i=N_control, sd1i=SD_control,
#' m1i=mean_control, n2i=N_treatment, sd2i=SD_treatment, m2i=mean_treatment,
#' data=eklof)
#' # Add the unit level predictor
#' eklof$Datapoint<-as.factor(seq(1, dim(eklof)[1], 1))
#' # fit a MLMR - accouting for some non-independence
#' eklof_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~ Grazer.type-1, random=list(~1|ExptID,
#' ~1|Datapoint), data=eklof)
#' results <- mod_results(eklof_MR, mod = "Grazer.type",  group = "First.author")
#' caterpillars(results, mod = "Grazer.type",
#' group = "First.author", xlab = "log(Response ratio) (lnRR)", g = FALSE)
#'
#' # Example 2
#' data(lim)
#' lim$vi<- 1/(lim$N - 3)
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article,
#' ~1|Datapoint), data=lim)
#' results_lim <- mod_results(lim_MR, mod = "Phylum", group = "Article")
#' caterpillars(results_lim, mod = "Phylum",
#' group = "Article", xlab = "Correlation coefficient", transfm = "tanh")
#' }
#' @export

caterpillars <- function(object, mod = "1",  group, xlab, overall = TRUE, transfm = c("none", "tanh", "invlogit", "percent", "percentr", "inv_ft"), n_transfm = NULL, colerrorbar = "#00CD00", colpoint = "#FFD700", colpoly = "red",  k = TRUE, g = TRUE, at = NULL, by = NULL, weights = "prop") {

  ## evaluate choices
  transfm <- match.arg(transfm) # if not specified it takes the first choice

  results <- .get_results(object, mod, group, N = NULL, by, at, weights)

  if (transfm != "none") {
    results <- transform_mod_results(results, transfm, n_transfm)
  }

  # Meta-analytic results
  mod_table <- results$mod_table
  # Data set
  data <- results$data
  data$lower <- data$yi - stats::qnorm(0.975)*sqrt(data$vi)
  data$upper <- data$yi + stats::qnorm(0.975)*sqrt(data$vi)

  label <- xlab


  if("Intrcpt" %in% mod_table$name){
    mod_table$name <- replace(as.vector(mod_table$name), which(mod_table$name == "Intrcpt"), "Overall")
  }

  # adding moderator names
  data$moderator <- factor(data$moderator, labels = mod_table$name)

  # data frame for the meta-analytic results
  mod_table$K <- as.vector(by(data, data[,"moderator"], function(x) length(x[,"yi"])))

  # Add in total levels of a grouping variable (e.g., study ID) within each moderator level.
  mod_table$g <- as.vector(num_studies(data, moderator, stdy)[,2])

  # the number of groups in a moderator & data points
  group_no <- nrow(mod_table)
  data_no <- nrow(data)

  # use dplyr here - need to change....
  # Dan can you make this basic R code - maybe I got it
  # data <- data[order(data$moderator, -data$yi),]
  data <- data %>% dplyr::group_by(moderator) %>% dplyr::arrange(moderator, dplyr::desc(yi)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Y = 1:data_no +
             unlist(mapply(function(x, y) rep(x*6 , y) , x = 1:group_no, y = mod_table$K))
    ) %>%
    data.frame()

  # mod ID
  mod_table$Y <- data %>% dplyr::group_by(moderator) %>%
    dplyr::summarise(Y = dplyr::first(Y)) %>%
    dplyr::select(Y) %>% t() %>% as.vector() -2

  # preparing for diamons for summary
  # modified from internal_viz_classicforest() from the R package, metaviz
  sum_data <- data.frame("x.diamond" = c(mod_table$lowerCL,
                                         mod_table$estimate ,
                                         mod_table$upperCL,
                                         mod_table$estimate ),
                         "y.diamond" = c(mod_table$Y,
                                         mod_table$Y + 1.2,
                                         mod_table$Y,
                                         mod_table$Y - 1.2),
                         "moderator" = rep(mod_table$name, times = 4)
  )

  # make caterpillars plot
  plot <- ggplot2::ggplot(data = data, ggplot2::aes(x = yi, y = Y)) +
    # 95 % CI
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lower, xmax = upper),
                           orientation = "y", colour = colerrorbar, height = 0, show.legend = FALSE, linewidth = 0.5, alpha = 0.6) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.5) +
    # creating dots for point estimates
    ggplot2::geom_point(colour = colpoint, size = 1) +
    # creating 95% prediction intervals
    ggplot2::geom_segment(data = mod_table, ggplot2::aes(x = lowerPR, y = Y, xend = upperPR, yend = Y, group = name)) +
    # creating diamonsts (95% CI)
    ggplot2::geom_polygon(data = sum_data, ggplot2::aes(x = x.diamond, y = y.diamond, group = moderator), fill = colpoly) +
    #ggplot2::facet_wrap(~moderator, scales = "free_y", nrow = GN,  strip.position = "left") + # using facet_wrap - does not really work well
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0, size = 8),# margin = margin(t=15, r=15, b=15, l=15)),
                   strip.background = ggplot2::element_rect(colour = NULL,
                                                   linetype = "blank",
                                                   fill = "gray90"),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::labs(x = label, y = "", parse = TRUE) +

    ggplot2::annotate('text', x = min(data$lower)*0.975, y = mod_table$Y,
                      label= mod_table$name, hjust = "left", size = 3.5) +
    ggplot2::coord_cartesian(xlim = c(min(data$lower)*1.05, max(data$upper)*1.05),
                    ylim = c((min(data$Y)-10), (max(data$Y)+4))
                    , expand = F)

  # putting k in
  if(k == TRUE && g == FALSE){
    plot <- plot +
      ggplot2::annotate('text', x = max(data$upper)*0.975, y = mod_table$Y-1.7,
                        label= paste("italic(k)==", mod_table$K), parse = TRUE, hjust = "right", size = 3.5)
  }

  # putting groups
  if(k == TRUE && g == TRUE){
    # get group numbers for moderator
    plot <- plot +
      ggplot2::annotate('text', x = max(data$upper)*0.975, y = mod_table$Y-1.7,
                        label= paste("italic(k)==", mod_table$K[1:group_no], "~~","(", mod_table$g[1:group_no], ")"), parse = TRUE, hjust = "right", size = 3.5)
  }


  return(plot)
}
