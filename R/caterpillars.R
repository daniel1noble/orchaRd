#' @title caterpillars
#' @description Using a metafor model object of class rma or rma.mv or a results table of class orchard, it creates a an caterpillars plot from mean effect size estimates for all levels of a given categorical moderator, their corresponding confidence intervals and prediction intervals
#' @param object Model object of class 'rma.mv', 'rma' or 'orchard' table of model results
#' @param mod The name of a moderator. Otherwise, "1" for intercept only model.
#' @param group The grouping variable that one wishes to plot beside total effect sizes, k. This could be study, species or whatever other grouping variable one wishes to present sample sizes.
#' @param data The data frame used to fit the rma.mv model object
#' @param xlab The effect size measure label.
#' @param overall Logical indicating whether to relabel "Intrcpt" (the default label from rma or rma.mv intercept only models or meta-analyses) to "Overall".
#' @param transfm If set to "tanh", a tanh transformation will be applied to effect sizes, converting Zr will to a correlation or pulling in extreme, values for other effect sizes (lnRR, lnCVR, SMD). If "none" is chosen then it will default.
#' @param k If TRUE, it displays k (number of effect sizes) on the plot
#' @param g If TRUE, it displays g (number of grouping levels for each level of the moderator) on the plot
#' @param by Used when one wants marginalised means. The 'condition' variable that one wishes to have the mean for the moderator vary.
#' @param at Used when one wants marginalised means. The 'condition' that one wishes to calculate the means at, but is not presented in output
#' @param weights Used when one wants marginalised means. How to marginalize categorical variables. The default is weights = "prop", which wights means for moderator levels based on their proportional representation in the data. For example, if "sex" is a moderator, and males have a larger sample size than females, then this will produce a weighted average, where males are weighted more towards the mean than females. This may not always be ideal. IN the case if sex, for example, males and females are roughly equally prevalent in a population. As such, you can give the moderator levels equal weight using weights = "equal".
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
#' results <- mod_results(eklof_MR, mod = "Grazer.type", data = eklof, group = "First.author")
#' caterpillars(results, mod = "Grazer.type", data = eklof,
#' group = "First.author", xlab = "log(Response ratio) (lnRR)", g = FALSE)
#'
#' # Example 2
#' data(lim)
#' lim$vi<- 1/(lim$N - 3)
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article,
#' ~1|Datapoint), data=lim)
#' results_lim <- mod_results(lim_MR, mod = "Phylum", data = lim, group = "Article")
#' caterpillars(results_lim, mod = "Phylum", data = lim,
#' group = "Article", xlab = "Correlaiton coefficent", transfm = "tanh")
#' }
#' @export

caterpillars <- function(object, mod = "1", data, group, xlab, overall = TRUE, transfm = c("none", "tanh"), k = TRUE, g = TRUE, at = NULL, by = NULL, weights = "prop") {

  if(any(class(object) %in% c("rma.mv", "rma"))){

    if(mod != "1"){
      results <-  orchaRd::mod_results(object, mod, group, data,
                                       by = by, at = at, weights = weights)
    } else {
      results <-  orchaRd::mod_results(object, mod = "1", group, data,
                                       by = by, at = at, weights = weights)
    }
  }

  if (any(class(object) %in% c("orchard"))) {results <- object}

       ## evaluate choices
  transfm <- match.arg(transfm) # if not sepcificed it takes the first choice

  # meta-analytic results
  mod_table <- results$mod_table

  # data set
  data <- results$data
  data$lower <- data$yi - stats::qnorm(0.975)*sqrt(data$vi)
  data$upper <- data$yi + stats::qnorm(0.975)*sqrt(data$vi)

  if(transfm == "tanh"){
    cols <- sapply(mod_table, is.numeric)
    mod_table[,cols] <- Zr_to_r(mod_table[,cols])
    data$yi <- Zr_to_r(data$yi)
    data$lower <- Zr_to_r(data$lower)
    data$upper <- Zr_to_r(data$upper)
    label <- xlab
  }else{
    label <- xlab
  }

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
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper),
                            colour = "#00CD00", height = 0, show.legend = FALSE, size = 0.5, alpha = 0.6) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.5) +
    # creating dots for point estimates
    ggplot2::geom_point(colour = "#FFD700", size = 1) +
    # creating 95% prediction intervals
    ggplot2::geom_segment(data = mod_table, ggplot2::aes(x = lowerPR, y = Y, xend = upperPR, yend = Y, group = name)) +
    # creating diamonsts (95% CI)
    ggplot2::geom_polygon(data = sum_data, ggplot2::aes(x = x.diamond, y = y.diamond, group = moderator), fill = "red") +
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
