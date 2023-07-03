#' @title orchard_plot
#' @description Using a \pkg{metafor} model object of class \code{rma} or \code{rma.mv}, or a results table of class \code{orchard}, it creates an orchard plot from mean effect size estimates for all levels of a given categorical moderator, and their corresponding confidence and prediction intervals.
#' @param object model object of class \code{rma.mv}, \code{rma}, or \code{orchard} table of model results.
#' @param mod the name of a moderator. Defaults to \code{"1"} for an intercept-only model. Not needed if an \code{orchard_plot} is provided with a \code{mod_results} object of class \code{orchard}.
#' @param group The grouping variable that one wishes to plot beside total effect sizes, k. This could be study, species, or any grouping variable one wishes to present sample sizes for. Not needed if an \code{orchard_plot} is provided with a \code{mod_results} object of class \code{orchard}.
#' @param by Character vector indicating the name that predictions should be conditioned on for the levels of the moderator.
#' @param at List of levels one wishes to predict at for the corresponding varaibles in 'by'. Used when one wants marginalised means. This argument can also be used to suppress levels of the moderator when argument \code{subset = TRUE}. Provide a list as follows: \code{list(mod = c("level1", "level2"))}.
#' @param weights Used when one wants marginalised means. How to marginalize categorical variables. The default is \code{weights = "prop"}, which weights moderator level means based on their proportional representation in the data. For example, if "sex" is a moderator, and males have a larger sample size than females, then this will produce a weighted average, where males are weighted more towards the mean than females. This may not always be ideal. In the case of sex, for example, males and females are roughly equally prevalent in a population. As such, you can give the moderator levels equal weight using \code{weights = "equal"}.
#' @param xlab The effect size measure label.
#' @param N The name of the column in the data specifying the sample size, N. Defaults to \code{NULL}, so that precision is plotted instead of sample size.
#' @param alpha The level of transparency for effect sizes represented in the orchard plot.
#' @param angle The angle of y labels. The default is 90 degrees.
#' @param cb If \code{TRUE}, it uses 20 colour blind friendly colors.
#' @param k If \code{TRUE}, it displays k (number of effect sizes) on the plot.
#' @param g If \code{TRUE}, it displays g (number of grouping levels for each level of the moderator) on the plot.
#' @param transfm If set to \code{"tanh"}, a tanh transformation will be applied to effect sizes, converting Zr to a correlation or pulling in extreme values for other effect sizes (lnRR, lnCVR, SMD). Defaults to \code{"none"}.
#' @param condition.lab Label for the condition being marginalized over.
#' @param trunk.size Size of the mean, or central point.
#' @param branch.size Size of the confidence intervals.
#' @param twig.size Size of the prediction intervals.
#' @param legend.pos Where to place the legend. To remove the legend, use \code{legend.pos = "none"}.
#' @param k.pos Where to put k (number of effect sizes) on the plot. Users can specify the exact position or they can use specify \code{"right"}, \code{"left"},  or \code{"none"}. Note that numeric values (0, 0.5, 1) can also be specified and this would give greater precision.
#' @param colour Colour of effect size shapes. By default, effect sizes are colored according to the \code{mod} argument. If \code{TRUE}, they are colored according to the grouping variable
#' @param fill If \code{TRUE}, effect sizes will be filled with colours. If \code{FALSE}, they will not be filled with colours.
#' @param weights Used when one wants marginalised means. How to marginalize categorical variables. The default is \code{weights = "prop"}, which weights moderator level means based on their proportional representation in the data. For example, if "sex" is a moderator, and males have a larger sample size than females, then this will produce a weighted average, where males are weighted more towards the mean than females. This may not always be ideal. In the case of sex, for example, males and females are roughly equally prevalent in a population. As such, you can give the moderator levels equal weight using \code{weights = "equal"}.
#' @param upper Logical, defaults to \code{TRUE}, indicating that the first letter of the character string for the moderator variable should be capitalized.
#' @param flip Logical, defaults to \code{TRUE}, indicating whether the plot should be flipped.
#' @return Orchard plot
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
#' # fit a MLMR - accounting for some non-independence
#' eklof_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~ Grazer.type-1,
#' random=list(~1|ExptID, ~1|Datapoint), data=eklof)
#' results <- mod_results(eklof_MR, mod = "Grazer.type", group = "ExptID")
#' orchard_plot(results, mod = "Grazer.type",
#' group = "ExptID", xlab = "log(Response ratio) (lnRR)")
#' # or
#' orchard_plot(eklof_MR, mod = "Grazer.type", group = "ExptID",
#' xlab = "log(Response ratio) (lnRR)")
#'
#' # Example 2
#' data(lim)
#' lim$vi<- 1/(lim$N - 3)
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article,
#' ~1|Datapoint), data=lim)
#' orchard_plot(lim_MR, mod = "Phylum", group = "Article",
#' xlab = "Correlation coefficient", transfm = "tanh", N = lim$N)
#' }
#' @export

orchard_plot <- function(object, mod = "1", group, xlab, N = NULL,
                         alpha = 0.5, angle = 90, cb = TRUE, k = TRUE, g = TRUE,
                         trunk.size = 3, branch.size = 1.2, twig.size = 0.5,
                         transfm = c("none", "tanh"), condition.lab = "Condition",
                         legend.pos = c("bottom.right", "bottom.left",
                                        "top.right", "top.left",
                                        "top.out", "bottom.out",
                                        "none"), # "none" - no legends
                         k.pos = c("right", "left", "none"),
                         colour = FALSE,
                         fill = TRUE,
                         weights = "prop", by = NULL, at = NULL, upper = TRUE, flip = TRUE)
{
  ## evaluate choices, if not specified it takes the first choice
     transfm <- match.arg(NULL, choices = transfm)
  legend.pos <- match.arg(NULL, choices = legend.pos)
       k.pos <- match.arg(NULL, choices = k.pos)

	if(any(class(object) %in% c("robust.rma", "rma.mv", "rma", "rma.uni"))){

	    if(mod != "1"){
	    results <-  orchaRd::mod_results(object, mod, group,  N,
	                                        by = by, at = at, weights = weights, upper = upper)
	  } else {
	    results <-  orchaRd::mod_results(object, mod = "1", group,  N,
	                                        by = by, at = at, weights = weights, upper = upper)
	    }
	  }

	if(any(class(object) %in% c("orchard"))) {
			results <- object
	}

	mod_table <- results$mod_table

  data_trim <- results$data
  # making sure factor names match
  data_trim$moderator <- factor(data_trim$moderator, levels = mod_table$name, labels = mod_table$name)

  data_trim$scale <- (1/sqrt(data_trim[,"vi"]))
	legend <- "Precision (1/SE)"

	if(any(N != "none")){
	  data_trim$scale <- data_trim$N
		  legend <- paste0("Sample Size ($\\textit{N}$)") # we want to use italic
		  #latex2exp::TeX()
	}

	if(transfm == "tanh"){
		                   cols <- sapply(mod_table, is.numeric)
		mod_table[,cols] <- Zr_to_r(mod_table[,cols])
		data_trim$yi <- Zr_to_r(data_trim$yi)
		                  label <- xlab
	}else{
		label <- xlab
	}

	# Add in total effect sizes for each level
	 mod_table$K <- as.vector(by(data_trim, data_trim[,"moderator"], function(x) length(x[,"yi"])))

	# Add in total levels of a grouping variable (e.g., study ID) within each moderator level.
	 mod_table$g <- as.vector(num_studies(data_trim, moderator, stdy)[,2])

	 # the number of groups in a moderator & data points
	 group_no <- length(unique(mod_table[, "name"]))

	 #data_no <- nrow(data)

	# colour blind friendly colours with grey
	 cbpl <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

	 # setting fruit colour
	 if(colour == TRUE){
	   color <- as.factor(data_trim$stdy)
	   color2 <- NULL
	 }else{
	   color <- data_trim$mod
	   color2 <- mod_table$name
	 }

	 # whether we fill fruit or not
	 if(fill == TRUE){
	   fill <- color
	 }else{
	     fill <- NULL
	   }

	 # whether marginal
	 if(names(mod_table)[2] == "condition"){

	   # the number of levels in the condition
	   condition_no <- length(unique(mod_table[, "condition"]))

	   plot <- ggplot2::ggplot() +
	     # pieces of fruit (bee-swarm and bubbles)
	     ggbeeswarm::geom_quasirandom(data = data_trim, ggplot2::aes(y = yi, x = moderator, size = scale, colour = color, fill = fill), alpha=alpha, shape = 21) +

	     ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "black", alpha = alpha) +
	     # creating CI
	     ggplot2::geom_linerange(data = mod_table, ggplot2::aes(x = name, ymin = lowerCL, ymax = upperCL),
	                             size = branch.size, position = ggplot2::position_dodge2(width = 0.3)) +
	     # drowning point estimate and PI
	     ggplot2::geom_pointrange(data = mod_table, ggplot2::aes(y = estimate, x = name, ymin = lowerPR, ymax = upperPR,  shape = as.factor(condition), fill = color2), size = twig.size, position = ggplot2::position_dodge2(width = 0.3), fatten = trunk.size) +
	     # this will only work for up to 5 different conditions
	     # flipping things around (I guess we could do use the same geoms but the below is the original so we should not change)
	     ggplot2::scale_shape_manual(values =  20 + (1:condition_no))  +
	     ggplot2::theme_bw() +
	     ggplot2::guides(fill = "none", colour = "none") +
	     ggplot2::theme(legend.position= c(0, 1), legend.justification = c(0, 1)) +
	     ggplot2::theme(legend.title = ggplot2::element_text(size = 9)) +
	     ggplot2::theme(legend.direction="horizontal") +
	     ggplot2::theme(legend.background = ggplot2::element_blank()) +
	     ggplot2::labs(y = label, x = "", size = latex2exp::TeX(legend)) +
	     ggplot2::labs(shape = condition.lab) +
	     ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, colour ="black",
	                                                        hjust = 0.5,
	                                                        angle = angle))
	if(flip){
		plot <- plot + ggplot2::coord_flip()
	}

	 } else {

	  plot <- ggplot2::ggplot() +
	    # pieces of fruit (bee-swarm and bubbles)
	    ggbeeswarm::geom_quasirandom(data = data_trim, ggplot2::aes(y = yi, x = moderator, size = scale, colour = color, fill = fill), alpha=alpha, shape = 21) +

	    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "black", alpha = alpha) +
	    # creating CI
	    ggplot2::geom_linerange(data = mod_table, ggplot2::aes(x = name, ymin = lowerCL, ymax = upperCL), size = branch.size) +
	    # drowning point estimate and PI
	    ggplot2::geom_pointrange(data = mod_table, ggplot2::aes(y = estimate, x = name,  ymin = lowerPR, ymax = upperPR, fill = color2), size = twig.size, fatten = trunk.size, shape = 21) +
	    ggplot2::theme_bw() +
	    ggplot2::guides(fill = "none", colour = "none") +
	    ggplot2::theme(legend.title = ggplot2::element_text(size = 9)) +
	    ggplot2::theme(legend.direction="horizontal") +
	    ggplot2::theme(legend.background = ggplot2::element_blank()) +
	    ggplot2::labs(y = label, x = "", size = latex2exp::TeX(legend)) +
	    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, colour ="black",
	                                                       hjust = 0.5,
	                                                       angle = angle)) #+
	    #ggplot2::theme(legend.position= c(1, 0), legend.justification = c(1, 0))
	 if(flip){
		plot <- plot + ggplot2::coord_flip()
	}
	 }

	   # adding legend
	 if(legend.pos == "bottom.right"){
	   plot <- plot + ggplot2::theme(legend.position= c(1, 0), legend.justification = c(1, 0))
	 } else if ( legend.pos == "bottom.left") {
	   plot <- plot + ggplot2::theme(legend.position= c(0, 0), legend.justification = c(0, 0))
	 } else if ( legend.pos == "top.right") {
	   plot <- plot + ggplot2::theme(legend.position= c(1, 1), legend.justification = c(1, 1))
	 } else if (legend.pos == "top.left") {
	   plot <- plot + ggplot2::theme(legend.position= c(0, 1), legend.justification = c(0, 1))
	 } else if (legend.pos == "top.out") {
	   plot <- plot + ggplot2::theme(legend.position="top")
	 } else if (legend.pos == "bottom.out") {
	   plot <- plot + ggplot2::theme(legend.position="bottom")
	 } else if (legend.pos == "none") {
	   plot <- plot + ggplot2::theme(legend.position="none")
	 }

	  # putting colors in
	  if(cb == TRUE){
	    plot <- plot +
	      ggplot2::scale_fill_manual(values = cbpl) +
	      ggplot2::scale_colour_manual(values = cbpl)
	  }

	  # putting k and g in
	  if(k == TRUE && g == FALSE && k.pos == "right"){
	    plot <- plot +
	      ggplot2::annotate('text', y = (max(data_trim$yi) + (max(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
	                        label= paste("italic(k)==", mod_table$K[1:group_no]), parse = TRUE, hjust = "right", size = 3.5)
	  } else if(k == TRUE && g == FALSE && k.pos == "left") {
	    plot <- plot +  ggplot2::annotate('text', y = (min(data_trim$yi) + (min(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
	                                      label= paste("italic(k)==", mod_table$K[1:group_no]), parse = TRUE, hjust = "left", size = 3.5)
	  } else if (k == TRUE && g == TRUE && k.pos == "right"){
	    # get group numbers for moderator
	    plot <- plot + ggplot2::annotate('text', y = (max(data_trim$yi) + (max(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
	                        label= paste("italic(k)==", mod_table$K[1:group_no], "~","(", mod_table$g[1:group_no], ")"),
	                        parse = TRUE, hjust = "right", size = 3.5)
	  } else if (k == TRUE && g == TRUE && k.pos == "left"){
	    # get group numbers for moderator
	    plot <- plot + ggplot2::annotate('text',  y = (min(data_trim$yi) + (min(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
	                        label= paste("italic(k)==", mod_table$K[1:group_no], "~","(", mod_table$g[1:group_no], ")"),
	                        parse = TRUE, hjust = "left", size = 3.5)
	  } else if (k == TRUE && g == FALSE && k.pos%in%c('right','left','none')==FALSE) {
	    # get group numbers for moderator
	    plot <- plot + ggplot2::annotate("text", y = k.pos, x = (seq(1, group_no,
	                                                                1) + 0.3), label = paste("italic(k)==", mod_table$K[1:group_no]),
	                                    parse = TRUE, size = 3.5)
	  } else if (k == TRUE && g == TRUE && k.pos%in%c('right','left','none')==FALSE) {
	    # get group numbers for moderator
	    plot <- plot + ggplot2::annotate("text", y = k.pos, x = (seq(1, group_no,
	                                                                1) + 0.3), label = paste("italic(k)==", mod_table$K[1:group_no],
	                                                                                         "~", "(", mod_table$g[1:group_no], ")"),
	                                    parse = TRUE, size = 3.5)
    }
	  return(plot)
}


#' @title Zr_to_r
#' @description Converts Zr back to r (Pearson's correlation coefficient)
#' @param df data frame of results of class 'orchard'
#' @return A data frame containing all the model results including mean effect size estimate, confidence and prediction intervals with estimates converted back to r
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @export

Zr_to_r <- function(df){
	return(sapply(df, tanh))
}

