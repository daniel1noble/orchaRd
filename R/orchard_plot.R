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

#' @title orchard_plot
#' @description Using a metafor model object of class rma or rma.mv or a results table of class orchard, it creates a an orchard plot from mean effect size estimates for all levels of a given categorical moderator, their corresponding confidence intervals and prediction intervals
#' @param object model object of class 'rma.mv', 'rma' or 'orchard' table of model results
#' @param mod the name of a moderator. Otherwise, "Int" for intercept only model.
#' @param xlab The effect size measure label.
#' @param N  The vector of sample size which an effect size is based on. If default, we use precision (the inverse of sampling standard error)
#' @param alpha The level of transparency for pieces of fruit (effect size)
#' @param angle The angle of y labels. The default is 90 degrees
#' @param cb If TRUE, it uses 8 colour blind friendly colors (7 colours plus grey)
#' @param k If TRUE, it displays k (number of effect sizes) on the plot
#' @param transfm If set to "tanh", a tanh transformation will be applied to effect sizes, converting Zr will to a correlation or pulling in extreme values for other effect sizes (lnRR, lnCVR, SMD). If "none" is chosen then it will default to
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
#' # fit a MLMR - accouting for some non-independence
#' eklof_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~ Grazer.type-1, random=list(~1|ExptID,
#' ~1|Datapoint), data=eklof)
#' results <- mod_results(eklof_MR, mod = "Grazer.type")
#' orchard_plot(results, mod = Grazer.type, xlab = "log(Response ratio) (lnRR)")
#' # or
#' orchard_plot(eklof_MR, mod = Grazer.type, xlab = "log(Response ratio) (lnRR)")
#'
#' # Example 2
#' data(lim)
#' lim$vi<- 1/(lim$N - 3)
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article,
#' ~1|Datapoint), data=lim)
#' orchard_plot(lim_MR, mod = "Phylum", xlab = "Correlaiton coefficent", transfm = "tanh", N = lim$N)
#' }
#' @export

# Note that cb = TRUE - we run out colour quite quickly
# TODO - we should say more than 8 groups (cb) we use cb = FALSE
# TODO- we can only do up to 5 conditions for drowing pruposes - we
# TODO - make it possible to where we want ot put legend - 1, 2, 3, 4 (top.right, top.left, bottom.right, bottom.left)
# TODO - making we can turn on and off legends too?? - I think we should
# TODO - we do not really need "Int" for marginal_means
# TODO - we can get ride of "mod" - I think discuss with Dan
orchard_plot <- function(object, mod = "Int", xlab, N = "none",
                         alpha = 0.5, angle = 90, cb = TRUE, k = TRUE,
                         transfm = c("none", "tanh"), condition.lab = "Condition")
                         #legend.pos = c("top.left", "", "", "", "top.out", "bottom.out"))
  {

  ## evaluate choices
  transfm <- match.arg(transfm) # if not sepcificed it takes the first choice

	if(any(class(object) %in% c("rma.mv", "rma"))){
		if(mod != "Int"){
			object <- mod_results(object, mod)
		} else{
			object <- mod_results(object, mod = "Int")
			}
	}

	mod_table <- object$mod_table

  data <- object$data
  data$moderator <- factor(data$moderator, levels = mod_table$name, labels = mod_table$name)

	data$scale <- (1/sqrt(data[,"vi"]))
	legend <- "Precision (1/SE)"

	if(any(N != "none")){
		  data$scale <- N
		  legend <- "Sample Size (N)"
	}

	if(transfm == "tanh"){
		                   cols <- sapply(mod_table, is.numeric)
		mod_table[,cols] <- Zr_to_r(mod_table[,cols])
		                data$yi <- Zr_to_r(data$yi)
		                  label <- xlab
	}else{
		label <- xlab
	}

	 mod_table$K <- as.vector(by(data, data[,"moderator"], function(x) length(x[,"yi"])))

	 # the number of groups in a moderator & data points
	 group_no <- length(unique(mod_table[, "name"]))

	 #data_no <- nrow(data)

	# colour blind friendly colours with grey
	 cbpl <- c("#E69F00","#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#56B4E9", "#999999")

	 # whether marginal
	 if(names(mod_table)[2] == "condition"){

	   # the number of levels in the condition
	   condition_no <- length(unique(mod_table[, "condition"]))

	   plot <- ggplot2::ggplot() +
	     # pieces of fruit (bee-swarm and bubbles)
	     ggbeeswarm::geom_quasirandom(data = data, ggplot2::aes(y = yi, x = moderator, size = scale, colour = moderator), alpha=alpha) +

	     ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "black", alpha = alpha) +
	     # creating CI
	     ggplot2::geom_linerange(data = mod_table, ggplot2::aes(x = name, ymin = lowerCL, ymax = upperCL),
	                             size = 1.2, position = ggplot2::position_dodge2(width = 0.3)) +
	     # drowning point estimate and PI
	     ggplot2::geom_pointrange(data = mod_table, ggplot2::aes(y = estimate, x = name, ymin = lowerPR, ymax = upperPR,  shape = as.factor(condition), fill = name),
	                              size = 0.5, position = ggplot2::position_dodge2(width = 0.3)) +
	     # this will only work for up to 5 different conditions
	     # flipping things around (I guess we could do use the same geoms but the below is the oiginal so we should not change)
	     ggplot2::scale_shape_manual(values =  20 + (1:condition_no)) + coord_flip() +
	     ggplot2::theme_bw() +
	     ggplot2::guides(fill = "none", colour = "none") +
	     ggplot2::theme(legend.position= c(0, 1), legend.justification = c(0, 1)) +
	     ggplot2::theme(legend.title = ggplot2::element_text(size = 9)) +
	     ggplot2::theme(legend.direction="horizontal") +
	     ggplot2::theme(legend.background = ggplot2::element_blank()) +
	     ggplot2::labs(y = label, x = "", size = legend) +
	     ggplot2::labs(shape = condition.lab) +
	     ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, colour ="black",
	                                               hjust = 0.5,
	                                               angle = angle))

	   plot <- plot +
	     ggplot2::annotate('text', y = (max(data$yi) + (max(data$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
	                       label= paste("italic(k)==", mod_table$K[1:group_no]), parse = TRUE, hjust = "right", size = 3.5)

	 }else{

	# Make the orchard plot
	  plot <- ggplot2::ggplot(data = mod_table, ggplot2::aes(x = estimate, y = name)) +
	    # pieces of fruit (bee-swarm and bubbles)
	  	ggbeeswarm::geom_quasirandom(data = data, ggplot2::aes(x = yi, y = moderator, size = scale, colour = moderator), groupOnX = FALSE, alpha=alpha) +
	  	# 95 %prediction interval (PI): twigs
	  	ggplot2::geom_errorbarh(ggplot2::aes(xmin = lowerPR, xmax = upperPR),  height = 0, show.legend = FALSE, size = 0.5, alpha = 0.6) +
	  	# 95 %CI: branches
	  	ggplot2::geom_errorbarh(ggplot2::aes(xmin = lowerCL, xmax = upperCL),  height = 0, show.legend = FALSE, size = 1.2) +
	  	ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = alpha) +
	  	# creating dots for truncks
	  	ggplot2::geom_point(ggplot2::aes(fill = name), size = 3, shape = 21) +
	  	# putting labels
	  	#ggplot2::annotate('text', x = (max(data$yi) + (max(data$yi)*0.10)), y = (seq(1, group_no, 1)+0.3),
	  	#                 label= paste("italic(k)==", mod_table$K), parse = TRUE, hjust = "right", size = 3.5) +
	  	ggplot2::theme_bw() +
	    ggplot2::guides(fill = "none", colour = "none") +
	    ggplot2::theme(legend.position= c(1, 0), legend.justification = c(1, 0)) +
	    ggplot2::theme(legend.title = ggplot2::element_text(size = 9)) +
	    ggplot2::theme(legend.direction="horizontal") +
	    ggplot2::theme(legend.background = ggplot2::element_blank()) +
	  	ggplot2::labs(x = label, y = "", size = legend) +
	    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, colour ="black",
	                                              hjust = 0.5,
	                                              angle = angle))

	  # putting k in
	  if(k == TRUE){
	    plot <- plot +
	      ggplot2::annotate('text', x = (max(data$yi) + (max(data$yi)*0.10)), y = (seq(1, group_no, 1)+0.3),
	                        label= paste("italic(k)==", mod_table$K), parse = TRUE, hjust = "right", size = 3.5)
	  }

	 }
	  # putting colors in
	  if(cb == TRUE){
	    plot <- plot +
	      ggplot2::scale_fill_manual(values=cbpl) +
	      ggplot2::scale_colour_manual(values=cbpl)
	  }

	  return(plot)
}



