

#' @title pub_bias_plot
#' @description This function adds the bias corrected mean and confidence intervals to an existing orchard plot that displays the overall meta-analytic mean effect size.
#' @param plot The orchard plot object to add the bias corrected mean and confidence intervals to. This plot needs to be a plot that displays the raw and a single meta-analytic mean (overall mean), confidence interval and predcition interval. 
#' @param fe_model The meta-analytic model (rma object) produced from a two-step correction (sensu Yang et al. 2023) (Step 1: Fixed effect model) with a robust corretion to correct the meta-analytic mean for publication bias (selection bias) and the dependency in the data. 
#' @param v_model An optional argument. The meta-analytic model (rma object) to deal with publication bias using the meta-regression approach proposed by Nakagawa et al. 2023. This model can have fixed and/or random effects, with one random effect being sampling error (se or v). The intercept from this model can be considred the corrected meta-analytic mean when sample size is infinite or sampling error is zero.
#' @param col The colour of the mean and confidence intervals.
#' @param plotadj The adjustment to the x-axis position of the mean and confidence intervals.
#' @param textadj The adjustment to the y-axis position of the mean and confidence intervals for the text displaying the type of correction. 
#' @param trunk.size Size of the mean, or central point.
#' @param branch.size Size of the confidence intervals.
#' @return An orchard plot with the corrected meta-analytic mean and confidence intervals added.
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' 	# Data
#' 		data(english)
#' 		# We need to calculate the effect sizes, in this case d
#' 		english <- escalc(measure = "SMD", n1i = NStartControl, sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E, m2i = MeanE, 
#' 						var.names=c("SMD","vSMD"),
#' 						data = english)
#' 
#' 	# Our MLMA model
#' 		english_MA1 <- rma.mv(yi = SMD, V = vSMD, random = list( ~ 1 | StudyNo, ~ 1 | EffectID),test = "t", data = english)
#' 
#' 	# Step 1: Fit the fixed effect model
#' 		english_MA2 <- rma.mv(yi = SMD, V = vSMD, data = english, test = "t")
#' 
#' 		english_MA3 <- rma(yi = SMD, vi = vSMD, data = english, test = "t", method = "FE")
#' 
#' 	# Step 2: Correct for dependency 
#' 		english_MA2_1 <- robust(english_MA2, cluster = english$StudyNo, clubSandwich=TRUE)
#' 
#' 	# Step 3: Testing modified eggers. Need intercept
#' 		english_MA4 <- rma.mv(yi = SMD, V = vSMD, mod = ~vSMD, random = list( ~ 1 | StudyNo, ~ 1 | EffectID),test = "t", data = english)
#' 
#' 	# Now plot the results
#' 
#' 		plot <- orchard_plot(english_MA1, group = "StudyNo",  xlab = "Standardized Mean Difference")
#' 		plot2 <- pub_bias_plot(plot, english_MA2_1)
#' 		plot3 <- pub_bias_plot(plot, english_MA2_1, english_MA4)
#' 		}
#' @export
  	pub_bias_plot <- function(plot, fe_model, v_model = NULL, col = c("red", "blue"), plotadj = -0.05, textadj = 0.05, branch.size = 1.2, trunk.size = 3){
		      
		# Add check to make sure it's an intercept ONLY model being added. Message to user if not.
		if(length(fe_model$b) > 1){
			stop("The model you are trying to add to the plot is not an intercept only model. Please ensure you have fit an intercept only meta-analysis. See vignette for details: https://daniel1noble.github.io/orchaRd/")
		}
		
		# Get the predictions from the final model and create a label for the plot
			pub_bias_data <- get_ints_dat(fe_model, type = "br")

		if(is.null(v_model)){
			
			# Add to Existing Orchard Plot
			plot + geom_pub_stats_yang(pub_bias_data, plotadj = plotadj, textadj = textadj, branch.size = branch.size, trunk.size = trunk.size) 
		
		} else{
			# Extract the corrected meta-analytic mean and CI
			pub_bias_data2 <- get_ints_dat(v_model, type = "bc")

			plot + geom_pub_stats_yang(pub_bias_data, plotadj = plotadj, textadj = textadj, branch.size = branch.size, trunk.size = trunk.size) + geom_pub_stats_naka(pub_bias_data2, plotadj = plotadj, textadj = textadj, branch.size = branch.size, trunk.size = trunk.size) 
		}		
	}


#####################
## Helper functioons
#####################


#' @title geom_pub_stats_yang
#' @description This function adds a corrected meta-analytic mean, sensu Yang et al. 2023, confidence interval and text annotation to an intercept only orchard plot.
#' @param data The data frame containing the corrected meta-analytic mean and confidence intervals.
#' @param col The colour of the mean and confidence intervals.
#' @param plotadj The adjustment to the x-axis position of the mean and confidence intervals.
#' @param textadj The adjustment to the y-axis position of the mean and confidence intervals for the text displaying the type of correction.
#' @param branch.size Size of the confidence intervals.
#' @param trunk.size Size of the mean, or central point.
#' @return A list of ggplot2 objects to be added to the orchard plot.
#' @author Daniel Noble - daniel.noble@anu.edu.au

	geom_pub_stats_yang <-  function(data, col = "red", plotadj = -0.05, textadj = 0.05, branch.size = 1.2, trunk.size = 3){
		list(ggplot2::geom_point(data = data[[1]], ggplot2::aes(x = name, y = pred), color = col, shape = "diamond", position = position_nudge(plotadj), size = trunk.size), 
				ggplot2::geom_linerange(data = data[[1]], ggplot2::aes(x = name, ymin = ci.lb, ymax = ci.ub), color = col, position = position_nudge(plotadj), size = branch.size),
					ggplot2::annotate("text", x = 1+plotadj-textadj, y = data[[1]]$pred+textadj, label = data[[2]], color = col, size = 4, hjust = data[[1]]$ci.ub -0.2)	
		)
	}

#' @title geom_pub_stats_naka
#' @description This function adds a corrected meta-analytic mean, sensu Nakagawa et al. 2022, confidence interval and text annotation to an intercept only orchard plot.
#' @param data The data frame containing the corrected meta-analytic mean and confidence intervals.
#' @param col The colour of the mean and confidence intervals.
#' @param plotadj The adjustment to the x-axis position of the mean and confidence intervals.
#' @param textadj The adjustment to the y-axis position of the mean and confidence intervals for the text displaying the type of correction.
#' @param branch.size Size of the confidence intervals.
#' @param trunk.size Size of the mean, or central point.
#' @return A list of ggplot2 objects to be added to the orchard plot.
#' @author Daniel Noble - daniel.noble@anu.edu.au

	geom_pub_stats_naka <- function(data, col = "blue", plotadj = -0.05, textadj = 0.05, branch.size = 1.2, trunk.size = 3) {
					list(ggplot2::geom_point(data = data[[1]], ggplot2::aes(x = name, y = pred), color = col, shape = "diamond", position = position_nudge(abs(plotadj)), size = trunk.size), 
					ggplot2::geom_linerange(data = data[[1]], ggplot2::aes(x = name, ymin = ci.lb, ymax = ci.ub), color = col, position = position_nudge(abs(plotadj)), size = branch.size), 
					ggplot2::annotate("text", x = 1+abs(plotadj)+textadj, y = data[[1]]$pred-textadj, label = data[[2]], color = col, size = 4, hjust = data[[1]]$ci.ub +0.2))
	}

#' @title get_ints_dat
#' @description This function extracts the corrected meta-analytic mean and confidence intervals from a model object.
#' @param model The rma model object containing the corrected meta-analytic mean and confidence intervals.
#' @param type The type of correction to extract the corrected meta-analytic mean and confidence intervals from. "br" (i.e., Bias Robust) for Yang et al. 2023, "bc" (i.e., Bias-Corrected) for Nakagawa et al. 2023.
#' @return A list containing the corrected meta-analytic mean and confidence intervals, and a label for the plot.
#' @author Daniel Noble - daniel.noble@anu.edu.au

	get_ints_dat <- function(model, type = c("bc", "br")){
			# Extract the corrected meta-analytic mean and CI
				type = match.arg(type)
				
				dat <- data.frame(name  =  "Intrcpt", 
								pred = model$b["intrcpt",], 
								ci.lb = model$ci.lb[1], 
								ci.ub = model$ci.ub[1])
				if(type == "bc"){
				lab <- paste0("Bias Corrected Estimate: ", round(dat$pred, 2), 
									", 95% CI [", round(dat$ci.lb, 2), ",", round(dat$ci.ub, 2), "]")
				}

				if(type == "br"){
				lab <- paste0("Bias Robust Estimate: ", round(dat$pred, 2), 
									", 95% CI [", round(dat$ci.lb, 2), ",", round(dat$ci.ub, 2), "]")
	}

	return(list(dat, lab))
	}