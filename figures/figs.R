

###################
## Setup
###################
	# Load packages
		devtools::install_github("daniel1noble/orchaRd", force = TRUE)
		library(orchaRd)
		library(metafor)
		library(tidyverse)
		library(latex2exp)
		library(patchwork)
		library(flextable)

	# Function for rounding  tables
		round_df <- function(df, digits) {
		  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))

		  df[,nums] <- round(df[,nums], digits = digits)

		  (df)
		}

	# Load the data
		data("pottier")

	# Load the phylogenetic correlation matrix
		data("phylo_matrix")

	# Load the sampling variance matrix
		data("VCV_dARR")

	# Heteroscedasticity modeled at the effect size level
		mod.habitat_het <-  rma.mv(yi=dARR, V=VCV_dARR,mods= ~habitat,method="REML",test="t",dfs="contain",
		                           random=list(~1|species_ID,~1|phylogeny,~habitat|es_ID),struct="HCS", rho=0, R =list(phylogeny = phylo_matrix),data=pottier)
		mod.habitat_het

		mod.habitat_hom <-  rma.mv(yi=dARR, V=VCV_dARR,mods= ~habitat,method="REML",test="t",dfs="contain",
		                           random=list(~1|species_ID,~1|phylogeny,~1|es_ID),R =list(phylogeny = phylo_matrix),data=pottier)
		mod.habitat_hom

###################
## Figure 1
###################

	# Figure for het model
	#dev.size() just find right sizing
	fig1a <-orchard_plot(mod.habitat_het, data = pottier, group = "species_ID", mod = "habitat", xlab = "Developmental Acclimation Response Ratio (dARR)", angle = 45)
	#ggsave(filename = "./figures/fig1A.pdf", width = 5.137255, height = 4.086274)

	# Homogen. Variance
	fig1b <- orchard_plot(mod.habitat_hom, data = pottier, group = "species_ID", mod = "habitat", xlab = "Developmental Acclimation Response Ratio (dARR)", angle = 45)
	#ggsave(filename = "./figures/fig1B.pdf", width = 5.137255, height = 4.086274)

	# Write the mod_table
	form_tab_a <- round_df(mod_results(mod.habitat_hom, data = pottier, group = "species_ID", mod = "habitat")$mod_table, 2) 
	form_tab_b <- round_df(mod_results(mod.habitat_het, data = pottier, group = "species_ID", mod = "habitat")$mod_table, 2) 

	fig1_tab_a <- flextable::flextable(form_tab_a) %>% flextable::mk_par(part = "header", value = flextable::as_paragraph(c("Name", "Mean", "L 95% CI", "U 95% CI", "L 95% PI", "U 95% PI"))) %>% autofit()
	#save_as_image(fig1_tab_a, "./figures/fig1a_tab.png")

	fig1_tab_b <- flextable::flextable(form_tab_b) %>% flextable::mk_par(part = "header", value = flextable::as_paragraph(c("Name", "Mean", "L 95% CI", "U 95% CI", "L 95% PI", "U 95% PI"))) %>% autofit()
	#save_as_image(fig1_tab_b, "./figures/fig1b_tab.png")

###################
## Figure 2
###################
	
	 prop <- mod_results(mod.habitat_hom, data = pottier, group = "species_ID", weights = "prop")
	equal <- mod_results(mod.habitat_hom, data = pottier, group = "species_ID", weights = "equal")

	fig2a <- orchard_plot(mod.habitat_hom, data = pottier, group = "species_ID", xlab = "Developmental Acclimation Response Ratio (dARR)", angle = 45, weights = "prop") + ylim(c(-1,1)) + annotate("text", x = 1.5, y = 0, label = TeX(paste0("$\\mu$ = ", round(prop$mod_table[1,2], 2), ", 95% CI = ", round(prop$mod_table[1,3], 2), " to ", round(prop$mod_table[1,4], 2)))) + #annotate("text", x = 1.5, y = 0, label = TeX("\\textbf{Proportional}")) + 
		ggtitle(TeX("\\textbf{Proportional}")) + 
		theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank())


	fig2b <- orchard_plot(mod.habitat_hom, data = pottier, group = "species_ID", xlab = "Developmental Acclimation Response Ratio (dARR)", angle = 45, weights = "equal") + ylim(c(-1,1)) + annotate("text", x = 1.5, y = -0.1, label = TeX(paste0("$\\mu$ = ", round(equal$mod_table[1,2], 2), ", 95% CI = ", round(equal$mod_table[1,3], 2), " to ", round(equal$mod_table[1,4], 2)))) + #annotate("text", x = 1.5, y = 0, label = TeX("\\textbf{Equal}")) 
		ggtitle(TeX("\\textbf{Equal}")) + 
		theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank())

	(fig2a | fig2b) + plot_annotation(tag_levels = "A", tag_suffix = ")") & theme(plot.tag = element_text(family = "Palatino", size = 24, face = "bold"))
	ggsave(filename = "./figures/fig2.pdf", width = 8.839216, height = 3.929412)

###################
## Figure 3
###################

	# Load the dataset that comes with orchaRd
	    data(fish)

	# Subset data for demonstration purposes.
	    warm_dat <- fish

	 # Fit the metaregerssion model
	model_fish <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, 
	                         mods = ~ experimental_design + trait.type + deg_dif + treat_end_days,
	                         method = "REML", test = "t", 
	                         random = list(~1 | group_ID, ~1 + trait.type| es_ID), 
	                                rho = 0, struc = "HCS", 
	                         data = warm_dat, 
	                         control=list(optimizer="optim", optmethod="Nelder-Mead"))

	fig3a <- orchaRd::orchard_plot(model_fish, group = "group_ID", mod = "trait.type", weights = "prop", data = warm_dat, xlab = "log Response Ratio (lnRR)", angle = 45, g = FALSE, legend.pos = "top.left", condition.lab = "Temperature Difference") + theme(legend.direction = "vertical", panel.grid = element_blank())
	#ggsave(filename = "./figures/fig4a.pdf", width = 4.337255, height = 4.258823)
	 
	fig3b <- orchaRd::orchard_plot(model_fish, group = "group_ID", mod = "trait.type", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif", weights = "prop", data = warm_dat, xlab = "log Response Ratio (lnRR)", angle = 45, g = FALSE, legend.pos = "top.left", condition.lab = "Temperature Difference") + theme(legend.direction = "vertical", panel.grid = element_blank(), axis.text.y = element_blank())
	#ggsave(filename = "./figures/fig4b.pdf", width = 4.337255, height = 4.258823)

	(fig3a | fig3b)  + plot_annotation(tag_levels = "A", tag_suffix = ")") & theme(plot.tag = element_text(family = "Palatino", size = 24, face = "bold"))
	ggsave(filename = "./figures/fig3.pdf", width = 10.22745, height = 4.92549)
 
###################
## Figure 4
###################
	data(lim)
	lim[, "year"] <- as.numeric(lim$year)
	lim$vi<- 1/(lim$N - 3)
	lim$si <- sqrt(lim$vi)

	model<-metafor::rma.mv(yi=yi, V=vi, mods= ~Environment*year,
	                       random=list(~1|Article,~1|Datapoint), data=na.omit(lim))

	lim2 <- lim %>% filter(!Reproduction == "?")
	lim2 <- lim2 %>% mutate(rep_prop = paste(Environment, Propagule, sep = "."))

	model2<-metafor::rma.mv(yi=yi, V=vi, mods= ~rep_prop,
	                       random=list(~1|Article,~1|Datapoint), data=na.omit(lim2))

	model3<-metafor::rma.mv(yi=yi, V=vi, mods= ~si*year,
	                       random=list(~1|Article,~1|Datapoint), data=na.omit(lim))


	fig4a <- orchaRd::orchard_plot(model2, group = "Article", mod = "rep_prop", legend.pos = "top.left", xlab = "Fisher's Z-transformed Correlation Coefficient (Zr)", data = lim2, angle = 45) +  ggtitle('Categorical x Categorical') + theme(plot.title = element_text(hjust = 0.5))


	lim_bubble <- orchaRd::mod_results(model, mod = "year", group = "Article",
	                    data = lim, weights = "prop", by = "Environment")

	fig4b <- orchaRd::bubble_plot(lim_bubble, data = lim, group = "Article", mod = "year", xlab = "Year", legend.pos = "top.left", ylab = "Fisher's Z-transformed Correlation Coefficient (Zr)") + ggtitle('Continuious x Categorical') + theme(plot.title = element_text(hjust = 0.5))

	lim_bubble2 <- orchaRd::mod_results(model3, mod = "si", group = "Article", at=list(year = c(1972, 2012)), data = lim, weights = "prop", by = "year")

	fig4c <- orchaRd::bubble_plot(lim_bubble2, data = lim, group = "Article", mod = "vi", legend.pos = "top.left", ylab = "Fisher's Z-transformed Correlation Coefficient (Zr)", xlab = "Sampling Standard Error", k = FALSE, g = FALSE) + ggtitle('Continuious x Continuious')+ theme(plot.title = element_text(hjust = 0.5))


	(fig4a | fig4b | fig4c) + plot_layout(widths = unit(c(10,12,7), c("cm", "cm", "cm"))) + plot_annotation(tag_levels = "A", tag_suffix = ")") & theme(plot.tag = element_text(family = "Palatino", size = 24, face = "bold"))
	ggsave(filename = "./figures/fig4.pdf", width = 15.968627, height = 6.509804)