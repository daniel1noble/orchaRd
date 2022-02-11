
#' @title eklof
#' @description Eklof et al. (2012) evaluated the effects of predation on benthic invertebrate communities. Using the log response ratio they quantified differences in abundance and/or biomass of gastropods and Amphipods in groups with and without predation in an experimental setting. Below we describe the variables used in our examples and analyses.
#' @format A data frame :
#' \describe{
#'   \item{ExptID}{Experiment / effect size ID}
#'   \item{First.author}{First author of publication}
#' 	 \item{Publication.year}{Year of publication}
#' 	 \item{Journal}{Journal published within}
#' 	 \item{Research.group}{Research group that data arose from}
#' 	 \item{Grazer.type}{Type pf grazer}
#' 	 \item{mean_treatment}{Mean abundance/biomass of invertebrates in treatment group}
#' 	 \item{SD_treatment}{Standard deviation in abundance/biomass of invertebrates in treatment group}
#' 	 \item{N_treatment}{Sample size of invertebrate sample in treatment group}
#' 	 \item{mean_control}{Mean abundance/biomass of invertebrates in control group}
#' 	 \item{SD_control}{Standard deviation in abundance/biomass of invertebrates in control group}
#' 	 \item{N_control}{Sample size of invertebrate sample in control group}
#'   ...
#' }
#' @references Eklof J.S., Alsterberg C., Havenhand J.N., Sundback K., Wood H.L., Gamfeldt L. 2012. Experimental climate change weakens the insurance effect of biodiversity. Ecology Letters, 15:864-872. https://doi.org/10.1111/j.1461-0248.2012.01810.x
#' @name eklof
#' @docType data
NULL

#' @title english
#' @description English and Uller (2016) performed a systematic review and meta-analysis on the effects of early life dietary restriction (a reduction in a major component of the diet without malnutrition; e.g. caloric restriction) on average at death, using the standardised mean difference (often called *d*). In a subsequent publication, Senior et al. (2017) analysed this dataset for effects of dietary-restriction on among-individual variation in the age at death using the log coefficient of variation ratio. A major prediction in both English & Uller (2016) and Senior et al. (2017) was that the type pf manipulation, whether the study manipulated quality of food versus the quantity of food would be important.
#' @format A data frame :
#' \describe{
#'   \item{StudyNo}{Study ID}
#'   \item{EffectID}{Effect size ID}
#'   \item{Author}{First author of study}
#' 	  \item{Year}{Year of study publication}
#'    \item{Journal}{Research journal study was published}
#' 	  \item{Species}{Common name of species}
#'    \item{Phylum}{Phylum of species}
#'    \item{ExptLifeStage}{Life stage when manipulation was undertaken}
#' 		\item{ManipType}{Type of food manipulation}
#' 		\item{CatchUp}{Whether species exhibits catchup growth}
#' 		\item{Sex}{Sex of organisms in sample}
#' 		\item{AdultDiet}{Diet adults were provided and whether it was restricted or control}
#' 		\item{NStartControl}{Sample size of the control group}
#' 		\item{NStartExpt}{Sample size of the treatment group}
#' 		\item{MeanC}{Mean of the control group}
#' 		\item{MeanE}{Mean of the treatment/experimental group}
#' 		\item{SD_C}{Standard deviation of the control group}
#' 		\item{SD_E}{Standard deviation of the treatment/experimental group}
#'   ...
#' }
#' @references English S, Uller T. 2016. Does early-life diet affect longevity? A meta-analysis across experimental studies. Biology Letters, 12: http://doi:10.1098/rsbl.2016.0291
#' @name english
#' @docType data
NULL

#' @title lim
#' @description Lim et al. (2014) meta-analysed the strength of correlation between maternal and offspring size within-species, across a very wide range of taxa. They found, that typically, there is a moderate positive correlation between maternal size and offspring size within species (i.e. larger mothers have larger offspring). Below we describe the variables use in our examples and or analyses.
#' @format A data frame :
#' \describe{
#' \item{Amniotes}{Amniote group}
#' \item{Article}{Study or research paper that data were collected}
#' \item{Author}{Authors of article}
#' \item{Class}{Class of organism}
#' \item{Common.Name}{Common name of species}
#' \item{Datapoint}{Observation level ID that identified each unique datapoint}
#' \item{Environment}{Environment organism is found, wild or captive}
#' \item{Family}{Family of species}
#' \item{Genus}{Genus of species}
#' \item{Maternal}{Mother length, weight or mass}
#' \item{N}{Sample size used to estaimted correlation}
#' \item{Offspring}{}
#' \item{Order}{Order of species}
#' \item{Phylum}{Phylum of the species}
#' \item{Propagule}{}
#' \item{RU}{}
#' \item{Reproduction}{Reproductive mode of species}
#' \item{animal}{Species name use for phylogenetic reconstruction}
#' \item{year}{year of study}
#' \item{yi}{Effect size correlation coefficient between maternal and offspring size within-species}
#'   ...
#' }
#' @references Lim J.N., Senior A.M., Nakagawa S. 2014. Heterogeneity in individual quality and reproductive trade-offs within species. Evolution. 68(8):2306-18. doi: 10.1111/evo.12446
#' @name lim
#' @docType data
NULL

#' @title fish
#' @description Fish data example.
#' @format A data frame
#' \describe{
#' \item{es_ID}{Effect size ID}
#' \item{species_ID}{Species ID}
#' \item{paper_ID}{Paper ID}
#' \item{group_ID}{Group ID}
#' \item{experimental_design}{Experimental design type}
#' \item{deg_dif}{Temperature difference between treatments}
#' \item{treat_end_days}{How long temperatures were applied}
#' \item{n_control}{Sample size for control}
#' \item{n_treat}{Sample size for treatment}
#' \item{mean_control}{Mean for control}
#' \item{mean_treat}{Mean for treatment}
#' \item{sd_control}{SD for control}
#' \item{sd_treat}{SD for treatment}
#' \item{lnrr}{Log Response ratio effect size}
#' \item{lnrr_vi}{Sampling variance for lnRR}
#' }
#' @references Lim J.N., Senior A.M., Nakagawa S. 2014. Heterogeneity in individual quality and reproductive trade-offs within species. Evolution. 68(8):2306-18. doi: 10.1111/evo.12446
#' @name fish
#' @docType data
NULL
