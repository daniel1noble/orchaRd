
#' @title eklof
#' @description Eklof et al. (2012) evaluated the effects of predation on benthic invertebrate communities. Using the log response ratio they quantified differences in abundance and/or biomass of gastropods and Amphipods in groups with and without predation in an experimental setting. Below we describe the variables used in our examples and analyses.
#' @format A data frame :
#' \describe{
#'   \item{ExptID}{Experiment / effect size ID}
#'   \item{First.author}{First author of publication}
#' 	 \item{Publication.year}{Year of publication}
#' 	 \item{Journal}{Journal published within}
#' 	 \item{Research.group}{Research group that data arose from}
#' 	 \item{Grazer.type}{Type of grazer}
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
#' @description English and Uller (2016) performed a meta-analysis on the effects of early life dietary restriction (a reduction in a major component of the diet without malnutrition; e.g. caloric restriction) on average age at death, using the standardised mean difference (often called *d*). In a subsequent publication, Senior et al. (2017) analysed this dataset for effects of dietary-restriction on among-individual variation in the age at death using the log coefficient of variation ratio. A major prediction in both English & Uller (2016) and Senior et al. (2017) was that the type of manipulation, whether the study manipulated quality of food versus the quantity of food, would be important.
#' @format A data frame :
#' \describe{
#'   \item{StudyNo}{Study ID}
#'   \item{EffectID}{Effect size ID}
#'   \item{Author}{First author of study}
#' 	  \item{Year}{Year of study publication}
#'    \item{Journal}{Research journal where study was published}
#' 	  \item{Species}{Common name of species}
#'    \item{Phylum}{Phylum of species}
#'    \item{ExptLifeStage}{Life stage when manipulation was undertaken}
#' 		\item{ManipType}{Type of food manipulation}
#' 		\item{CatchUp}{Whether species exhibits catchup growth}
#' 		\item{Sex}{Sex of organisms in sample}
#' 		\item{AdultDiet}{Diet adults were provided and whether it was restricted (treatment) or the control}
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
#' @description Lim et al. (2014) meta-analysed the strength of correlation between maternal and offspring size within-species, across a very wide range of taxa. They found that typically there is a moderate positive correlation between maternal size and offspring size within species (i.e. larger mothers have larger offspring). Below we describe the variables use in our examples and analyses.
#' @format A data frame :
#' \describe{
#' \item{Amniotes}{Amniote group}
#' \item{Article}{Study or research paper that data were collected from}
#' \item{Author}{Authors of article}
#' \item{Class}{Class of organism}
#' \item{Common.Name}{Common name of species}
#' \item{Datapoint}{Observation level ID that identified each unique datapoint}
#' \item{Environment}{Environment organism is found, wild or captive}
#' \item{Family}{Family of species}
#' \item{Genus}{Genus of species}
#' \item{Maternal}{Mother length, or mass}
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
#' \item{deg_dif}{Temperature difference between treatments}
#' \item{es_ID}{Effect size ID}
#' \item{experimental_design}{Experimental design type}
#' \item{group_ID}{Group ID}
#' \item{lnrr}{Log Response ratio effect size}
#' \item{lnrr_vi}{Sampling variance for lnRR}
#' \item{mean_control}{Mean for control}
#' \item{mean_treat}{Mean for treatment}
#' \item{n_control}{Sample size for control}
#' \item{n_treat}{Sample size for treatment}
#' \item{paper_ID}{Paper ID}
#' \item{sd_control}{SD for control}
#' \item{sd_treat}{SD for treatment}
#' \item{species_ID}{Species ID}
#' \item{trait.type}{Type of trait effect size was calculated from}
#' \item{treat_end_days}{Duration of time developmental treatmemt was applied}
#' }
#' @references Lim J.N., Senior A.M., Nakagawa S. 2014. Heterogeneity in individual quality and reproductive trade-offs within species. Evolution. 68(8):2306-18. doi: 10.1111/evo.12446
#' @name fish
#' @docType data
NULL

#' @title pottier
#' @description Pottier et al. (2021) meta-analysis on sex-differences in acclimation
#' @format A data frame :
#' \describe{
#'   \item{initials}{}
#'   \item{es_ID}{}
#'   \item{study_ID}{}
#'   \item{species_ID}{}
#'   \item{population_ID}{}
#'   \item{family_ID}{}
#'   \item{shared_trt_ID}{}
#'   \item{cohort_ID}{}
#'   \item{note_ID}{}
#'   \item{data_source}{}
#'   \item{data_url}{}
#'   \item{fig_file_name}{}
#'   \item{data_type}{}
#'   \item{data_file_name}{}
#'   \item{peer.reviewed}{}
#'   \item{ref}{}
#'   \item{title}{}
#'   \item{pub_year}{}
#'   \item{journal}{}
#'   \item{thesis_chapter}{}
#'   \item{doi}{}
#'   \item{citation}{}
#'   \item{phylum}{}
#'   \item{class}{}
#'   \item{order}{}
#'   \item{family}{}
#'   \item{genus}{}
#'   \item{species}{}
#'   \item{genus_species}{}
#'   \item{habitat}{}
#'   \item{taxonomic_group}{}
#'   \item{reproduction_mode}{}
#'   \item{life_stage_manip}{}
#'   \item{life_stage_tested}{}
#'   \item{brought_common_temp}{}
#'   \item{mobility_life_stage_manip}{}
#'   \item{time_common_temp}{}
#'   \item{common_temp}{}
#'   \item{exp_design}{}
#'   \item{origin_hatching}{}
#'   \item{latitude}{}
#'   \item{longitude}{}
#'   \item{elevation}{}
#'   \item{season}{}
#'   \item{year}{}
#'   \item{body_length}{}
#'   \item{body_mass}{}
#'   \item{age_tested}{}
#'   \item{sex}{}
#'   \item{housing_temp}{}
#'   \item{incubation_independent}{}
#'   \item{metric}{}
#'   \item{endpoint}{}
#'   \item{acc_temp_low}{}
#'   \item{acc_temp_high}{}
#'   \item{acc_temp_var}{}
#'   \item{is_acc_temp_fluctuating}{}
#'   \item{acc_duration}{}
#'   \item{ramping}{}
#'   \item{set_time}{}
#'   \item{n_test_temp}{}
#'   \item{n_replicates_per_temp}{}
#'   \item{n_animals_per_replicate}{}
#'   \item{humidity}{}
#'   \item{oxygen}{}
#'   \item{salinity}{}
#'   \item{pH}{}
#'   \item{photoperiod}{}
#'   \item{gravidity}{}
#'   \item{starved}{}
#'   \item{minor_concerns}{}
#'   \item{major_concerns}{}
#'   \item{notes_moderators}{}
#'   \item{mean_HT_low}{}
#'   \item{sd_HT_low}{}
#'   \item{n_HT_low}{}
#'   \item{mean_HT_high}{}
#'   \item{sd_HT_high}{}
#'   \item{n_HT_high}{}
#'   \item{error_type}{}
#'   \item{notes_es}{}
#'   \item{exclude}{}
#'   \item{n_trt}{}
#'   \item{n_cohort}{}
#'   \item{within_study_mean_low}{}
#'   \item{within_study_mean_high}{}
#'   \item{within_study_sd_low}{}
#'   \item{within_study_sd_high}{}
#'   \item{between_study_mean_low}{}
#'   \item{between_study_mean_high}{}
#'   \item{between_study_sd_low}{}
#'   \item{between_study_sd_high}{}
#'   \item{imputed}{}
#'   \item{imputed_sd_low}{}
#'   \item{imputed_sd_high}{}
#'   \item{dARR}{}
#'   \item{Var_dARR}{}
#'   \item{precision}{}
#'   \item{is_concern}{}
#'   \item{search_string}{}
#'   \item{unique_name}{}
#'   \item{ott_id}{}
#'   \item{phylogeny}{}
#'   ...
#' }
#' @references Pottier et al. 2022. Functional Ecology
#' @name pottier
#' @docType data
NULL

#' @title phylo_matrix
#' @description Phylogenetic matrix used in Pottier et al. (2021) meta-analysis on sex-differences in acclimation
#' @format A phylogenetic correlation matrix :
#' @references Pottier et al. 2022. Functional Ecology
#' @name phylo_matrix
#' @docType data
NULL

#' @title VCV_dARR
#' @description Sampling error matrix used in Pottier et al. (2021) meta-analysis on sex-differences in acclimation
#' @format A sampling variance matrix :
#' @references Pottier et al. 2022. Functional Ecology
#' @name VCV_dARR
#' @docType data
NULL
