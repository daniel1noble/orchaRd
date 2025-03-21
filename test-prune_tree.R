library(devtools)
library(metafor)
library(ape)
library(rotl)

load_all()

spp_names_1 <- c("Festuca arundinacea",
               "Lolium perenne",
               "Solanum tuberosum",
               "Nassella neesiana")

# get phylo matrix
taxa1 <- tnrs_match_names(spp_names_1, context_name = "Land plants")
tree1 <- tol_induced_subtree(taxa1$ott_id)
tree1$tip.label <- strip_ott_ids(tree1$tip.label, remove_underscores = TRUE)

# Compute phylo matrix
tree1 <- ape::compute.brlen(tree1)
phylo_matrix_1 <- ape::vcv(tree1, cor = TRUE)

phylo_matrix_1

# Species 2

spp_names_2 <- c("Festuca arundinacea",
               "Lolium perenne",
               "Nassella neesiana")

# get phylo matrix
taxa2 <- tnrs_match_names(spp_names_2, context_name = "Land plants")
tree2 <- tol_induced_subtree(taxa2$ott_id)
tree2$tip.label <- strip_ott_ids(tree2$tip.label, remove_underscores = TRUE)

# Compute phylo matrix
tree2 <- ape::compute.brlen(tree2)
phylo_matrix_2 <- ape::vcv(tree2, cor = TRUE)


phylo_matrix_2

# What about .create_tmp_phylo_matrix?

dat <- data.frame(spp_names = spp_names_2, stringsAsFactors = FALSE)
phylo_args <- list(tree = tree1, species_colname = "spp_names")

phylo_matrix_3 <- .create_tmp_phylo_matrix(dat, phylo_args)

phylo_matrix_3

identical(phylo_matrix_3, phylo_matrix_1)

identical(phylo_matrix_3, phylo_matrix_2)

