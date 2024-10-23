



vert.invert <- "invertebrates"

# load library
library(ape)
# source functions
source("functions.R")
# load in tree
tree <- read.tree(paste0("../data/", vert.invert, "/pruned_tree.nwk"))
# load in results
final.results <- read.csv(paste0("../results/", vert.invert, "/final_results.csv"))
# rename results as species
rownames(final.results) <- final.results$species
# remove non-numeric data
dat <- na.omit(final.results[sapply(final.results, is.numeric)])
# vector of species in data frame
species <- gsub(" ", "_", rownames(dat))
# re-add species
dat$species <- species
# intersection between tip label and species
species <- intersect(tree$tip.label, species)
# reorder and maybe subset dataframe
dat <- dat[match(species, dat$species), ]
# species
species <- dat$species
# find remove species column
dat <- dat[sapply(dat, is.numeric)]
# prune tree based on species
pruned.tree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% intersect(tree$tip.label, species))])
# calculate pics
pics <- as.data.frame(sapply(dat, calcPic, species = species, tree = pruned.tree))
# write pics
write.csv(pics, paste0("../results/", vert.invert, "/final_results_pics.csv"), row.names = FALSE)
# pca
pca_result <- prcomp(pics, center = TRUE, scale. = TRUE)
# get loading scores
loadings <- as.data.frame(pca_result$rotation)
# write pic loadings
write.csv(loadings, paste0("../results/", vert.invert, "/final_results_pic_loadings.csv"))


