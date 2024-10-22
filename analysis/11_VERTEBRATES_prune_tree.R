



# load stuff in
library(ape)
source("../analysis/functions.R")
dat <- read.csv("../results/vertebrates/final_results.csv")
tree <- read.tree("../data/vertebrates/chordates_species.nwk")

# prune and format tree
sp <- dat$species
spf <- sub("^([^_]*_[^_]*)_.*", "\\1", gsub(" ", "_", sp))
sp.intersect <- intersect(tree$tip.label, spf)
pruned.tree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% sp.intersect)])
spn <- sp[match(spf[match(pruned.tree$tip.label, spf)], spf)]
pruned.tree$tip.label <- spn
write.tree(pruned.tree, file = paste0("../data/vertebrates/pruned_tree.nwk"))







