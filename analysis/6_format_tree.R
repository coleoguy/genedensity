



# load stuff in
library(ape)
library(dispRity)
source("../analysis/functions.R")
dat <- read.csv("../data/data.csv")
tree <- read.tree("../data/chordates_species.nwk")

# format and prune tree
sp <- unique(dat$species)
spf <- sub("^([^_]*_[^_]*)_.*", "\\1", gsub(" ", "_", sp))
sp.intersect <- intersect(tree$tip.label, spf)
pruned.tree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% sp.intersect)])
spn <- sp[match(spf[match(pruned.tree$tip.label, spf)], spf)]
pruned.tree$tip.label <- spn

write.tree(pruned.tree, file = paste0("../data/formatted_tree.nwk"))







