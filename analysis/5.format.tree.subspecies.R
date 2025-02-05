



# load stuff in
library(ape)
library(dispRity)
dat <- read.csv("../data/data.csv")
tree <- read.tree("../data/chordates.species.nwk")

# format and prune tree
sp <- unique(dat$species)
sp.underscore <- sub("^([^_]*_[^_]*)_.*", "\\1", gsub(" ", "_", sp))
sp.intersect <- intersect(tree$tip.label, sp.underscore)
pruned.tree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% sp.intersect)])
sp.match.subsp <- sp[match(sp.underscore[match(pruned.tree$tip.label, sp.underscore)], 
                           sp.underscore)]
pruned.tree$tip.label <- sp.match.subsp

write.tree(pruned.tree, file = paste0("../data/formatted.tree.nwk"))







