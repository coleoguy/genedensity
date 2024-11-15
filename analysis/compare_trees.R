

library(ape)
library(phangorn)
library(dendextend)
source("../analysis/functions.R")
dat <- read.csv("../results/vertebrates/unparsed.csv")
dat <- dat[dat$clade == "Mammalia", ]
timetree <- read.tree("../data/vertebrates/chordates_species.nwk")
goodtree <- read.nexus("../data/vertebrates/mammal.nex")

# format and prune good tree
sp <- unique(dat$species)
spf <- sub("^([^_]*_[^_]*)_.*", "\\1", gsub(" ", "_", sp))
intersect.timetree <- intersect(timetree$tip.label, spf)
pruned.timetree <- keep.tip(timetree, intersect.timetree)
spn <- sp[match(spf[match(pruned.timetree$tip.label, spf)], spf)]
pruned.timetree$tip.label <- spn

# format and prune timetree tree
spf <- tolower(spf)
intersect.goodtree <- intersect(goodtree$tip.label, spf)
pruned.goodtree <- keep.tip(goodtree, intersect.goodtree)
spn <- sp[match(spf[match(pruned.goodtree$tip.label, spf)], spf)]
pruned.goodtree$tip.label <- spn

# prune by intersection
intersect.both <- intersect(pruned.timetree$tip.label, pruned.goodtree$tip.label)
pruned.goodtree <- keep.tip(pruned.goodtree, intersect.both)
pruned.timetree <- keep.tip(pruned.timetree, intersect.both)

#compare

#summary
comparePhylo(unroot(pruned.goodtree), unroot(pruned.timetree))
#distances
treedist(unroot(pruned.goodtree), unroot(pruned.timetree))
#correlation
cor(pruned.goodtree$edge.length, pruned.timetree$edge.length)
#agreement
agree.tree <- consensus(pruned.goodtree, pruned.timetree, p = 1.0)

#write
write.tree(pruned.timetree, file = paste0("../data/vertebrates/timetree.nwk"))
write.tree(pruned.goodtree, file = paste0("../data/vertebrates/goodtree.nwk"))
write.tree(agree.tree, file = paste0("../data/vertebrates/agree.nwk"))







