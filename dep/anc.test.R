



library(phytools)
tree <- read.tree("../data/formatted_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
dat <- read.csv("../results/parsed.csv")
dat <- dat[dat$thrs == 0.8, ]
dat <- dat[, c("species", "rsq", "clade"), ]
dat <- daht[!duplicated(dat$species), ]
int <- intersect(tree$tip.label, dat$species)
dat <- dat[dat$species %in% int, ]
pruned.tree <- keep.tip(tree, int)
dat <- dat[match(pruned.tree$tip.label, dat$species), ]
dat <- setNames(dat$rsq, dat$species)

ml <- fastAnc(pruned.tree, dat, CI = TRUE)
bayes <- anc.Bayes(pruned.tree, dat, ngen=10000)
