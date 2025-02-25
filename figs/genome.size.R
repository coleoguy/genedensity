

# largest repeats in each clade
library(phytools)
library(caper)
dat <- read.csv("../results/parsed.csv")
dat <- dat[!duplicated(dat$species) & 
             !is.na(dat$chromnum.1n) & 
             !is.na(dat$est.gnsz.Mb), ]
dat <- na.omit(dat[, c("species", "clade", "rsq", "est.gnsz.Mb")])

tree <- read.tree("../data/formatted.tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)

int <- intersect(dat$species, tree$tip.label)
pruned.tree <- keep.tip(tree, int)
dat <- dat[dat$species %in% int, ]

cd <- comparative.data(pruned.tree, dat, names.col = "species", vcv = T)
model <- pgls(rsq ~ est.gnsz.Mb, data = cd)

cols <- c("#7570b3", "#d95f02", "#e7298a", "#1b9e77")[as.factor(dat$clade)]
plot(dat$est.gnsz.Mb / 100, 
     dat$rsq, 
     pch = 16, 
     ylim = c(0, 1.1), 
     cex = 0.8, 
     col = cols)
abline(model, lwd = 2, )

