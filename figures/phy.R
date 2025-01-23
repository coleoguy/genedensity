


# load stuff in
library(phytools)
library(viridis)
dat <- read.csv("../results/parsed.csv")
dat <- dat[dat$thrs == 0.8, ]
dat <- dat[!duplicated(dat$species), ]
dat <- dat[, c("species", "clade", "rsq")]
tree <- read.tree("../data/formatted_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
int <- intersect(tree$tip.label, dat$species)
pruned.tree <- keep.tip(tree, int)
pruned.tree <- ladderize(pruned.tree)
dat <- dat[dat$species %in% int, ]
dat <- dat[match(pruned.tree$tip.label, dat$species), ]

rsq <- setNames(dat$rsq, dat$species)
p <- contMap(pruned.tree, rsq, plot = FALSE)
p <- setMap(p, viridis(n=300))



plot(p, 
     res = 300,
     lwd = c(1.5, 5),
     fsize = c(0.3, 0.5), 
     outline = FALSE,
     sig = 2,
     xlim = c(0, 115), 
     ylim = c(-150, 450),  #450
     direction = "downwards", 
     ftype = "off", 
     plot = TRUE)


#plot(p, 
#     res = 300,
#     fsize = c(0.00000000001, 0.5), 
#     lwd = c(1.5, 5),
#     outline = FALSE,
#     sig = 2,
#     type = "fan",
#     xlim = c(-500, 500), 
#     ylim = c(-500, 500), 
#     plot = TRUE)