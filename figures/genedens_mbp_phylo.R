

# load stuff in
library(phytools)
dat <- read.csv("../results/vertebrates/final_results.csv")
tree <- read.tree("../data/vertebrates/pruned_tree.nwk")

# subset data for beta
dat <- na.omit(dat[, c("species", "beta", "clade")])
dat$beta <- dat$beta * 1000000
dat <- dat[dat$beta > 0, ]
tree$tip.label <- gsub("_", " ", tree$tip.label)
sp.intersect <- intersect(dat$species, tree$tip.label)
dat <- dat[dat$species %in% sp.intersect, ]
pruned.tree <- keep.tip(tree, sp.intersect)
b <- setNames(dat$beta, dat$species)
map <- contMap(pruned.tree, b, plot = FALSE)

# which orders to label on the phylogeny
clades <- c("Actinopterygii", "Mammalia", "Sauria")

# function to add clade labels
cladeLabel <- function(i, dat) {
  clade.sp <- dat[dat$clade == i, ]$species
  node <- findMRCA(pruned.tree, tips = clade.sp)
  arc.cladelabels(node = node, 
                  text = i, 
                  ln.offset = 1.03, 
                  lab.offset = 1.11,
                  mark.node = F) 
}

# plot total chromosome number
plot(map, 
     res = 600,
     fsize = c(0.000000000001, 0.8), 
     lwd = c(1.7, 5),
     outline = FALSE,
     sig = 2,
     type = "fan")

# apply function
sapply(clades, cladeLabel, dat = dat)



