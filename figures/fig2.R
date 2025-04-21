

# r2 beeswarm
library(beeswarm)
dat <- read.csv("../results/rsq.csv")
dat <- dat[dat$clade %in% c("Mammalia", "Sauropsida", "Actinopterygii"), ]
dat$clade <- factor(dat$clade, levels = c("Mammalia", "Actinopterygii", "Sauropsida"))
map <- c("Mammalia" = "#d95f02", "Actinopterygii" = "#7570b3", "Sauropsida" = "#1b9e77")
cols <- map[dat$clade]
nmam <- nrow(dat[dat$clade == "Mammalia", ])
nfish <- nrow(dat[dat$clade == "Actinopterygii", ])
nrep <- nrow(dat[dat$clade == "Sauropsida", ])
beeswarm(rsq ~ clade,
xlab = NA,
ylab = "R-squared",
data = dat,
pwcol = cols,
pch = 16,
spacing = 1.2,
labels = NA)
text(1, 0.06, "Mammals", xpd = NA, adj = c(0.5, 0.5))
text(2, 0.06, "Ray-finned fish", xpd = NA, adj = c(0.5, 0.5))
text(3, 0.06, "Reptiles", xpd = NA, adj = c(0.5, 0.5))
text(1, -0.02, paste0("n=", nmam), xpd = NA, adj = c(0.5, 0.5))
text(2, -0.02, paste0("n=", nfish), xpd = NA, adj = c(0.5, 0.5))
text(3, -0.02, paste0("n=", nrep), xpd = NA, adj = c(0.5, 0.5))


# phylogenetic anova
library(phytools)
tree <- read.tree("../data/formatted.tree.nwk")
int <- intersect(dat$species, tree$tip.label)
pruned.tree <- keep.tip(tree, int)
sig <- phylosig(pruned.tree,
                setNames(dat$rsq, dat$species),
                method = "lambda",
                test = TRUE,
                niter = 10000)[[4]]
if (sig < 0.05) {
  x <- setNames(dat$clade, dat$species)
  y <- setNames(dat$rsq, dat$species)
  phylANOVA(pruned.tree, x, y, nsim = 10000, posthoc = TRUE)
}

