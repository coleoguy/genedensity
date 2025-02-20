
library(phytools)
library(geiger)
library(caper)
library(geiger)
library(mvMORPH)
library(bayou)

# subset results
dat <- read.csv("../results/parsed.csv")
dat <- dat[!is.na(dat$chromnum.1n) & !duplicated(dat$species), ]
dat <- na.omit(dat[, c("species", "clade", "rsq")])


# prune tree
tree <- read.tree("../data/formatted.tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
int <- intersect(dat$species, tree$tip.label)
dat <- dat[dat$species %in% int, ]
pruned.tree <- keep.tip(tree, int)

rsq.vec <- setNames(dat$rsq, dat$species)
rsq.mtx <- as.matrix(rsq.vec)

fivenum(rsq.vec)
var(rsq.vec)



prior <- make.prior(tree = pruned.tree, 
                    dists = list(alpha = "exp", sigma2 = "exp", theta = "dnorm"), 
                    param = list(alpha = 1, sigma2 = 1, theta = c(0, 1)))










bm.nodrift <- mvBM(pruned.tree, rsq.mtx, model = "BM1")
bm.multirate <- mvBM(pruned.tree, rsq.mtx, model = "BMM")
bm.timedependent <- evo.model(rsq.vec, pruned.tree, model = "BMrate")












l <- fitContinuous(phy = pruned.tree, 
                   dat = rsq,
                   model = c("BM","OU","EB","rate_trend","lambda","kappa","delta","mean_trend","white"), 
                   niter = 10000)






bm <- fitContinuous(tree, rsq, model = "BM", niter = 10000)
ou <- fitContinuous(tree, rsq, model = "OU", niter = 10000)
eb <- fitContinuous(tree, rsq, model = "EB", niter = 10000)
rt <- fitContinuous(tree, rsq, model = "rate_trend", niter = 10000)
lambda <- fitContinuous(tree, rsq, model = "lambda", niter = 10000)
kappa <- fitContinuous(tree, rsq, model = "kappa", niter = 10000)
delta <- fitContinuous(tree, rsq, model = "delta", niter = 10000)
mt <- fitContinuous(tree, rsq, model = "mean_trend", niter = 10000)
white <- fitContinuous(tree, rsq, model = "white", niter = 10000)









