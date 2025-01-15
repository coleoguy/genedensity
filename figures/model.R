
# mammals; all repeats

library(phytools)
library(caper)
dat <- read.csv("../results/parsed.csv")
dat <- dat[dat$thrs == 0.9, ]
dat <- dat[!duplicated(dat$species), ]
dat$median.trans <- 1 - (dat[["total.rep.median"]]/70)
dat <- na.omit(dat[, c("species", "rsq", "clade", "median.trans")])
dat <- dat[dat$clade == "Mammalia", ]
dat <- dat[dat$species != "Callithrix jacchus", ]
plot(dat$median.trans, dat$rsq)
model <- glm(rsq ~ median.trans, data = dat)
res <- setNames(resid(model), dat$species)

tree <- read.tree("../data/formatted_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
int <- intersect(tree$tip.label, dat$species)
pruned.tree <- keep.tip(tree, int)

if (phylosig(pruned.tree, res, method = "lambda", test = TRUE)[[4]] < 0.05) {
  cd <- comparative.data(pruned.tree, dat, names.col = "species")
  model <- pgls(rsq ~ median.trans, data = cd)
}

plot(dat$median.trans, 
     dat$rsq, 
     xlab = "expansion recency",
     ylab = "rsq")
title("placeholder")
slope <- summary(model)$coefficients[2, 1]
int <- summary(model)$coefficients[1, 1]
abline(int, slope)















# all clades; all repeats

library(phytools)
library(caper)
dat <- read.csv("../results/parsed.csv")
dat <- dat[dat$thrs == 0.9, ]
dat <- dat[!duplicated(dat$species), ]
dat$median.trans <- 1 - (dat[["total.rep.median"]]/70)
dat <- na.omit(dat[, c("species", "rsq", "clade", "median.trans")])
dat <- dat[dat$species != "Callithrix jacchus", ]
plot(dat$median.trans, dat$rsq)
model <- glm(rsq ~ median.trans, data = dat)
res <- setNames(resid(model), dat$species)

tree <- read.tree("../data/formatted_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
int <- intersect(tree$tip.label, dat$species)
pruned.tree <- keep.tip(tree, int)

if (phylosig(pruned.tree, res, method = "lambda", test = TRUE)[[4]] < 0.05) {
  cd <- comparative.data(pruned.tree, dat, names.col = "species")
  model <- pgls(rsq ~ median.trans, data = cd)
}

plot(dat$median.trans, 
     dat$rsq, 
     xlab = "expansion recency",
     ylab = "rsq")
title("placeholder")
slope <- summary(model)$coefficients[2, 1]
int <- summary(model)$coefficients[1, 1]
abline(int, slope)























# all clades; LINEs

library(phytools)
library(caper)
dat <- read.csv("../results/parsed.csv")
dat <- dat[dat$thrs == 0.9, ]
dat <- dat[!duplicated(dat$species), ]
dat$median.trans <- 1 - (dat[["ltr.rep.median"]]/70)
dat <- na.omit(dat[, c("species", "rsq", "clade", "median.trans")])
dat <- dat[dat$species != "Callithrix jacchus", ]
plot(dat$median.trans, dat$rsq)
model <- glm(rsq ~ median.trans, data = dat)
res <- setNames(resid(model), dat$species)

tree <- read.tree("../data/formatted_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
int <- intersect(tree$tip.label, dat$species)
pruned.tree <- keep.tip(tree, int)

if (phylosig(pruned.tree, res, method = "lambda", test = TRUE)[[4]] < 0.05) {
  cd <- comparative.data(pruned.tree, dat, names.col = "species")
  model <- pgls(rsq ~ median.trans, data = cd)
}

plot(dat$median.trans, 
     dat$rsq, 
     xlab = "expansion recency",
     ylab = "rsq")
title("placeholder")
slope <- summary(model)$coefficients[2, 1]
int <- summary(model)$coefficients[1, 1]
abline(int, slope)

