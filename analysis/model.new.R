
# before pgls: new parsing method kept 21 mammals and rsq~transformed median has slope -1.6988 and p value 0.016207; no need to exponentiate weights
# after pgls: 20 species, beta = 1.747489, p = 0.0183, r2 = 0.2731145, predicted rsq diff between highest and lowest medians: 0.37446190
# what about other clades?


# new model
dat <- read.csv("../results/parsed.csv")
dat <- dat[!duplicated(dat$species), ]
dat <- dat[!is.na(dat$chromnum.1n), ]
dat$median.trans <- 1 - (dat$total.rep.median/70)
dat <- na.omit(dat[, c("species", "rsq", "clade", "median.trans")])
dat <- dat[dat$clade == "Mammalia", ]
dat <- dat[dat$species != "Callithrix jacchus", ]
library(phytools)
tree <- read.tree("../data/formatted_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
int <- intersect(tree$tip.label, dat$species)
pruned.tree <- keep.tip(tree, int)
dat1 <- dat[dat$species %in% int, ]
dat1 <- dat1[match(pruned.tree$tip.label, dat1$species), ]
model <- glm(rsq ~ median.trans, data = dat1)
res <- setNames(resid(model), dat1$species)
phylosig(pruned.tree, res, method="lambda", test=TRUE)
library(nlme)
summary(gls(rsq ~ median.trans, 
            data = dat1))
library(piecewiseSEM)
rsquared(gls(rsq ~ median.trans, 
             data = dat1))


qwert <- list()
classes <- c("total", "line", "sine", "ltr", "dna", "rc")
for (qwe in classes) {
  dat <- read.csv("../results/parsed.csv")
  dat <- dat[!duplicated(dat$species), ]
  dat <- dat[!is.na(dat$chromnum.1n), ]
  dat$median.trans <- 1 - (dat[[paste0(qwe, ".rep.median")]]/70)
  dat <- na.omit(dat[, c("species", "rsq", "clade", "median.trans")])
  dat <- dat[dat$clade == "Mammalia", ]
  dat <- dat[dat$species != "Callithrix jacchus", ]
  library(phytools)
  tree <- read.tree("../data/formatted_tree.nwk")
  tree$tip.label <- gsub("_", " ", tree$tip.label)
  int <- intersect(tree$tip.label, dat$species)
  pruned.tree <- keep.tip(tree, int)
  dat1 <- dat[dat$species %in% int, ]
  dat1 <- dat1[match(pruned.tree$tip.label, dat1$species), ]
  model <- glm(rsq ~ median.trans, data = dat1)
  res <- setNames(resid(model), dat1$species)
  signal <- phylosig(pruned.tree, res, method="lambda", test=TRUE)[4]
  if (signal < 0.05) {
    library(nlme)
    qwert[[paste0(qwe, ".slope")]] <- summary(gls(rsq ~ median.trans, data = dat1))$tTable[2, 1]
    qwert[[paste0(qwe, ".p")]] <- summary(gls(rsq ~ median.trans, data = dat1))$tTable[2, 4]
    library(piecewiseSEM)
    qwert[[paste0(qwe, ".r2")]] <- rsquared(gls(rsq ~ median.trans, data = dat1))[[5]]
  } else {
    qwert[[paste0(qwe, ".slope")]] <- summary(glm(rsq ~ median.trans, data = dat))$coefficients[2, 1]
    qwert[[paste0(qwe, ".p")]] <- summary(glm(rsq ~ median.trans, data = dat))$coefficients[2, 4]
    library(piecewiseSEM)
    qwert[[paste0(qwe, ".r2")]] <- rsquared(glm(rsq ~ median.trans, data = dat))[[5]]
  }
}
qwert <- as.data.frame(qwert)








