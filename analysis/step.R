
# transform results
dat <- read.csv("../results/vertebrates/parsed.csv")
dat <- dat[!duplicated(dat$species), ]
dat$median.trans <- 1 - (dat$median/70)
dat$totalrep.prop <- dat$totalrep.pct * 0.01

# subset data and perform stepwise model selection
d <- na.omit(dat[, c("species", "rsq", "clade", "median.trans", "totalrep.prop")])
summary(step(glm(d$rsq ~ d$median.trans * d$totalrep.prop))) # d$rsq ~ d$median.trans * d$totalrep.prop

# stepwise selection for each clade
m <- d[d$clade == "Mammalia", ]
summary(step(glm(m$rsq ~ m$median.trans * m$totalrep.prop))) # m$rsq ~ m$median.trans
f <- d[d$clade == "Actinopterygii", ]
summary(step(glm(f$rsq ~ f$median.trans * f$totalrep.prop))) # f$rsq ~ f$median.trans * f$totalrep.prop
r <- d[d$clade == "Sauria", ]
summary(step(glm(r$rsq ~ r$median.trans * r$totalrep.prop))) # r$rsq ~ 1

# prune tree
library(phytools)
tree <- read.tree("../data/vertebrates/formatted_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
int <- intersect(d$species, tree$tip.label)
pruned.tree <- keep.tip(tree, int)

# test for phylogenetic signal in overall model
res <- setNames(resid(step(glm(d$rsq ~ d$median.trans * d$totalrep.prop))), d$species)
phylosig(pruned.tree, res, method = "lambda", test = TRUE)

# phylogenetic signals in each clade
res <- setNames(resid(glm(m$rsq ~ m$median.trans)), m$species)
phylosig(pruned.tree, res, method = "lambda", test = TRUE) # mammals
res <- setNames(resid(glm(f$rsq ~ f$median.trans * f$totalrep.prop)), f$species)
phylosig(pruned.tree, res, method = "lambda", test = TRUE) # fish
res <- setNames(resid(glm(r$rsq ~ 1)), r$species)
phylosig(pruned.tree, res, method = "lambda", test = TRUE) # reptiles

# PGLS for mammals
library(caper)
cd <- comparative.data(tree, m, species)
summary(pgls(rsq ~ median.trans, data = cd))

# did step() overlook this model for mammals...
summary(pgls(rsq ~ totalrep.prop * median.trans, data = cd))

# ...and for fish?
cd <- comparative.data(tree, f, species)
summary(pgls(rsq ~ totalrep.prop * median.trans, data = cd))


