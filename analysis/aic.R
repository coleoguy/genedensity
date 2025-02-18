

library(MuMIn)
library(phytools)
library(caper)

# identify repeats of interest
terms <- c("dna", "line", "ltr", "sine", "others", "unknown")

combs <- unlist(lapply(1:length(terms), function(x) {
  apply(combn(terms, x), 2, paste, collapse = ".")
}))
dat <- read.csv("../results/parsed.csv")
dat <- dat[!is.na(dat$chromnum.1n) & !duplicated(dat$species), ]
dat <- na.omit(dat[, c("species", "clade", "rsq", "rep.prop.total", paste0("rep.prop.", terms), "rep.age.total", paste0("rep.age.", combs))])
dat <- dat[dat$clade == "Mammalia", ]

# prune tree
tree <- read.tree("../data/formatted.tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
int <- intersect(dat$species, tree$tip.label)
dat <- dat[dat$species %in% int, ]
pruned.tree <- keep.tip(tree, int)

# fit PGLS models
pgls.models <- data.frame()
for (i in combs) {
  
  # get repeat age and proportion
  rep <- unlist(strsplit(i, "\\."))
  prop <- as.numeric(rowSums(as.data.frame(dat[, c(paste0("rep.prop.", rep))])))
  age <- dat[, c(paste0("rep.age.", i))]
  
  # create new dataframe to make cd object for PGLS; normalize data
  sub <- dat[, c("species", "clade", "rsq")]
  sub$age.norm <- (age - range(age)[1]) / diff(range(age))
  sub$prop.norm <- (prop - range(prop)[1]) / diff(range(prop))
  
  # fit model
  cd <- comparative.data(pruned.tree, sub, names.col = "species", vcv = TRUE)
  model <- dredge(pgls(rsq ~ age.norm*prop.norm, data = cd), subset = dc(x1, x2, x1:x2))
  pgls.models <- rbind(pgls.models, model)
}
# record model terms
pgls.models$model <- rep(combs, each = 5)

write.csv(as.data.frame(pgls.models), "../results/aic.csv", row.names = FALSE)
pgls.models <- read.csv("../results/aic.csv")














################ for GLS models without phylogenetic correction ################

library(MuMIn)
library(phytools)

terms <- c("dna", "line", "ltr", "sine", "others", "unknown")
combs <- unlist(lapply(1:length(terms), function(x) {
  apply(combn(terms, x), 2, paste, collapse = ".")
}))
dat <- read.csv("../results/parsed.csv")
dat <- dat[!is.na(dat$chromnum.1n) & !duplicated(dat$species), ]
dat <- na.omit(dat[, c("species", "clade", "rsq", "rep.prop.total", paste0("rep.prop.", terms), "rep.age.total", paste0("rep.age.", combs))])
dat <- dat[dat$clade == "Mammalia", ]

# fit GLS models
gls.models <- data.frame()
for (i in combs) {
  rep <- unlist(strsplit(i, "\\."))
  prop <- as.numeric(rowSums(as.data.frame(dat[, c(paste0("rep.prop.", rep))])))
  age <- dat[, c(paste0("rep.age.", i))]
  age.norm <- (age - range(age)[1]) / diff(range(age))
  prop.norm <- (prop - range(prop)[1]) / diff(range(prop))
  model <- dredge(glm(dat$rsq ~ age.norm*prop.norm, na.action = na.fail), subset = dc(x1, x2, x1:x2))
  gls.models <- rbind(gls.models, model)
}
gls.models$model <- rep(combs, each = 5)



