
# This model does not work with mammals but works for
# other clades as well as overall


# transform results
dat <- read.csv("../results/vertebrates/parsed.csv")
dat <- dat[!duplicated(dat$species), ]
dat <- dat[!is.na(dat$chromnum.1n), ]
dat$median.trans <- 1 - (dat$median/70)
dat$totalrep.prop <- dat$totalrep.pct * 0.01

# subset data
dat <- na.omit(dat[, c("species", "rsq", "clade", "median.trans", "totalrep.prop", "est.gnsz.Mbp", "asmblysize.Mbp")])
# dat <- dat[dat$clade == "Mammalia", ]

# remove assembly with bloated assembly size
dat <- dat[dat$species != "Callithrix jacchus", ]

# assign weights; weights approach 1 as assembly sizes approach genome sizes. weights tend toward 0 as assembly sizes deviate from genome sizes
dat$w <- 1 - (abs(dat$asmblysize.Mbp - dat$est.gnsz.Mbp) / dat$est.gnsz.Mbp)

# add an exponent to emphasize high quality genomes
dat$w <- dat$w^20

# visualize weights
hist(dat$w)

# model
model <- glm(rsq ~ totalrep.prop + totalrep.prop:median.trans, weights = dat$w, data = dat)

library(plotly)
x <- seq(min(dat$totalrep.prop), max(dat$totalrep.prop), length.out = 50)
y <- seq(min(dat$median.trans), max(dat$median.trans), length.out = 50)
grid <- expand.grid(totalrep.prop = x, median.trans = y)
grid$rsq <- predict(model, newdata = grid, type = "response")
z <- matrix(grid$rsq, nrow = 50, ncol = 50)
plot <- plot_ly(x = ~x, y = ~y, z = ~z) %>%
  add_surface() %>%
  layout(scene = list(
    xaxis = list(title = "Repeat Content"),
    yaxis = list(title = "Repeat Recency"),
    zaxis = list(title = "GD Variation")
  ))
plot


# phylogenetic signal test (no phylogenetic signals within clades and none overall)
library(phytools)
tree <- read.tree("../data/vertebrates/formatted_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
int <- intersect(dat$species, tree$tip.label)
pruned.tree <- keep.tip(tree, int)
pruned.dat <- dat[dat$species %in% int, ]
pruned.model <- glm(pruned.dat$rsq ~ pruned.dat$totalrep.prop + pruned.dat$totalrep.prop:pruned.dat$median.trans, weights = pruned.dat$w)
res <- setNames(resid(pruned.model), pruned.dat$species)
phylosig.p <- phylosig(pruned.tree, res, method = "lambda", test = TRUE)[4]$P


