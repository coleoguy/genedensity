

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
dat <- dat[dat$clade == "Actinopterygii", ]

# remove assembly with bloated assembly size
dat <- dat[dat$species != "Callithrix jacchus", ]

# assign weights; weights approach 1 as assembly sizes approach genome sizes. weights tend toward 0 as assembly sizes deviate from genome sizes
dat$w <- 1 - (abs(dat$asmblysize.Mbp - dat$est.gnsz.Mbp) / dat$est.gnsz.Mbp)

# add an exponent to emphasize high quality genomes
dat$w <- dat$w^20

# model
model <- glm(rsq ~ totalrep.prop + totalrep.prop:median.trans, weights = dat$w, data = dat)

# convert to plotting format
x <- seq(min(dat$totalrep.prop), max(dat$totalrep.prop), length.out = 100)
y <- seq(min(dat$median.trans), max(dat$median.trans), length.out = 100)
grid <- expand.grid(totalrep.prop = x, median.trans = y)
grid$rsq <- predict(model, newdata = grid, type = "response")
z <- matrix(grid$rsq, nrow = length(x), ncol = length(y))

# plot
original <- par(no.readonly = TRUE)
par(mar = c(4, 4, 3, 8) + 0.1)
image(x = x, 
      y = y, 
      z = z, 
      col = viridis(100), 
      xlab = "",
      ylab = "",
      main = "Title")
mtext("Repeat content", side=1, line=2.5)
mtext("Expansion recency", side=2, line=2.5)
contour(x, 
        y, 
        z, 
        add = TRUE,
        col = "black", 
        lwd = 1,
        drawlabels = FALSE)
par(new = TRUE)
par(mar = c(4, 25, 3, 6))
z_range <- seq(min(z), max(z), length.out = 100)
image(1, 
      z_range, 
      t(matrix(z_range)), 
      col = viridis(100), 
      xaxt = "n",
      yaxt = "n", 
      xlab = "", 
      ylab = "")
#axis(4, 
#     at = c(round(min(z), 2), 0.4, 0.6, 0.8, round(max(z), 2)), 
#     labels = c(round(min(z), 2), 0.4, 0.6, 0.8, round(max(z), 2)))
axis(4, at = pretty(z_range), labels = round(pretty(z_range), 2))
mtext("Predicted variation", side=4, line=2.5)
par(original)

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

# color by points
library(viridis)
cols <- viridis(length(unique(dat$w)), alpha = 0.45)[as.factor(dat$w)]
plot(dat$totalrep.prop, 
     dat$rsq,
     col = cols,
     pch = 16)
abline(glm(dat$rsq ~ dat$totalrep.prop, weights = dat$w), col = "blue")
abline(glm(dat$rsq ~ dat$totalrep.prop))
legend("bottomright", 
       legend = round(seq(min(dat$totalrep.prop), max(dat$totalrep.prop), length.out = 5), 2), 
       fill = viridis(5), 
       title = "Legend")

plot(dat$median.trans, 
     dat$rsq,
     col = cols,
     pch = 16)
abline(glm(dat$rsq ~ dat$median.trans, weights = dat$w), col = "blue")
abline(glm(dat$rsq ~ dat$median.trans))
legend("bottomright", 
       legend = round(seq(min(dat$median.trans), max(dat$median.trans), length.out = 5), 2), 
       fill = viridis(5), 
       title = "Legend")

plot(dat$median.trans * dat$totalrep.prop, 
     dat$rsq,
     col = cols,
     pch = 16)
term <- dat$median.trans * dat$totalrep.prop
abline(glm(dat$rsq ~ term, weights = dat$w), col = "blue")
abline(glm(dat$rsq ~ term))
legend("bottomright", 
       legend = round(seq(min(dat$median.trans * dat$totalrep.prop), max(dat$median.trans * dat$totalrep.prop), length.out = 5), 2), 
       fill = viridis(5), 
       title = "Legend")



