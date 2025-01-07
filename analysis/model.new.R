
# SINEs: 

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















#best model for SINEs according to step(); only the median is significant

# model
model <- glm(formula = rsq ~ median.trans + chromnum.1n + median.trans:chromnum.1n, data = dat)

# convert to plotting format
x <- seq(min(dat$median.trans), max(dat$median.trans), length.out = 100)
y <- seq(min(dat$chromnum.1n), max(dat$chromnum.1n), length.out = 100)
grid <- expand.grid(median.trans = x, chromnum.1n = y)
grid$rsq <- predict(model, newdata = grid, type = "response")
z <- matrix(grid$rsq, nrow = length(x), ncol = length(y))

# plot
par(mar = c(4, 4, 3, 8) + 0.1)
image(x = x, 
      y = y, 
      z = z, 
      col = viridis(100), 
      xlab = "",
      ylab = "",
      main = "Title")
mtext("expansion recency", side=1, line=2.5)
mtext("chromosome number", side=2, line=2.5)
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
mtext("Predicted consistency", side=4, line=2.5)
par(mar = c(5, 4, 4, 2))










