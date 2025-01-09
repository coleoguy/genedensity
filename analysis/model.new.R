

classes <- c("total", "line", "sine", "ltr", "dna", "rc")
for (qwe in classes) {
  dat <- read.csv("../results/parsed.csv")
  dat <- dat[!duplicated(dat$species), ]
  dat <- dat[!is.na(dat$chromnum.1n), ]
  dat$median.trans <- 1 - (dat[[paste0(qwe, ".rep.median")]]/70)
  dat$rep.prop <- dat[[paste0(qwe, ".rep.pct")]] / 100
  dat <- na.omit(dat[, c("species", "rsq", "clade", "median.trans", "rep.prop", "chromnum.1n")])
  dat <- dat[dat$clade == "Mammalia", ]
  dat <- dat[dat$species != "Callithrix jacchus", ]
  library(phytools)
  tree <- read.tree("../data/formatted_tree.nwk")
  tree$tip.label <- gsub("_", " ", tree$tip.label)
  int <- intersect(tree$tip.label, dat$species)
  pruned.tree <- keep.tip(tree, int)
  dat1 <- dat[dat$species %in% int, ]
  dat1 <- dat1[match(pruned.tree$tip.label, dat1$species), ]
  formula <- rsq ~ rep.prop * median.trans * chromnum.1n
  sink("NUL")
  model <- step(glm(formula, data = dat1))
  sink()
  res <- setNames(resid(model), dat1$species)
  signal <- phylosig(pruned.tree, res, method="lambda", test=TRUE)[4]
  if (signal < 0.05) {
    library(nlme)
    library(MASS)
    sink("NUL")
    model <- stepAIC(gls(formula, data = dat1, method = "ML"), direction = "both")
    sink()
    print(paste0("============== ", qwe, " =============="))
    #print(summary(model)$call)
    print(summary(model)$tTable) 
  } else {
    sink("NUL")
    model <- step(glm(formula, data = dat)) 
    sink()
    print(paste0("============== ", qwe, " =============="))
    #print(summary(model)$call)  
    print(summary(model)$coefficients) 
  }
}
































classes <- c("total", "line", "sine", "ltr", "dna", "rc")
for (qwe in classes) {
  dat <- read.csv("../results/parsed.csv")
  dat <- dat[!duplicated(dat$species), ]
  dat <- dat[!is.na(dat$chromnum.1n), ]
  dat$median.trans <- 1 - (dat[[paste0(qwe, ".rep.median")]]/70)
  dat$rep.prop <- dat[[paste0(qwe, ".rep.pct")]] / 100
  dat <- na.omit(dat[, c("species", "rsq", "clade", "median.trans", "rep.prop", "chromnum.1n")])
  dat <- dat[dat$clade == "Mammalia", ]
  dat <- dat[dat$species != "Callithrix jacchus", ]
  library(phytools)
  tree <- read.tree("../data/formatted_tree.nwk")
  tree$tip.label <- gsub("_", " ", tree$tip.label)
  int <- intersect(tree$tip.label, dat$species)
  pruned.tree <- keep.tip(tree, int)
  dat1 <- dat[dat$species %in% int, ]
  dat1 <- dat1[match(pruned.tree$tip.label, dat1$species), ]
  formula <- rsq ~ rep.prop * median.trans * chromnum.1n
  sink("NUL")
  model <- step(glm(formula, data = dat1))
  sink()
  res <- setNames(resid(model), dat1$species)
  signal <- phylosig(pruned.tree, res, method="lambda", test=TRUE)[4]
  if (signal < 0.05) {
    library(nlme)
    library(MASS)
    sink("NUL")
    model <- stepAIC(gls(formula, data = dat1, method = "ML"), direction = "both")
    sink()
    print(paste0("============== ", qwe, " =============="))
    #print(summary(model)$call)
    print(summary(model)$tTable) 
    
    
    terms <- c(
      "chromnum.1n", 
      "median.trans", 
      "rep.prop", 
      "chromnum.1n:median.trans", 
      "chromnum.1n:rep.prop", 
      "median.trans:rep.prop", 
      "chromnum.1n:median.trans:rep.prop"
    )
    
    effects <- data.frame(summary(model)$tTable)[terms, ]
    # clade, threshold, repeat, p or beta, results
    
    p <- c("p", effects$p.value)
    beta <- c("beta", effects$Value)
    
    
  } else {
    sink("NUL")
    model <- step(glm(formula, data = dat)) 
    sink()
    
    effects <- data.frame(summary(model)$coefficients)[terms, ]
    # clade, threshold, repeat, p or beta, results
    
    p <- c("p", effects$Pr...t..)
    beta <- c("beta", effects$Estimate)
  }
}


































































# total
# model
library(phytools)
dat <- read.csv("../results/parsed.csv")
dat <- dat[!duplicated(dat$species), ]
dat <- dat[!is.na(dat$chromnum.1n), ]
dat$rep.prop <- dat$total.rep.pct / 100
dat <- na.omit(dat[, c("species", "rsq", "clade", "rep.prop", "chromnum.1n")])
dat <- dat[dat$clade == "Mammalia", ]
dat <- dat[dat$species != "Callithrix jacchus", ]
tree <- read.tree("../data/formatted_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
int <- intersect(tree$tip.label, dat$species)
pruned.tree <- keep.tip(tree, int)
dat <- dat[dat$species %in% int, ]
dat <- dat[match(pruned.tree$tip.label, dat$species), ]
formula <- rsq ~ rep.prop * chromnum.1n
library(nlme)
model <- gls(formula, data = dat, method = "ML")

# predict
x <- seq(min(dat$rep.prop), max(dat$rep.prop), length.out = 100)
y <- seq(min(dat$chromnum.1n), max(dat$chromnum.1n), length.out = 100)
grid <- expand.grid(rep.prop= x, chromnum.1n = y)
grid$rsq <- predict(model, newdata = grid, type = "response")
z <- matrix(grid$rsq, nrow = length(x), ncol = length(y))
z.range <- seq(min(z), max(z), length.out = 100)

# plot
library(viridis)
z.labels <- seq(min(z), max(z), length.out = 7)
par(mar = c(4, 4, 3, 8))
image(x = x, 
      y = y, 
      z = z, 
      col = viridis(100), 
      xlab = "",
      ylab = "",
      main = "Total")
mtext("repeat content", side=1, line=2.5)
mtext("chromosome number", side=2, line=2.5)
contour(x, 
        y, 
        z, 
        levels = z.labels, 
        add = TRUE,
        col = "black", 
        lwd = 1,
        drawlabels = FALSE)
par(new = TRUE)
par(mar = c(4, 25, 3, 6))
image(1, 
      z.range, 
      t(matrix(z.range)), 
      col = viridis(100), 
      xaxt = "n",
      yaxt = "n", 
      xlab = "", 
      ylab = "")
#axis(4, 
#     at = c(round(min(z), 2), 0.4, 0.6, 0.8, round(max(z), 2)), 
#     labels = c(round(min(z), 2), 0.4, 0.6, 0.8, round(max(z), 2)))
axis(4, at = z.labels, labels = round(z.labels, 2), cex.axis = 0.9)
mtext("Predicted consistency", side=4, line=2.5)
par(mar = c(5, 4, 4, 2))
































# line: nothing significant
# best model for SINEs according to step(); only the median is significant

# model
dat <- read.csv("../results/parsed.csv")
dat <- dat[!duplicated(dat$species), ]
dat <- dat[!is.na(dat$chromnum.1n), ]
dat$median.trans <- 1 - (dat[["sine.rep.median"]]/70)
dat <- na.omit(dat[, c("species", "rsq", "clade", "median.trans", "chromnum.1n")])
dat <- dat[dat$clade == "Mammalia", ]
dat <- dat[dat$species != "Callithrix jacchus", ]
model <- glm(formula = rsq ~ median.trans + chromnum.1n + median.trans:chromnum.1n, data = dat)

# predict
x <- seq(min(dat$median.trans), max(dat$median.trans), length.out = 100)
y <- seq(min(dat$chromnum.1n), max(dat$chromnum.1n), length.out = 100)
grid <- expand.grid(median.trans = x, chromnum.1n = y)
grid$rsq <- predict(model, newdata = grid, type = "response")
z <- matrix(grid$rsq, nrow = length(x), ncol = length(y))
z.range <- seq(min(z), max(z), length.out = 100)

# plot
library(viridis)
z.labels <- seq(min(z), max(z), length.out = 7)
par(mar = c(4, 4, 3, 8))
image(x = x, 
      y = y, 
      z = z, 
      col = viridis(100), 
      xlab = "",
      ylab = "",
      main = "SINEs")
mtext("expansion recency", side=1, line=2.5)
mtext("chromosome number", side=2, line=2.5)
contour(x, 
        y, 
        z, 
        levels = z.labels, 
        add = TRUE,
        col = "black", 
        lwd = 1,
        drawlabels = FALSE)
par(new = TRUE)
par(mar = c(4, 25, 3, 6))
image(1, 
      z.range, 
      t(matrix(z.range)), 
      col = viridis(100), 
      xaxt = "n",
      yaxt = "n", 
      xlab = "", 
      ylab = "")
#axis(4, 
#     at = c(round(min(z), 2), 0.4, 0.6, 0.8, round(max(z), 2)), 
#     labels = c(round(min(z), 2), 0.4, 0.6, 0.8, round(max(z), 2)))
axis(4, at = z.labels, labels = round(z.labels, 2), cex.axis = 0.9)
mtext("Predicted consistency", side=4, line=2.5)
par(mar = c(5, 4, 4, 2))









# ltr: at high chromosome numbers, increase in expansion recency causes rsq to decrease and vice versa. at low expansion recency, rsq is consistent across chromosome numbers. at high expansion recency, rsq is 0.98 at n=39 and 0.12 at n=9
dat <- read.csv("../results/parsed.csv")
dat <- dat[!duplicated(dat$species), ]
dat <- dat[!is.na(dat$chromnum.1n), ]
dat$median.trans <- 1 - (dat[["ltr.rep.median"]]/70)
dat <- na.omit(dat[, c("species", "rsq", "clade", "median.trans", "chromnum.1n")])
dat <- dat[dat$clade == "Mammalia", ]
dat <- dat[dat$species != "Callithrix jacchus", ]
model <- glm(rsq ~ median.trans + chromnum.1n + median.trans:chromnum.1n, data = dat)

# predict
library(viridis)
x <- seq(min(dat$median.trans), max(dat$median.trans), length.out = 100)
y <- seq(min(dat$chromnum.1n), max(dat$chromnum.1n), length.out = 100)
grid <- expand.grid(median.trans = x, chromnum.1n = y)
grid$rsq <- predict(model, newdata = grid, type = "response")
z <- matrix(grid$rsq, nrow = length(x), ncol = length(y))
z.range <- seq(min(z), max(z), length.out = 100)

# plot
z.labels <- seq(min(z), max(z), length.out = 7)
par(mar = c(4, 4, 3, 8))
image(x = x, 
      y = y, 
      z = z, 
      col = viridis(100), 
      xlab = "",
      ylab = "",
      main = "LTRs")
mtext("expansion recency", side=1, line=2.5)
mtext("chromosome number", side=2, line=2.5)
contour(x, 
        y, 
        z, 
        levels = z.labels, 
        add = TRUE,
        col = "black", 
        lwd = 1,
        drawlabels = FALSE)
par(new = TRUE)
par(mar = c(4, 25, 3, 6))
image(1, 
      z.range, 
      t(matrix(z.range)), 
      col = viridis(100), 
      xaxt = "n",
      yaxt = "n", 
      xlab = "", 
      ylab = "")
#axis(4, 
#     at = c(round(min(z), 2), 0.4, 0.6, 0.8, round(max(z), 2)), 
#     labels = c(round(min(z), 2), 0.4, 0.6, 0.8, round(max(z), 2)))
axis(4, at = z.labels, labels = round(z.labels, 2), cex.axis = 0.9)
mtext("Predicted consistency", side=4, line=2.5)
par(mar = c(5, 4, 4, 2))

