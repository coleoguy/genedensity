


























# mammal total
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
      main = "mammals; total; high contiguity (cutoff = 0.9)")
mtext("repeat content (prop. of assembly size)", side=1, line=2.5)
mtext("chromosome number (1n)", side=2, line=2.5)
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
mtext("Predicted consistency (R2)", side=4, line=2.5)
par(mar = c(5, 4, 4, 2))






























# mammal ltr: at high chromosome numbers, increase in expansion recency causes rsq to decrease and vice versa. at low expansion recency, rsq is consistent across chromosome numbers. at high expansion recency, rsq is 0.98 at n=39 and 0.12 at n=9
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
      main = "mammal; ltr; high contiguity (cutoff = 0.8)")
mtext("expansion recency (1 - [(divergence of repeat w/median length)/70])", side=1, line=2.5)
mtext("chromosome number (1n)", side=2, line=2.5)
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
mtext("Predicted consistency (R2)", side=4, line=2.5)
par(mar = c(5, 4, 4, 2))
















































# LINEs for all clades
library(nlme)
library(phytools)
library(plotly)
library(viridis)
dat <- read.csv("../results/parsed.csv")
dat <- dat[!duplicated(dat$species), ]
dat <- dat[!is.na(dat$chromnum.1n), ]
dat <- dat[dat$clade %in% c("Mammalia", "Sauria", "Actinopterygii"), ]
dat$median.trans <- 1 - (dat[["line.rep.median"]]/70)
dat$rep.prop <- dat$line.rep.pct / 100
dat <- na.omit(dat[, c("species", "rsq", "clade", "rep.prop", "median.trans", "chromnum.1n")])
dat <- dat[dat$species != "Callithrix jacchus", ]
tree <- read.tree("../data/formatted_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
int <- intersect(tree$tip.label, dat$species)
pruned.tree <- keep.tip(tree, int)
dat <- dat[dat$species %in% int, ]
dat <- dat[match(pruned.tree$tip.label, dat$species), ]
model <- gls(rsq ~ median.trans + chromnum.1n + rep.prop + median.trans:chromnum.1n + rep.prop:chromnum.1n, data = dat, method = "ML")
dat$predicted <- as.numeric(predict(model, newdata = dat, type = "response"))
dat$residuals <- dat$rsq - dat$predicted

# Plot predicted
plot_ly(
  data = dat,
  x = ~chromnum.1n,
  y = ~median.trans,
  z = ~rep.prop,
  color = ~predicted,
  colors = viridis(3),
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3)
) %>%
  layout(
    title = "predicted",
    scene = list(
      xaxis = list(title = "chromosome number (1n)"),
      yaxis = list(title = "expansion recency (1 - [(divergence of repeat w/median length)/70])"),
      zaxis = list(title = "repeat content (prop. of assembly size)")
    ),
    colorbar = list(title = "Gene Density Variation")
  )

a <- setNames(factor(dat$clade, levels = unique(dat$clade)), dat$species)
b <- setNames(dat$rsq, dat$species)
anova <- aov.phylo(a ~ b, pruned.tree)
anova <- phylANOVA(pruned.tree, b, a)

# plot prediction grid
x <- seq(min(dat$chromnum.1n), max(dat$chromnum.1n), length.out = 50)
y <- seq(min(dat$median.trans), max(dat$median.trans), length.out = 50)
z <- seq(min(dat$rep.prop), max(dat$rep.prop), length.out = 50)
grid <- expand.grid(chromnum.1n = x, median.trans = y, rep.prop = z)
grid$rsq <- predict(model, newdata = grid, type = "response")
plot_ly(
  data = grid,
  x = ~chromnum.1n,
  y = ~median.trans,
  z = ~rep.prop,
  color = ~rsq,
  colors = viridis(100),
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3)
) %>%
  layout(
    title = "predicted",
    scene = list(
      xaxis = list(title = "chromosome number (1n)"),
      yaxis = list(title = "expansion recency (1 - [(divergence of repeat w/median length)/70])"),
      zaxis = list(title = "repeat content (prop. of assembly size)")
    ),
    colorbar = list(title = "gene density consistency (R2)")
  )

# Plot observed
plot_ly(
  data = vis,
  x = ~chromnum.1n,
  y = ~median.trans,
  z = ~rep.prop,
  color = ~rsq,
  visible = ~visible,
  colors = viridis(100),
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3)
) %>%
  layout(
    title = "observed",
    scene = list(
      xaxis = list(title = "chromosome number (1n)"),
      yaxis = list(title = "expansion recency (1 - [(divergence of repeat w/median length)/70])"),
      zaxis = list(title = "repeat content (prop. of assembly size)")
    ),
    colorbar = list(title = "Gene Density Variation")
  )

# Plot residuals
plot_ly(
  data = dat,
  x = ~chromnum.1n,
  y = ~median.trans,
  z = ~rep.prop,
  color = ~residuals,
  colors = viridis(3),
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3)
) %>%
  layout(
    title = "residuals",
    scene = list(
      xaxis = list(title = "chromosome number (1n)"),
      yaxis = list(title = "expansion recency (1 - [(divergence of repeat w/median length)/70])"),
      zaxis = list(title = "repeat content (prop. of assembly size)")
    ),
    colorbar = list(title = "Gene Density Variation")
  )

# Plot clades
plot_ly(
  data = dat,
  x = ~chromnum.1n,
  y = ~median.trans,
  z = ~rep.prop,
  color = ~clade,
  colors = viridis(3),
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3)
) %>%
  layout(
    title = "clades",
    scene = list(
      xaxis = list(title = "chromosome number (1n)"),
      yaxis = list(title = "expansion recency (1 - [(divergence of repeat w/median length)/70])"),
      zaxis = list(title = "repeat content (prop. of assembly size)")
    ),
    colorbar = list(title = "Gene Density Variation")
  )






























terms <- c(
  "intercept", 
  "chromnum.1n", 
  "median.trans", 
  "rep.prop", 
  "chromnum.1n:median.trans", 
  "chromnum.1n:rep.prop", 
  "median.trans:rep.prop", 
  "chromnum.1n:median.trans:rep.prop"
)
models <- read.csv("../results/models.csv")
names(models) <- c("clade", "thrs", "repeat", "phylosig", "stat", terms)
models$mean <- apply(models[, terms[-1]], 1, function(x) mean(x, na.rm = T))
models <- models[models$stat == "beta", ]
models <- models[models$clade == "Mammalia", ]
models <- models[models$`repeat` == "total", ]

