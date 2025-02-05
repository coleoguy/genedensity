


library(phytools)
library(caper)
library(viridis)
dat <- read.csv("../results/parsed.csv")
dat <- dat[!duplicated(dat$species), ]
dat <- dat[!is.na(dat$chromnum.1n), ]
dat <- dat[dat$clade == "Mammalia", ]
dat <- na.omit(dat[, c("species", "rsq", "total.rep.median", "total.rep.pct")])
dat$total.rep.median <- 1 - (dat$total.rep.median/70)
dat$total.rep.pct <- dat$total.rep.pct / 100
tree <- read.tree("../data/formatted_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
int <- intersect(tree$tip.label, dat$species)
pruned.tree <- keep.tip(tree, int)
model <- glm(rsq ~ total.rep.pct + total.rep.median, data = dat)
sig <- phylosig(pruned.tree,
                setNames(resid(model), dat$species),
                method = "lambda",
                test = TRUE,
                nsim = 10000)[[4]]
if (sig < 0.05) {
  cd <- comparative.data(pruned.tree, dat, names.col = "species", vcv = TRUE)
  model <- pgls(rsq ~ total.rep.pct + total.rep.median, data = cd)
}
x <- seq(min(cd$data$total.rep.pct), max(cd$data$total.rep.pct), length.out = 100)
y <- seq(min(cd$data$total.rep.median), max(cd$data$total.rep.median), length.out = 100)
grid <- expand.grid(
  total.rep.pct = x,
  total.rep.median = y
)
grid$rsq <- predict(model, newdata = grid)
z <- matrix(grid$rsq, nrow = length(x), ncol = length(y))
image(x = x, 
      y = y, 
      z = z, 
      col = viridis(100), 
      xlab = "repeat proportion",
      ylab = "repeat recency",
      main = "mammals all repeats")

res <- 1000
norm <- (dat$rsq - min(dat$rsq)) / (max(dat$rsq) - min(dat$rsq))
palette <- viridis(res)
idx <- round(norm * res-1) + 1
cols <- palette[idx]

plot(dat$total.rep.pct, 
     dat$total.rep.median, 
     col = cols, 
     pch = 16)
plot(dat$total.rep.median, dat$rsq)
plot(dat$total.rep.pct, dat$rsq)


