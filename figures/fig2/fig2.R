
# pics

library(beeswarm)
dat <- read.csv("../../results/rsq.csv")
clades <- c("Mammalia", "Actinopterygii", "Sauropsida")
dat <- dat[dat$clade %in% clades, ]
dat$clade <- factor(dat$clade, levels = clades)
map <- setNames(c("#d95f02", "#7570b3", "#1b9e77"), clades)
cols <- map[dat$clade]
num.sp <- sapply(clades, function(cl) {
  nrow(dat[dat$clade == cl, ])})
# plot
beeswarm(rsq ~ clade, 
         xlab = NA, 
         ylab = NA, 
         data = dat,
         pwcol = cols,
         pch = 16,
         spacing = 1.2,
         xaxt = "n", yaxt = "n")
# xaxt
axis(side   = 1,
     at     = seq(num.sp),
     labels = paste0("n=", num.sp), 
     mgp = c(3, 0.9, 0))
# yaxt
axis(side = 2, mgp = c(3, 0.9, 0))
mtext(expression(R^2), side = 2, line = 2.4)








# no pics

library(beeswarm)
dat <- read.csv("../../results/rsq.csv")
clades <- c("Mammalia", "Actinopterygii", "Sauropsida")
dat <- dat[dat$clade %in% clades, ]
dat$clade <- factor(dat$clade, levels = clades)
map <- setNames(c("#d95f02", "#7570b3", "#1b9e77"), clades)
cols <- map[dat$clade]
num.sp <- sapply(clades, function(cl) {
  nrow(dat[dat$clade == cl, ])})
# plot
beeswarm(rsq ~ clade, 
         xlab = NA, 
         ylab = NA, 
         data = dat,
         pwcol = cols,
         pch = 16,
         spacing = 1.4,
         xaxt = "n", yaxt = "n")
# xaxt
axis(side   = 1,
     at     = seq(num.sp),
     labels = c("Mammals", "Ray-finned fish", "Reptiles"), 
     mgp = c(3, 0.9, 0))
axis(side   = 1,
     at     = seq(num.sp),
     labels = paste0("n=", num.sp), 
     mgp = c(3, 2.1, 0))
# yaxt
axis(side = 2, mgp = c(3, 0.9, 0))
mtext(expression(R^2), side = 2, line = 2.4)


# phylogenetic anova
library(phytools)
tree <- read.tree("../../data/formatted-tree.nwk")
int <- intersect(dat$species, tree$tip.label)
pruned.tree <- keep.tip(tree, int)
sig <- phylosig(pruned.tree,
                setNames(dat$rsq, dat$species),
                method = "lambda",
                test = TRUE,
                niter = 10000)[[4]]
if (sig < 0.05) {
  x <- setNames(dat$clade, dat$species)
  y <- setNames(dat$rsq, dat$species)
  phylANOVA(pruned.tree, x, y, nsim = 10000, posthoc = TRUE)
}











