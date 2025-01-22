


# rsq scatter
dat <- read.csv("../results/parsed.csv")
dat <- dat[dat$thrs == 0.8, ]
dat <- dat[dat$clade %in% c("Mammalia", "Sauria", "Actinopterygii"), ]
dat <- dat[!is.na(dat$chromnum.1n), ]
sub <- dat[dat$species %in% c("Bos taurus", "Peromyscus maniculatus bairdii"), ]
cols <- c("#af8dc3", "#7fbf7b")[as.factor(sub$species)]
par(mar = c(6, 4, 2, 2)+0.1)
plot(x = sub$size.Mbp,
y = sub$genecount,
xlab = NA,
ylab = "Gene Count",
col = cols,
pch = 16)
text(91, -150, "Chromosome Size (Mb)", xpd = NA, adj = 0)
text(60, -600, "Bos taurus", xpd = NA, adj = 0)
points(58, -610, pch = 16, col = "#af8dc3", xpd = NA, cex = 1.5)
text(130, -600, "Peromyscus maniculatus bairdii", xpd = NA, adj = 0)
points(128, -600, pch = 16, col = "#7fbf7b", xpd = NA, cex = 1.5)
x1 <- sub[sub$species == "Bos taurus", ]$size.Mbp
y1 <- sub[sub$species == "Bos taurus", ]$genecount
model1 <- lm(y1 ~ x1)
abline(model1,
col = "#af8dc3",
lwd = 1.5)
x2 <- sub[sub$species == "Peromyscus maniculatus bairdii", ]$size.Mbp
y2 <- sub[sub$species == "Peromyscus maniculatus bairdii", ]$genecount
model2 <- lm(y2 ~ x2)
abline(model2,
col = "#7fbf7b",
lwd = 1.5)
par(mar = c(5.1, 4.1, 4.1, 2.1))







# rsq barplot
rsq1 <- summary(lm(genecount ~ size.Mbp, data = sub[sub$species == "Bos taurus", ]))$r.squared
rsq2 <- summary(lm(genecount ~ size.Mbp, data = sub[sub$species == "Peromyscus maniculatus bairdii", ]))$r.squared
barplot(height = setNames(c(rsq1, rsq2), c("Domestic cattle", "North American deer mouse")),
col = c("#af8dc3", "#7fbf7b"),
ylim = c(0, 1),
ylab = "R-squared")










library(beeswarm)
dat <- read.csv("../results/parsed.csv")
dat <- dat[dat$thrs == 0.8, ]
dat <- dat[!duplicated(dat$species), ]
dat <- dat[dat$clade %in% c("Mammalia", "Sauria", "Actinopterygii"), ]
dat <- dat[!is.na(dat$chromnum.1n), ]
dat$clade <- factor(dat$clade, levels = c("Mammalia", "Actinopterygii", "Sauria"))
map <- c("Mammalia" = "#d95f02", "Actinopterygii" = "#7570b3", "Sauria" = "#1b9e77")
cols <- map[dat$clade]
nmam <- nrow(dat[dat$clade == "Mammalia", ])
nfish <- nrow(dat[dat$clade == "Actinopterygii", ])
nrep <- nrow(dat[dat$clade == "Sauria", ])
beeswarm(rsq ~ clade,
xlab = NA,
ylab = "R-squared",
data = dat,
pwcol = cols,
pch = 16,
spacing = 1.4,
labels = NA)
text(1, 0.05, "Mammals", xpd = NA, adj = c(0.5, 0.5))
text(2, 0.05, "Ray-finned fish", xpd = NA, adj = c(0.5, 0.5))
text(3, 0.05, "Reptiles", xpd = NA, adj = c(0.5, 0.5))
text(1, -0.025, paste0("n=", nmam), xpd = NA, adj = c(0.5, 0.5))
text(2, -0.025, paste0("n=", nfish), xpd = NA, adj = c(0.5, 0.5))
text(3, -0.025, paste0("n=", nrep), xpd = NA, adj = c(0.5, 0.5))
# phylogenetic anova
library(phytools)
dat <- na.omit(dat[, c("species", "clade", "rsq")])
tree <- read.tree("../data/formatted_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
int <- intersect(dat$species, tree$tip.label)
pruned.tree <- keep.tip(tree, int)
sig <- phylosig(pruned.tree,
setNames(dat$rsq, dat$species),
method = "lambda",
test = TRUE,
nsim = 10000)[[4]]
if (sig < 0.05) {
x <- setNames(dat$clade, dat$species)
y <- setNames(dat$rsq, dat$species)
phylANOVA(pruned.tree, x, y, nsim = 10000, posthoc = TRUE)
}


# estimated genome size
library(phytools)
library(caper)
dat <- read.csv("../results/parsed.csv")
dat <- dat[dat$thrs == 0, ]
dat <- dat[!duplicated(dat$species), ]
dat$clade <- factor(dat$clade, levels = c("Mammalia", "Actinopterygii", "Sauria"))
dat <- na.omit(dat[, c("species", "clade", "rsq", "est.gnsz.Mbp")])
map <- c("Mammalia" = "#d95f02", "Actinopterygii" = "#7570b3", "Sauria" = "#1b9e77")
cols <- map[dat$clade]
cols <- adjustcolor(cols, alpha.f = 0.5)
plot(rsq ~ est.gnsz.Mbp,
     data = dat,
     pch = 16,
     col = cols)
tree <- read.tree("../data/formatted_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
int <- intersect(dat$species, tree$tip.label)
pruned.tree <- keep.tip(tree, int)
model <- glm(rsq ~ est.gnsz.Mbp, data = dat)
res <- resid(model)
phylosig <- phylosig(pruned.tree,
                     setNames(res, dat$species),
                     method = "lambda",
                     test = TRUE,
                     nsim = 10000)[[4]]
if (phylosig < 0.05) {
  cd <- comparative.data(pruned.tree, dat, names.col = "species", vcv = TRUE)
  model <- pgls(rsq ~ est.gnsz.Mbp, data = cd)
}
abline(model,
       lwd = 1.5)
dat.mam <- dat[dat$clade == "Mammalia", ]
model.mam <- glm(rsq ~ est.gnsz.Mbp, data = dat.mam)
res.mam <- resid(model.mam)
phylosig.mam <- phylosig(pruned.tree,
                         setNames(res.mam, dat.mam$species),
                         method = "lambda",
                         test = TRUE,
                         nsim = 10000)[[4]]
if (phylosig.mam < 0.05) {
  cd.mam <- comparative.data(pruned.tree, dat.mam, names.col = "species", vcv = TRUE)
  model.mam <- pgls(rsq ~ est.gnsz.Mbp, data = cd.mam)
}
abline(model.mam,
       col = adjustcolor("#d95f02", alpha.f = 0.5),
       lwd = 1.5)
dat.fish <- dat[dat$clade == "Actinopterygii", ]
model.fish <- glm(rsq ~ est.gnsz.Mbp, data = dat.fish)
res.fish <- resid(model.fish)
phylosig.fish <- phylosig(pruned.tree,
                          setNames(res.fish, dat.fish$species),
                          method = "lambda",
                          test = TRUE,
                          nsim = 10000)[[4]]
if (phylosig.fish < 0.05) {
  cd.fish <- comparative.data(pruned.tree, dat.fish, names.col = "species", vcv = TRUE)
  model.fish <- pgls(rsq ~ est.gnsz.Mbp, data = cd.fish)
}
abline(model.fish,
       col = adjustcolor("#7570b3", alpha.f = 0.5),
       lwd = 1.5)
dat.rep <- dat[dat$clade == "Sauria", ]
model.rep <- glm(rsq ~ est.gnsz.Mbp, data = dat.rep)
res.rep <- resid(model.rep)
phylosig.rep <- phylosig(pruned.tree,
                         setNames(res.rep, dat.rep$species),
                         method = "lambda",
                         test = TRUE,
                         nsim = 10000)[[4]]
if (phylosig.rep < 0.05) {
  cd.rep <- comparative.data(pruned.tree, dat.rep, names.col = "species", vcv = TRUE)
  model.rep <- pgls(rsq ~ est.gnsz.Mbp, data = cd.rep)
}
abline(model.rep,
       col = adjustcolor("#1b9e77", alpha.f = 0.5),
       lwd = 1.5)

