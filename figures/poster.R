


# rsq scatter
dat <- read.csv("../results/parsed.csv")
dat <- dat[dat$thrs == 0.8, ]
dat <- dat[dat$clade %in% c("Mammalia", "Sauria", "Actinopterygii"), ]
dat <- dat[!is.na(dat$chromnum.1n), ]
sub <- dat[dat$species %in% c("Bos taurus", "Peromyscus maniculatus bairdii"), ]
sub$genecount <- sub$genecount / 1000
cols <- c("#f4a582", "#0571b0")[as.factor(sub$species)]
par(mar = c(6, 4, 2, 2)+0.1)
plot(x = sub$size.Mbp,
     y = sub$genecount,
     xlab = NA,
     ylab = "Gene count (thousands)",
     ylim = c(0.3, 3), 
     xlim = c(30, 200), 
     col = cols,
     cex = 1.5, 
     pch = 16)
text(91, -150, "Chromosome size (Mb)", xpd = NA, adj = 0)
text(60, -600, "Bos taurus", xpd = NA, adj = 0)
points(58, -610, pch = 16, col = "#af8dc3", xpd = NA, cex = 1.5)
text(130, -600, "Peromyscus maniculatus bairdii", xpd = NA, adj = 0)
points(128, -600, pch = 16, col = "#7fbf7b", xpd = NA, cex = 1.5)
x1 <- sub[sub$species == "Bos taurus", ]$size.Mbp
y1 <- sub[sub$species == "Bos taurus", ]$genecount
model1 <- lm(y1 ~ x1)
abline(model1,
col = "#f4a582",
lwd = 1.5)
x2 <- sub[sub$species == "Peromyscus maniculatus bairdii", ]$size.Mbp
y2 <- sub[sub$species == "Peromyscus maniculatus bairdii", ]$genecount
model2 <- lm(y2 ~ x2)
abline(model2,
col = "#0571b0",
lwd = 1.5)
par(mar = c(5.1, 4.1, 4.1, 2.1))







# rsq barplot
rsq1 <- summary(lm(genecount ~ size.Mbp, data = sub[sub$species == "Bos taurus", ]))$r.squared
rsq2 <- summary(lm(genecount ~ size.Mbp, data = sub[sub$species == "Peromyscus maniculatus bairdii", ]))$r.squared
barplot(height = setNames(c(rsq1, rsq2), c("Domestic cattle", "North American deer mouse")),
col = c("#af8dc3", "#7fbf7b"),
ylim = c(0, 1),
ylab = "R-squared")









# r2 beeswarm
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














# largest repeats in each clade
library(data.table)
library(beeswarm)
library(viridis)
name <- cl <- r2 <- largest <- c()
dat <- read.csv("../data/data.csv")
dir <- "../results/individual_species_results"
files <- paste0(dir, "/",  list.files(dir))
contigs <- lapply(files, fread)
contigs <- as.data.frame(rbindlist((contigs), fill = TRUE))
contigs <- contigs[contigs$size.Mbp >= 10, ]
rm <- names(table(contigs$species)[table(contigs$species) < 3])
contigs <- contigs[!(contigs$species %in% rm), ]
df <- data.frame()
parsed <- data.frame()
for (z in unique(contigs$species)) {
  sub <- contigs[contigs$species == z, ]
  cont <- sum(sub$size.Mbp)
  total <- contigs[contigs$species == z, ]$asmblysize[1]
  if (cont/total >= 0) {
    parsed <- rbind(parsed, sub)
  }
}
sp <- unique(parsed$species)
for (species in sp) {
  sub <- parsed[which(parsed$species == species), ]
  # fit <- summary(glm(sub$genecount ~ sub$size.Mbp))
  # beta <- fit$coefficients[2, 1]
  # pval.beta <- fit$coefficients[2, 4]
  rsq <- summary(lm(sub$genecount ~ sub$size.Mbp))$r.squared
  # weightmean <- sum(sub$genedens * sub$size.Mbp) / sum(sub$size.Mbp)
  # weightsd <- sqrt(sum(sub$size.Mbp * (sub$genedens - weightmean)^2) / sum(sub$size.Mbp))
  # weightcv <- weightsd / weightmean
  stats <- data.frame(species, rsq, 0)
  sub <- merge(sub, stats, by = "species", all = TRUE)
  sub <- merge(dat[dat$species == species, ], sub, by = "species", all = TRUE)
  df <- rbind(df, sub)
}
df$clade <- df$class
df[df$clade %in% "Aves", ]$clade <- "Sauria"
df[df$clade %in% "Reptilia", ]$clade <- "Sauria"
df[!(df$clade %in% c("Actinopterygii", "Mammalia", "Sauria")), ]$clade <- "Others"

files <- list.files("../results/divsums")
sp <- gsub("_", " ", gsub(".divsum$", "", files))
asmbsz <- df[!duplicated(df$species), ]
asmbsz <- asmbsz[asmbsz$species %in% sp, ]
asmbsz <- setNames(asmbsz$asmblysize.Mbp*1000000, asmbsz$species)
repstats <- data.frame()
for (i in 1:length(sp)) {
  species <- sp[i]
  divsum <- readLines(paste0("../results/divsums/", files[i]))
  phrase <- "Coverage for each repeat class and divergence (Kimura)"
  start.index <- match(phrase, divsum) + 1
  divsum <- divsum[start.index:length(divsum)]
  divsum <- read.table(textConnection(divsum), 
                       sep = " ", 
                       header = TRUE)
  divsum <- divsum[-c(which(sapply(divsum, function(col) all(is.na(col)))))]
  divsum <- colSums(divsum)
  name <- c(name, species)
  cl <- c(cl, unique(df[df$species == species, ]$clade))
  r2 <- c(r2, unique(df[df$species == species, ]$rsq))
  largest <- c(largest, names(tail(sort(divsum), 1)))
}
df <- data.frame(name, cl, r2, largest)
df <- df[df$cl %in% c("Mammalia", "Actinopterygii", "Sauria"), ]
df$cl <- factor(df$cl, levels = c("Mammalia", "Actinopterygii", "Sauria"))
df$largest[df$largest %in% c("DNA.TcMar.Tc1", "LINE.L2", "SINE.tRNA.Core")] <- "Others"
map <- c("Mammalia" = "#d95f02", "Actinopterygii" = "#7570b3", "Sauria" = "#1b9e77")
cols <- map[df$cl]
par(mar = c(5, 4, 4, 7)+0.1)
beeswarm(r2 ~ largest, 
         data = df, 
         xlab = NA, 
         ylab = "R2", 
         pch = 16, 
         pwcol = cols, 
         spacing = 1.1, 
         at = c(1, 1.92, 3, 4.08, 5), 
         labels = c("LINE/CR1", "LINE/L1", "LINE/RTE-BovB", "Others", "Unclassified"))
text(1.5, -0.32, "Mammals", xpd = NA, adj = c(0.5, 0.5))
text(3.166667, -0.32, "Ray-finned fish", xpd = NA, adj = c(0.5, 0.5))
text(4.833334, -0.32, "Reptiles", xpd = NA, adj = c(0.5, 0.5))
points(1, -0.326, cex = 1.5, pch = 16, col = "#d95f02", xpd = NA, adj = 0)
points(2.46, -0.326, cex = 1.5, pch = 16, col = "#7570b3", xpd = NA)
points(4.39, -0.326, cex = 1.5, pch = 16, col = "#1b9e77", xpd = NA)
par(mar = c(5.1, 4.1, 4.1, 2.1))




# repeat landscape
dat <- read.csv("../results/parsed.csv")
dat <- dat[dat$thrs == 0, ]
dat <- dat[!duplicated(dat$species), ]
dat <- na.omit(dat[, c("species", "rsq", "clade", "asmblysize.Mbp")])
files <- list.files("../results/divsums")
sp <- gsub("_", " ", gsub(".divsum$", "", files))
asmbsz <- dat[dat$species %in% sp, ]
asmbsz <- setNames(asmbsz$asmblysize.Mbp*1000000, asmbsz$species)
cols <- c("#e31a1c", "#fb9a99", "#33a02c", "#b2df8a", "#1f78b4", "#a6cee3")
par(mar = c(3, 4, 1, 7)+0.1)
for (i in c(19, 138)) {
  species <- sp[i]
  divsum <- readLines(paste0("../results/divsums/", files[i]))
  phrase <- "Coverage for each repeat class and divergence (Kimura)"
  start.index <- match(phrase, divsum) + 1
  divsum <- divsum[start.index:length(divsum)]
  divsum <- read.table(textConnection(divsum), 
                       sep = " ", 
                       header = TRUE)
  divsum <- divsum[-c(which(sapply(divsum, function(col) all(is.na(col)))))]
  classes <- c("LINE", "SINE", "LTR", "DNA", "Div", "Unknown")
  for (j in classes) {
    pat <- paste0("^", j, "(\\.|$)")
    headers <- grep(pat, names(divsum), value = TRUE)
    sub <- divsum[, headers]
    sums <- rowSums(as.matrix(sub)) / asmbsz[i] * 100
    divsum <- divsum[, !names(divsum) %in% headers]
    assign(j, sums)
  }
  Others <- rowSums(as.matrix(divsum)) / asmbsz[i] * 100
  divsum <- data.frame(LINE, SINE, LTR, DNA, Others, Unknown)
  divsum <- t(as.matrix(divsum))
  divsum <- divsum[, 1:51]
  divsum <- divsum[nrow(divsum):1, ]
  barplot(divsum, 
          col = cols, 
          space = 0, 
          border = NA, 
          xlab = NA, 
          ylab = NA, 
          ylim = c(0, 1.1*3.049599))
  axis(1)
  pos <- seq(0.2 *3.049599, 0.85 * 3.049599, length.out = 6)
  text(55, pos[6], "LINE", xpd = NA, adj = c(0, 0.5))
  text(55, pos[5], "SINE", xpd = NA, adj = c(0, 0.5))
  text(55, pos[4], "LTR", xpd = NA, adj = c(0, 0.5))
  text(55, pos[3], "DNA", xpd = NA, adj = c(0, 0.5))
  text(55, pos[2], "Others", xpd = NA, adj = c(0, 0.5))
  text(55, pos[1], "Unclassified", xpd = NA, adj = c(0, 0.5))
  points(53, pos[6], pch = 15, col = cols[6], xpd = NA, cex = 1.5)
  points(53, pos[5], pch = 15, col = cols[5], xpd = NA, cex = 1.5)
  points(53, pos[4], pch = 15, col = cols[4], xpd = NA, cex = 1.5)
  points(53, pos[3], pch = 15, col = cols[3], xpd = NA, cex = 1.5)
  points(53, pos[2], pch = 15, col = cols[2], xpd = NA, cex = 1.5)
  points(53, pos[1], pch = 15, col = cols[1], xpd = NA, cex = 1.5)
}
par(mar = c(5.1, 4.1, 4.1, 2.1))





# phylogeny
library(phytools)
library(viridis)
dat <- read.csv("../results/parsed.csv")
dat <- dat[dat$thrs == 0.8, ]
dat <- dat[!duplicated(dat$species), ]
dat <- dat[, c("species", "clade", "rsq")]
tree <- read.tree("../data/formatted_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
int <- intersect(tree$tip.label, dat$species)
pruned.tree <- keep.tip(tree, int)
pruned.tree <- ladderize(pruned.tree)
dat <- dat[dat$species %in% int, ]
dat <- dat[match(pruned.tree$tip.label, dat$species), ]
rsq <- setNames(dat$rsq, dat$species)
p <- contMap(pruned.tree, rsq, plot = FALSE)
p <- setMap(p, viridis(n=300))
plot(p, 
     res = 300,
     lwd = c(1.5, 5),
     fsize = c(0.3, 0.5), 
     outline = FALSE,
     sig = 2,
     xlim = c(0, 115), 
     ylim = c(-150, 450),  #450
     direction = "downwards", 
     ftype = "off", 
     plot = TRUE)

