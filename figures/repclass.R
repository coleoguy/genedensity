
# get the biggest class of repeats for each species and see if any is 
# associated with higher or lower r2 values

terms <- c("line", "sine", "ltr", "dna", "rc")

dat <- read.csv("../results/parsed.csv")
dat <- dat[dat$thrs == 0.9, ]
df <- unique(dat[, c("species", "clade", "rsq")])
largest <- c()
for (i in df$species) {
  rep <- unique(dat[dat$species == i, ][, paste0(terms, ".rep.pct")])
  if (any(is.na(rep)) == TRUE) {
    largest <- c(largest, NA)
  } else {
    largest <- c(largest, strsplit(names(rep[which.max(rep)]), "\\.")[[1]][1])
  }
}
df$largest <- largest

library(beeswarm)
library(viridis)
cols <- viridis(4)
levels <- factor(as.factor(df$clade), levels = c("Others", "Actinopterygii", "Sauria", "Mammalia"))
par(mar = c(5, 4, 4, 7)+0.1)
beeswarm(rsq ~ largest, 
         data = df, 
         xlab = "largest repeat class", 
         ylab = "R2", 
         pch = 16, 
         pwcol = cols[levels])
text(5, (0.25+((0.8-0.25))), "Mammalia", xpd = NA, adj = 0)
text(5, (0.25+(2*(0.8-0.25)/3)), "Actinopterygii", xpd = NA, adj = 0)
text(5, (0.25+((0.8-0.25)/3)), "Sauria", xpd = NA, adj = 0)
text(5, 0.25, "Others", xpd = NA, adj = 0)
points(4.9, (0.25+((0.8-0.25))), pch = 16, col = "#FDE725FF", xpd = NA)
points(4.9, (0.25+(2*(0.8-0.25)/3)), pch = 16, col = "#31688EFF", xpd = NA)
points(4.9, (0.25+((0.8-0.25)/3)), pch = 16, col = "#35B779FF", xpd = NA)
points(4.9, 0.25, pch = 16, col = "#440154FF", xpd = NA)
par(mar = c(5.1, 4.1, 4.1, 2.1))
viridis(4)


# do clades differ in mean rsq after applying the new parsing step? NO
library(phytools)
tree <- read.tree("../data/formatted_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
int <- intersect(df$species, tree$tip.label)
pruned.tree <- keep.tip(tree, int)
rsq <- setNames(df$rsq, df$species)
signal <- phylosig(pruned.tree, rsq, method = "lambda", test = TRUE)[[4]]
if (signal < 0.05) {
  x <- setNames(df$clade, df$species)
  y <- setNames(df$rsq, df$species)
  phylANOVA(pruned.tree, x, y, nsim = 10000)
} else {
  summary(aov(rsq ~ clade, data = df))
}
