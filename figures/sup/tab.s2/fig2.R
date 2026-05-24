# beeswarm plot for each predictor
# also do tests

library(beeswarm)

dat <- read.csv("../../../results/rsq.csv")
clades <- c("Mammalia", "Actinopterygii", "Sauropsida")
dat <- dat[dat$clade %in% clades, ]
dat$clade <- factor(dat$clade, levels = clades)
map <- setNames(c("#d95f02", "#7570b3", "#1b9e77"), clades)
cols <- map[dat$clade]
num.sp <- sapply(clades, function(cl) nrow(dat[dat$clade == cl, ]))

# flip gini and cv so higher = more homogeneous, matching main figures
dat$gini.flip <- -dat$gini
dat$cv.flip <- -dat$cv

responses <- c("rsq", "gini.flip", "cv.flip")
resp.names <- c("rsq", "gini", "cv")
resp.ylabs <- list(
  expression(R^2), 
  "-(Gini)", 
  "-(CV)"
)

# plot
par(mfrow = c(1, 1))
for (ri in seq_along(responses)) {
  resp <- responses[ri]

  par(mar = c(5, 5, 2, 2) + 0.1)

  beeswarm(dat[[resp]] ~ dat$clade, 
           xlab = NA, ylab = NA, 
           pwcol = cols, pch = 16, 
           spacing = 1.4, cex = 1.2, 
           xaxt = "n", yaxt = "n")

  axis(side = 1, at = seq(num.sp), 
       labels = c("Mammals", "Ray-finned fish", "Reptiles"), 
       mgp = c(3, 0.9, 0))
  axis(side = 1, at = seq(num.sp), 
       labels = paste0("n=", num.sp), mgp = c(3, 2.1, 0))
  axis(side = 2, mgp = c(3, 0.9, 0))
  mtext(resp.ylabs[[ri]], side = 2, line = 2.4)
}




# phylo vs non-phylo
library(phytools)
tree <- read.tree("../../../data/formatted-tree.nwk")
tree$tip.label <- tolower(tree$tip.label)
int <- intersect(dat$species, tree$tip.label)
pruned.tree <- keep.tip(tree, int)
sub.dat <- dat[dat$species %in% int, ]

anova.results <- data.frame()

for (i in seq_along(responses)) {
  resp <- responses[i]
  
  # standard anova
  aov.fit <- aov(sub.dat[[resp]] ~ sub.dat$clade)
  aov.summ <- summary(aov.fit)[[1]]
  
  # tukey hsd
  tukey <- TukeyHSD(aov.fit)
  tukey.df <- as.data.frame(tukey$`sub.dat$clade`)
  tukey.df$comparison <- rownames(tukey.df)
  
  # phylo anova
  y <- setNames(sub.dat[[resp]], sub.dat$species)
  x <- setNames(sub.dat$clade, sub.dat$species)
  pa <- phylANOVA(pruned.tree, x, y, nsim = 1000, posthoc = FALSE)
  
  for (j in seq_len(nrow(tukey.df))) {
    anova.results <- rbind(anova.results, data.frame(
      response = resp, 
      anova.F = aov.summ[["F value"]][1], 
      anova.P = aov.summ[["Pr(>F)"]][1], 
      tukey.comparison = tukey.df$comparison[j], 
      tukey.diff = tukey.df$diff[j], 
      tukey.lwr = tukey.df$lwr[j], 
      tukey.upr = tukey.df$upr[j], 
      tukey.P = tukey.df$`p adj`[j], 
      phylo.F = pa$F, 
      phylo.P = pa$Pf
    ))
  }
}

write.csv(anova.results, "anova-results.csv", row.names = FALSE)
