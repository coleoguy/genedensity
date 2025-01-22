
rep <- read.csv("../results/vertebrates/all_repeat_landscapes.csv")
sp <- unique(rep$species)
df <- data.frame()
for (i in sp) {
  sub <- rep[rep$species == i, ]
  group <- unique(sub$repeat.class)
  group <- group[!(group %in% c("Other", "Unknown"))]
  for (j in group) {
    subsub <- sub[sub$repeat.class == j, ]
    median.bin <- which(cumsum(subsub$percent.cvrg) > sum(subsub$percent.cvrg)/2)[1]
    lower <- cumsum(subsub$percent.cvrg)[median.bin-1]
    if (median.bin == 1) {
      lower <- 0
    }
    upper <- cumsum(subsub$percent.cvrg)[median.bin+1]
    mid <- sum(subsub$percent.cvrg)/2
    mediandvg <- median.bin + (mid-lower)/(upper-lower)
    species <- i
    repeatclass <- j
    df <- rbind(df, data.frame(species, mediandvg, repeatclass))
  }
}
top4 <- head(names(sort(table(df$repeatclass), decreasing = TRUE)), 4)
dat <- read.csv("../results/vertebrates/parsed.csv")
dat <- dat[!duplicated(dat$species), ]

first <- setNames(df[df$repeatclass == top4[1], 1:2], c("species", top4[1]))
second <- setNames(df[df$repeatclass == top4[2], 1:2], c("species", top4[2]))
third <- setNames(df[df$repeatclass == top4[3], 1:2], c("species", top4[3]))
last <- setNames(df[df$repeatclass == top4[4], 1:2], c("species", top4[4]))

dat <- merge(dat, first, by = "species", all.x = TRUE)
dat <- merge(dat, second, by = "species", all.x = TRUE)
dat <- merge(dat, third, by = "species", all.x = TRUE)
dat <- merge(dat, last, by = "species", all.x = TRUE)





# load stuff in
packages <- c("ape", "ggplot2", "nlme")
lapply(packages, library, character.only = TRUE)
source("../analysis/functions.R")
tree <- read.tree(paste0("../data/vertebrates/formatted_tree.nwk"))
tree$tip.label <- gsub("_", " ", tree$tip.label)

# gather and subset relevant results
dat <- dat[!is.na(dat$chromnum.1n), ]
dat <- na.omit(dat[, c("species", "rsq", "DNA", "clade")])
sp.intersect <- intersect(dat$species, tree$tip.label)
dat <- dat[dat$species %in% sp.intersect, ]

# prune tree
pruned.tree <- keep.tip(tree, sp.intersect)

# create PGLS object for trendline
pgls.model <- gls(rsq ~ DNA, 
                  data = dat, 
                  correlation = corBrownian(phy = pruned.tree, form = ~species),
                  method = "ML")
summary <- summary(pgls.model)
intercept <- signif(summary$tTable[1, 1], 3)
slope <- signif(summary$tTable[2, 1], 3)
slope.p <- signif(summary$tTable[2, 4], 3)

# pgls for mammals
mam <- dat[dat$clade %in% "Mammalia", ]
mam.model <- gls(rsq ~ DNA, 
                 data = mam, 
                 correlation = corBrownian(phy = pruned.tree, form = ~species),
                 method = "ML")
summary <- summary(mam.model)
mam.int <- signif(summary$tTable[1, 1], 3)
mam.slope <- signif(summary$tTable[2, 1], 3)
mam.slope.p <- signif(summary$tTable[2, 4], 3)


# pgls for fish
fish <- dat[dat$clade %in% "Actinopterygii", ]
fish.model <- gls(rsq ~ DNA, 
                  data = fish, 
                  correlation = corBrownian(phy = pruned.tree, form = ~species),
                  method = "ML")
summary <- summary(fish.model)
fish.int <- signif(summary$tTable[1, 1], 3)
fish.slope <- signif(summary$tTable[2, 1], 3)
fish.slope.p <- signif(summary$tTable[2, 4], 3)

# pgls for reptiles
rep <- dat[dat$clade %in% "Sauria", ]
rep.model <- gls(rsq ~ DNA, 
                 data = rep, 
                 correlation = corBrownian(phy = pruned.tree, form = ~species),
                 method = "ML")
summary <- summary(rep.model)
rep.int <- signif(summary$tTable[1, 1], 3)
rep.slope <- signif(summary$tTable[2, 1], 3)
rep.slope.p <- signif(summary$tTable[2, 4], 3)


# calculate PICs for permutation test of pearson correlation coefficient
y <- pic(setNames(dat$rsq, dat$species), pruned.tree)
x <- pic(setNames(dat$DNA, dat$species), pruned.tree)
perm.p <- signif(permTest(x, y, 100000, "pearson"), 3)

# set factors for figure legend
dat$clade <- factor(dat$clade, levels = c("Mammalia", "Actinopterygii", "Sauria", "Others"))

# graph
ggplot(dat, aes(x = DNA, y = rsq, color = clade)) +
  geom_point(shape = 16, alpha = 0.4, size = 2.3) +
  scale_color_manual(labels = c(
    paste0("Mammals\n(n = ", sum(dat$clade == "Mammalia"), ")"),
    paste0("Ray-finned fish\n(n = ", sum(dat$clade == "Actinopterygii"), ")"), 
    paste0("Reptiles\n(n = ", sum(dat$clade == "Sauria"), ")"),
    paste0("Others\n(n = ", sum(dat$clade == "Others"), ")")
  ), values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3"))+
  theme(plot.title = element_text(hjust = 0.475),
        plot.subtitle = element_text(hjust = 0.475), 
        axis.line = element_line(color = "black"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "black", linetype = "dotted", size = 0.25),
        legend.position = c(0.86, 0.69),
        legend.key.size = unit(21, "points"),
        legend.margin = margin(r = 5, l = 1, t = 2, b = 4))+
  #geom_abline(intercept = intercept, slope = slope, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_abline(intercept = mam.int, slope = mam.slope, color = "#e41a1c", linetype = "dashed", linewidth = 0.5) +
  geom_abline(intercept = fish.int, slope = fish.slope, color = "#377eb8", linetype = "dashed", linewidth = 0.5) +
  geom_abline(intercept = rep.int, slope = rep.slope, color = "#4daf4a", linetype = "dashed", linewidth = 0.5) +
  labs(title = bquote(italic(r)^2~"vs Median Divergence"), 
       subtitle = bquote(italic(β) * "-coefficient" == .(slope) * "," ~~ italic(β) ~ italic(p) * "-value" == .(slope.p) * "," ~~ "permutation" ~ italic(p) * "-value" == .(perm.p)),
       x = "Median Divergence", 
       y = bquote(italic(r)^2))
