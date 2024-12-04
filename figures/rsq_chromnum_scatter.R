





library(ggplot2)
library(ggbeeswarm)
library(ape)
library(nlme)
dat <- read.csv("../results/vertebrates/parsed.csv")
dat <- dat[!duplicated(dat$species), ]
dat <- dat[dat$clade %in% c("Actinopterygii", "Mammalia", "Sauria"), ]
dat <- na.omit(dat[, c(which(colnames(dat) %in% c("species", "clade", "chromnum.1n", "rsq")))])
dat <- dat[!is.na(dat$chromnum.1n), ]


tree <- read.tree("../data/vertebrates/formatted_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
sp.intersect <- intersect(dat$species, tree$tip.label)
dat <- dat[dat$species %in% sp.intersect, ]
pruned.tree <- keep.tip(tree, sp.intersect)


# create PGLS object for trendline
pgls.model <- gls(rsq ~ chromnum.1n, 
                  data = dat, 
                  correlation = corBrownian(phy = pruned.tree, form = ~species),
                  method = "ML")
summary <- summary(pgls.model)
intercept <- signif(summary$tTable[1, 1], 3)
slope <- signif(summary$tTable[2, 1], 3)
slope.pval <- signif(summary$tTable[2, 4], 3)

# pgls for mammals
mam <- dat[dat$clade %in% "Mammalia", ]
mam.model <- gls(rsq ~ chromnum.1n, 
                 data = mam, 
                 correlation = corBrownian(phy = pruned.tree, form = ~species),
                 method = "ML")
summary <- summary(mam.model)
mam.int <- signif(summary$tTable[1, 1], 3)
mam.slope <- signif(summary$tTable[2, 1], 3)
mam.slope.p <- signif(summary$tTable[2, 4], 3)


# pgls for fish
fish <- dat[dat$clade %in% "Actinopterygii", ]
fish.model <- gls(rsq ~ chromnum.1n, 
                  data = fish, 
                  correlation = corBrownian(phy = pruned.tree, form = ~species),
                  method = "ML")
summary <- summary(fish.model)
fish.int <- signif(summary$tTable[1, 1], 3)
fish.slope <- signif(summary$tTable[2, 1], 3)
fish.slope.p <- signif(summary$tTable[2, 4], 3)

# pgls for reptiles
rep <- dat[dat$clade %in% "Sauria", ]
rep.model <- gls(rsq ~ chromnum.1n, 
                 data = rep, 
                 correlation = corBrownian(phy = pruned.tree, form = ~species),
                 method = "ML")
summary <- summary(rep.model)
rep.int <- signif(summary$tTable[1, 1], 3)
rep.slope <- signif(summary$tTable[2, 1], 3)
rep.slope.p <- signif(summary$tTable[2, 4], 3)

dat$clade <- factor(dat$clade, levels = c("Mammalia", "Actinopterygii", "Sauria"))

# graph
ggplot(dat, aes(x = chromnum.1n, y = rsq, color = clade)) +
  geom_point(shape = 16, alpha = 0.4, size = 2.3) +
  scale_color_manual(labels = c(
    paste0("Mammals\n(n = ", sum(dat$clade == "Mammalia"), ")"),
    paste0("Ray-finned fish\n(n = ", sum(dat$clade == "Actinopterygii"), ")"),
    paste0("Reptiles\n(n = ", sum(dat$clade == "Sauria"), ")")
  ), values = c("#e41a1c", "#377eb8", "#4daf4a"))+
  theme(plot.title = element_text(hjust = 0.475),
        plot.subtitle = element_text(hjust = 0.475), 
        axis.line = element_line(color = "black"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "black", linetype = "dotted", size = 0.25),
        legend.position = c(0.87, 0.35),
        legend.key.size = unit(21, "points"),
        legend.margin = margin(r = 5, l = 1, t = 2, b = 4))+
  xlim(c(7, 59)) +
  geom_abline(intercept = intercept, slope = slope, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_abline(intercept = mam.int, slope = mam.slope, color = "#e41a1c", linetype = "dashed", linewidth = 0.5) +
  geom_abline(intercept = fish.int, slope = fish.slope, color = "#377eb8", linetype = "dashed", linewidth = 0.5) +
  geom_abline(intercept = rep.int, slope = rep.slope, color = "#4daf4a", linetype = "dashed", linewidth = 0.5) +
  labs(title = bquote(italic(R)^2 ~ "vs Chromosome Count"), 
       x = bquote("Chromosome Count"), 
       y = bquote(italic(R)^2))
ggsave(filename = paste0("rsq_chromnum_scatter.jpg"), 
       plot = last_plot(), 
       width = 7680, 
       height = 4320, 
       units = "px", 
       dpi = 1100)

