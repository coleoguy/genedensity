



# load stuff in
packages <- c("ape", "ggplot2", "nlme")
lapply(packages, library, character.only = TRUE)
source("../analysis/functions.R")
results <- read.csv("../results/vertebrates/parsed.csv")
tree <- read.tree(paste0("../data/vertebrates/formatted_tree.nwk"))
tree$tip.label <- gsub("_", " ", tree$tip.label)

dat <- results[!duplicated(results$species), ]
dat <- na.omit(dat[, c("species", "clade", "rsq")])
meandens <- c()
for (i in unique(results$species)) {
  sub <- subset(results, species == i)
  meandens <- c(meandens, mean(sub$genedens))
}
dat$meandens <- meandens
dat$normrsq <- dat$rsq / dat$meandens

sp.intersect <- intersect(dat$species, tree$tip.label)
dat <- dat[dat$species %in% sp.intersect, ]
pruned.tree <- keep.tip(tree, sp.intersect)

# create PGLS object for trendline
pgls.model <- gls(rsq ~ meandens, 
                  data = dat, 
                  correlation = corBrownian(phy = pruned.tree, form = ~species),
                  method = "ML")
summary <- summary(pgls.model)
slope <- signif(summary$tTable[2, 1], 3)
intercept <- signif(summary$tTable[1, 1], 3)
slope.pval <- signif(summary$tTable[2, 4], 3)


# pgls for mammals
mam <- dat[dat$clade %in% "Mammalia", ]
mam.model <- gls(rsq ~ meandens, 
                 data = mam, 
                 correlation = corBrownian(phy = pruned.tree, form = ~species),
                 method = "ML")
summary <- summary(mam.model)
mam.int <- signif(summary$tTable[1, 1], 3)
mam.slope <- signif(summary$tTable[2, 1], 3)
mam.slope.p <- signif(summary$tTable[2, 4], 3)


# pgls for fish
fish <- dat[dat$clade %in% "Actinopterygii", ]
fish.model <- gls(rsq ~ meandens, 
                  data = fish, 
                  correlation = corBrownian(phy = pruned.tree, form = ~species),
                  method = "ML")
summary <- summary(fish.model)
fish.int <- signif(summary$tTable[1, 1], 3)
fish.slope <- signif(summary$tTable[2, 1], 3)
fish.slope.p <- signif(summary$tTable[2, 4], 3)

# pgls for reptiles
rep <- dat[dat$clade %in% "Sauria", ]
rep.model <- gls(rsq ~ meandens, 
                 data = rep, 
                 correlation = corBrownian(phy = pruned.tree, form = ~species),
                 method = "ML")
summary <- summary(rep.model)
rep.int <- signif(summary$tTable[1, 1], 3)
rep.slope <- signif(summary$tTable[2, 1], 3)
rep.slope.p <- signif(summary$tTable[2, 4], 3)

# set factors for figure legend
dat$clade <- factor(dat$clade, levels = c("Mammalia", "Actinopterygii", "Sauria", "Others"))

# graph
ggplot(dat, aes(x = meandens, y = rsq, color = clade)) +
  coord_trans(x = "log10") +
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
        legend.position = c(0.6, 0.2),
        legend.key.size = unit(21, "points"),
        legend.margin = margin(r = 5, l = 1, t = 2, b = 4)) +
  geom_abline(intercept = intercept, slope = slope, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_abline(intercept = mam.int, slope = mam.slope, color = "#e41a1c", linetype = "dashed", linewidth = 0.5) +
  geom_abline(intercept = fish.int, slope = fish.slope, color = "#377eb8", linetype = "dashed", linewidth = 0.5) +
  geom_abline(intercept = rep.int, slope = rep.slope, color = "#4daf4a", linetype = "dashed", linewidth = 0.5) +
  labs(title = bquote(italic(R)^2 ~ "vs Mean Gene Density"), 
       subtitle = bquote(italic(β) == .(slope) * "," ~~ italic(β) ~ italic(p) * "-value" == .(slope.pval)),
       x = bquote("Mean Gene Density (Genes per Mb)"), 
       y = bquote(italic(R)^2)) +
  xlim(c(min(meandens), 70))
ggsave(filename = paste0("rsq_meandens_scatter.jpg"), 
       plot = last_plot(), 
       width = 7680, 
       height = 4320, 
       units = "px", 
       dpi = 1100)
