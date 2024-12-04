

packages <- c("ggplot2", "ggbeeswarm", "phytools")
lapply(packages, library, character.only = TRUE)
dat <- read.csv("../results/vertebrates/parsed.csv")
dat <- unique(dat[1:34])
dat <- dat[!is.na(dat$chromnum.1n), ]
dat <- dat[, c("species", "rsq", "clade")]
dat <- na.omit(dat[dat$clade != "Others", ])
num.mammals <- sum(dat$clade == "Mammalia")
num.fish <- sum(dat$clade == "Actinopterygii")
num.reptiles <- sum(dat$clade == "Sauria")

# set factors for figure legend
dat$clade <- factor(dat$clade, levels = c("Mammalia", "Actinopterygii", "Sauria"))
# anova
aov <- aov(rsq ~ clade, data = dat)
anova <- summary(aov)[[1]][1, 5]
if (anova < 0.05) {
  tukey <- TukeyHSD(aov)$clade
}

# read and prune tree
tree <- read.tree("../data/vertebrates/formatted_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
sp.intersect <- intersect(dat$species, tree$tip.label)
datsub <- dat[dat$species %in% sp.intersect, ]
pruned.tree <- keep.tip(tree, sp.intersect)

# phylogenetic anova
anova <- phylANOVA(pruned.tree, setNames(datsub$clade, datsub$species), setNames(datsub$rsq, datsub$species))


ggplot(dat, aes(x = clade, y = rsq, fill = clade)) +
  ggtitle(bquote(italic(r)^2 ~ "Across Clades"))+
  theme(plot.title = element_text(hjust = 0.45), 
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")),
        panel.grid.major = element_line(color = "black", linetype = "dotted", size = 0.25))+
  scale_x_discrete(labels = c("Mammalia" = bquote("Mammals"~~(n==.(num.mammals))), 
                              "Actinopterygii" = bquote("Ray-finned fish"~~(n==.(num.fish))), 
                              "Sauria" = bquote("Reptiles"~~(n==.(num.reptiles))))) +
  labs(x = "", 
       y = bquote(italic(r)^2)) +
  guides(fill = "none") +
  geom_violin() +
  geom_boxplot(width = 0.05, outliers = FALSE) +
  ylim(c(0, 1)) +
  geom_beeswarm(shape = 16, size = 1.5, cex = 1.75, alpha = 0.4, fill = "black", color = "black")
ggsave(filename = paste0("rsq_clades_violin.jpg"), 
       plot = last_plot(), 
       width = 7680, 
       height = 4320, 
       units = "px", 
       dpi = 1100)
