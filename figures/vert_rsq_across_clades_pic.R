

packages <- c("ape", "ggplot2", "ggbeeswarm")
lapply(packages, library, character.only = TRUE)
source("../analysis/functions.R")
clades.gnsz <- read.csv("../data/vertebrates/clades_gnsz.csv")
parsed.results <- read.csv("../results/vertebrates/parsed_results.csv")
tree <- read.tree("../data/vertebrates/chordates_species.nwk")
species <- unique(parsed.results$species)
rsq <- as.numeric(lapply(species, getRsq))
gnsz_Gbp <- getGnszEst(species) / 1000000000
class <- getClass(species)
custom.clade <- class
custom.clade[custom.clade == "Actinopterygii"] <- "Ray-finned fish"
custom.clade[custom.clade == "Aves"] <- "Reptiles"
custom.clade[custom.clade == "Mammalia"] <- "Mammals"
custom.clade[custom.clade == "Reptilia"] <- "Reptiles"
others <- !(class %in% c("Actinopterygii", "Aves", "Mammalia", "Reptilia"))
custom.clade <- factor(custom.clade, levels = c("Mammals", "Ray-finned fish", "Reptiles", "Others"))
custom.clade[which(others)] <- "Others"
dat <- na.omit(data.frame(species, rsq, custom.clade))
rsq.pic <- getPIC(dat[1:2], tree)
rsq.clades <- data.frame(
  dat[dat$species %in% intersect(dat$species, names(rsq.pic)), ],
  rsq.pic)
rsq.clades <- rsq.clades[!(rsq.clades$custom.clade %in% "Others"), ]
num.mammals <- length(which(rsq.clades$custom.clade == "Mammals"))
num.fish <- length(which(rsq.clades$custom.clade == "Ray-finned fish"))
num.reptiles <- length(which(rsq.clades$custom.clade == "Reptiles"))
anova <- aov(rsq.clades$rsq.pic ~ rsq.clades$custom.clade)
aov.result <- summary(anova)[[1]][1, 5]
if (aov.result < 0.05) {
  tukey <- TukeyHSD(anova)
}

ggplot(rsq.clades, aes(x = custom.clade, y = rsq.pic, fill = custom.clade)) +
  ggtitle(bquote("PIC(" * italic(r)^2 * ") Across Clades"))+
  theme(plot.title = element_text(hjust = 0.45), 
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")),
        panel.grid.major = element_line(color = "black", linetype = "dotted", size = 0.25))+
  scale_x_discrete(labels = c("Mammals" = bquote("Mammals"~~(n==.(num.mammals))), 
                              "Ray-finned fish" = bquote("Ray-finned fish"~~(n==.(num.fish))), 
                              "Reptiles" = bquote("Reptiles"~~(n==.(num.reptiles))))) +
  labs(x = "", 
       y = bquote("PIC(" * italic(r)^2 * ")")) +
  guides(fill = "none") +
  geom_violin() +
  geom_boxplot(width = 0.05, outliers = FALSE) +
  ylim(c(-0.14, 0.07)) +
  geom_beeswarm(shape = 16, size = 1.5, cex = 1.75, alpha = 0.4, fill = "black", color = "black")
ggsave(filename = "vert_rsq_across_clades_pic.jpg", 
       plot = last_plot(), 
       width = 7680, 
       height = 4320, 
       units = "px", 
       dpi = 1100)
