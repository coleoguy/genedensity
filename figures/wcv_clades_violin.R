


packages <- c("ggplot2", "ggbeeswarm")
lapply(packages, library, character.only = TRUE)
results <- read.csv("../results/vertebrates/parsed.csv")
results <- unique(results[1:22])
results <- results[!is.na(results$chromnum.1n), ]
dat <- results[, c("species", "weightcv", "clade")]
dat <- na.omit(dat[dat$clade != "Others", ])
num.mammals <- sum(dat$clade == "Mammalia")
num.fish <- sum(dat$clade == "Actinopterygii")
num.reptiles <- sum(dat$clade == "Sauria")

# set factors for figure legend
dat$clade <- factor(dat$clade, levels = c("Mammalia", "Actinopterygii", "Sauria"))


aovresult <- aov(weightcv ~ clade, data = dat)
anova <- summary(aovresult)[[1]][1, 5]
if (anova < 0.05) {
  tukey <- TukeyHSD(aovresult)$clade
}

ggplot(dat, aes(x = clade, y = weightcv, fill = clade)) +
  ggtitle(bquote("Coefficient of Variation Across Clades"))+
  theme(plot.title = element_text(hjust = 0.45), 
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")),
        panel.grid.major = element_line(color = "black", linetype = "dotted", size = 0.25))+
  scale_x_discrete(labels = c("Mammalia" = bquote("Mammals"~~(n==.(num.mammals))), 
                              "Actinopterygii" = bquote("Ray-finned fish"~~(n==.(num.fish))), 
                              "Sauria" = bquote("Reptiles"~~(n==.(num.reptiles))))) +
  labs(x = "", 
       y = bquote("Coefficient of Variation")) +
  guides(fill = "none") +
  geom_violin() +
  geom_boxplot(width = 0.05, outliers = FALSE) +
  ylim(c(0.04, 0.76)) +
  geom_beeswarm(shape = 16, size = 1.5, cex = 1.75, alpha = 0.4, fill = "black", color = "black")
ggsave(filename = paste0("wcv_clades_violin.jpg"), 
       plot = last_plot(), 
       width = 7680, 
       height = 4320, 
       units = "px", 
       dpi = 1100)
