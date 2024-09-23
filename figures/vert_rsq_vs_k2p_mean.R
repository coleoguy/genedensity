


library(ggplot2)
source("../analysis/functions.R")
parsed.results <- read.csv("../results/vertebrates/parsed_results.csv")
clades.gnsz <- read.csv("../data/vertebrates/clades_gnsz.csv")
files <- list.files("../results/vertebrates/repeat_landscape")
vert.invert <- "vertebrates"
k2p.mean <- sapply(files, getK2pMean)
species.lower <- gsub("_", " ", gsub("_summary\\.divsum$", "", files))
species <- gsub("^(\\w)(.*)", "\\U\\1\\L\\2", species.lower, perl = TRUE)
rsq <- sapply(species, getRsq)
class <- sapply(species, getClass)
custom.clade <- class
custom.clade[custom.clade == "Actinopterygii"] <- "Ray-finned fish"
custom.clade[custom.clade == "Aves"] <- "Reptiles"
custom.clade[custom.clade == "Mammalia"] <- "Mammals"
custom.clade[custom.clade == "Reptilia"] <- "Reptiles"
other.clades <- !custom.clade %in% c("Ray-finned fish", "Reptiles", "Mammals")
custom.clade[other.clades] <- "Others"
custom.clade <- factor(custom.clade, levels = c("Mammals", "Ray-finned fish", "Reptiles", "Others"))
rsq.vs.k2p.mean <- na.omit(data.frame(species, k2p.mean, rsq, class, custom.clade))
fit <- lm(rsq.vs.k2p.mean$rsq ~ rsq.vs.k2p.mean$k2p.mean)
slope <- signif(summary(fit)$coefficients[2, 1], 3)
intercept <- signif(summary(fit)$coefficients[1, 1], 3)
slope.pvalue <- signif(summary(fit)$coefficients[2, 4], 3)
slope.rsquared <- signif(summary(fit)$adj.r.squared, 3)

ggplot(rsq.vs.k2p.mean, aes(x = k2p.mean, y = rsq, color = custom.clade)) +
  geom_point(shape = 16, alpha = 0.4, size = 2.3) +
  scale_color_manual(labels = c(
    paste0("Mammals\n(n = ", sum(custom.clade == "Mammals"), ")"),
    paste0("Ray-finned fish\n(n = ", sum(custom.clade == "Ray-finned fish"), ")"), 
    paste0("Reptiles\n(n = ", sum(custom.clade == "Reptiles"), ")"),
    paste0("Others\n(n = ", sum(custom.clade == "Others"), ")")
  ), values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3"))+
  ggtitle(bquote(italic(r)^2~"vs Estimated Genome Size"))+
  theme(plot.title = element_text(hjust = 0.475),
        axis.line = element_line(color = "black"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "#f2f2f2", color = "black", linewidth = 0.5),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "black", linetype = "dotted", size = 0.25),
        legend.position = c(0.86, 0.69),
        legend.key.size = unit(21, "points"))+
  xlim(c(10, 26)) +
  ylim(c(0, 1))+
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5, fullrange = TRUE)+
  labs(x = "Mean K2P Distance (subs per site)", y = bquote(italic(r)^2))
ggsave(filename = "vert_rsq_vs_k2p_mean.jpg", 
       plot = last_plot(), 
       width = 7680, 
       height = 4320, 
       units = "px", 
       dpi = 1100)

