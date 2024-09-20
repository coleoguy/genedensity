

# load package
library(ggplot2)



parsed.results <- read.csv("../results/vertebrates/parsed_results.csv")
clades.gnsz <- read.csv("../data/vertebrates/clades_gnsz.csv")
table <- data.frame(parsed.results$species, parsed.results$species.rsquared)
table <- unique(table)
species <- table[[1]]
rsq <- table[[2]]
dat <- clades.gnsz[!is.na(clades.gnsz$class), ]
dat <- data.frame(dat$species, dat$class)
dat <- dat[dat[[1]] %in% species, ]
species <- dat[[1]]
class <- dat[[2]]
rsq <- parsed.results$species.rsquared
table <- data.frame(parsed.results$species, rsq)
table <- unique(na.omit(table))
table <- table[which(species %in% table[[1]]), ]
rsq <- table$rsq
custom.clade <- class
custom.clade[custom.clade == "Actinopterygii"] <- "Ray-finned fish"
custom.clade[custom.clade == "Aves"] <- "Reptiles"
custom.clade[custom.clade == "Mammalia"] <- "Mammals"
custom.clade[custom.clade == "Reptilia"] <- "Reptiles"
other.clades <- !custom.clade %in% c("Ray-finned fish", "Reptiles", "Mammals")
custom.clade[other.clades] <- "Others"
custom.clade <- factor(custom.clade, levels = c("Mammals", "Ray-finned fish", "Reptiles", "Others"))
clade.rsq <- data.frame(species, class, rsq, custom.clade)
num.mammals <- sum(custom.clade == "Mammals")
num.fish <- sum(custom.clade == "Ray-finned fish")
num.reptiles <- sum(custom.clade == "Reptiles")
num.others <- sum(custom.clade == "Others")
anova <- aov(rsq ~ custom.clade, data = clade.rsq)
pvalue <- (summary(anova)[[1]])[1, 5]
if (pvalue < 0.05) {
  tukey <- TukeyHSD(anova)
}
clade.rsq.filtered <- clade.rsq[clade.rsq$custom.clade != "Others", ]

ggplot(clade.rsq.filtered, aes(x = custom.clade, y = rsq, fill = custom.clade)) +
  ggtitle(bquote(italic(r)^2~"Across Clades"))+
  theme(plot.title = element_text(hjust = 0.45), 
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")),
        panel.grid.major = element_line(color = "black", linetype = "dotted", size = 0.25))+
  scale_x_discrete(labels = c("Mammals" = bquote("Mammals"~~(n==.(num.mammals))), 
                              "Ray-finned fish" = bquote("Ray-finned fish"~~(n==.(num.fish))), 
                              "Reptiles" = bquote("Reptiles"~~(n==.(num.reptiles))))) +
  labs(x = "", y = bquote(italic(r)^2)) +
  guides(fill = "none") +
  geom_violin() +
  geom_boxplot(width = 0.05)+
  geom_jitter(shape = 16, size = 1.5, position = position_jitter(0.23), alpha = 0.4, fill = "black", color = "black")
ggsave(filename = "rsq_across_clades.jpg", 
       plot = last_plot(), 
       width = 7680, 
       height = 4320, 
       units = "px", 
       dpi = 1100)
