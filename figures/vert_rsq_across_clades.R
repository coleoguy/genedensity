

packages <- c("ape", "ggplot2", "ggbeeswarm")
lapply(packages, library, character.only = TRUE)
source("../analysis/functions.R")
dat <- read.csv("../data/vertebrates/chromnum_clade_gnsz_rsq.csv")
dat <- dat[dat$clade %in% c("Mammalia", "Actinopterygii", "Sauria"), ]
num.mammals <- sum(dat$clade == "Mammalia")
num.fish <- sum(dat$clade == "Actinopterygii")
num.reptiles <- sum(dat$clade == "Sauria")
ggplot(dat, aes(x = clade, y = rsq, fill = clade)) +
  ggtitle(bquote(italic(r)^2 * "Across Clades"))+
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
ggsave(filename = "vert_rsq_across_clades.jpg", 
       plot = last_plot(), 
       width = 7680, 
       height = 4320, 
       units = "px", 
       dpi = 1100)
