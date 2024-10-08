


library(ggplot2)
parsed.results <- read.csv("../results/vertebrates/parsed_results.csv")
clades.gnsz <- read.csv("../data/vertebrates/clades_gnsz.csv")
dat <- clades.gnsz[!is.na(clades.gnsz$genome.size.est_bp), ]
dat <- dat[!is.na(dat$class), ]
species <- dat$species
table <- unique(data.frame(parsed.results$species, parsed.results$species.rsquared))
table <- table[table[[1]] %in% species, ]
rsq <- table[[2]]
gnsz_bp <- dat$genome.size.est_bp
gnsz_Gbp <- gnsz_bp / 1000000000
class <- dat$class
custom.clade <- class
custom.clade[custom.clade == "Actinopterygii"] <- "Ray-finned fish"
custom.clade[custom.clade == "Aves"] <- "Reptiles"
custom.clade[custom.clade == "Mammalia"] <- "Mammals"
custom.clade[custom.clade == "Reptilia"] <- "Reptiles"
other.clades <- !custom.clade %in% c("Ray-finned fish", "Reptiles", "Mammals")
custom.clade[other.clades] <- "Others"
custom.clade <- factor(custom.clade, levels = c("Mammals", "Ray-finned fish", "Reptiles", "Others"))
rsq.vs.gnsz <- data.frame(species, rsq, gnsz_Gbp, class, custom.clade)
fit <- lm(rsq.vs.gnsz$rsq ~ rsq.vs.gnsz$gnsz_Gbp)
slope <- signif(summary(fit)$coefficients[2, 1], 3)
intercept <- signif(summary(fit)$coefficients[1, 1], 3)
fit.pval <- signif(summary(fit)$coefficients[2, 4], 3)
fit.rsq <- signif(summary(fit)$adj.r.squared, 3)

ggplot(rsq.vs.gnsz, aes(x = gnsz_Gbp, y = rsq, color = custom.clade)) +
  geom_point(shape = 16, alpha = 0.4, size = 2.3) +
  scale_color_manual(labels = c(
    paste0("Mammals\n(n = ", sum(rsq.vs.gnsz$custom.clade == "Mammals"), ")"),
    paste0("Ray-finned fish\n(n = ", sum(rsq.vs.gnsz$custom.clade == "Ray-finned fish"), ")"), 
    paste0("Reptiles\n(n = ", sum(rsq.vs.gnsz$custom.clade == "Reptiles"), ")"),
    paste0("Amphibians\n(n = ", sum(rsq.vs.gnsz$custom.clade == "Others"), ")")
  ), values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3"))+
  theme(plot.title = element_text(hjust = 0.475),
        plot.subtitle = element_text(hjust = 0.475), 
        axis.line = element_line(color = "black"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "#f2f2f2", color = "black", linewidth = 0.5),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "black", linetype = "dotted", size = 0.25),
        legend.position = c(0.86, 0.69),
        legend.key.size = unit(21, "points"))+
  xlim(c(0, 6.8)) +
  ylim(c(0.09, 1.02))+
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5, fullrange = TRUE)+
  labs(title = bquote(italic(r)^2~"vs Estimated Genome Size"), 
       subtitle = bquote(italic(β) * "-coefficient" == .(slope) * "," ~~ italic(p) * "-value" == .(fit.pval) * "," ~~ italic(r)^2 == .(fit.rsq)),
       x = "Estimated Genome Size (Gbp)", 
       y = bquote(italic(r)^2))
  ggsave(filename = "vert_rsq_vs_gnsz.jpg", 
       plot = last_plot(), 
       width = 7680, 
       height = 4320, 
       units = "px", 
       dpi = 1100)
