


packages <- c("ape", "ggplot2")
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
dat <- na.omit(data.frame(species, rsq, gnsz_Gbp, custom.clade))
rsq.pic <- getPIC(dat[1:2], tree)
gnsz.pic <- getPIC(dat[c(1, 3)], tree)
rsq.vs.gnsz <- data.frame(
  dat[dat$species %in% names(rsq.pic), ],
  rsq.pic,
  gnsz.pic)
fit <- summary(lm(rsq.vs.gnsz$rsq.pic ~ rsq.vs.gnsz$gnsz.pic))
slope <- signif(fit$coefficients[2, 1], 3)
intercept <- signif(fit$coefficients[1, 1], 3)
slope.pval <- signif(fit$coefficients[2, 4], 3)
fit.rsq <- signif(fit$adj.r.squared, 3)

ggplot(rsq.vs.gnsz, aes(x = gnsz.pic, y = rsq.pic, color = custom.clade)) +
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
        legend.position = c(0.84, 0.65),
        legend.key.size = unit(21, "points"))+
  xlim(c(-1.1, 1.4)) +
  ylim(c(-0.28, 0.15))+
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5, fullrange = TRUE)+
  labs(title = bquote("PIC(" * italic(r)^2 * ") vs PIC(Estimated Genome Size)"), 
       subtitle = bquote(italic(β) * "-coefficient" == .(slope) * "," ~~ italic(p) * "-value" == .(slope.pval) * "," ~~ italic(r)^2 == .(fit.rsq)),
       x = "PIC(Estimated Genome Size (Gbp))", 
       y = bquote("PIC(" * italic(r)^2 * ")"))
ggsave(filename = "vert_rsq_vs_gnsz_pic.jpg", 
       plot = last_plot(), 
       width = 7680, 
       height = 4320, 
       units = "px", 
       dpi = 1100)


