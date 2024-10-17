

# load stuff in
packages <- c("ape", "ggplot2", "nlme")
lapply(packages, library, character.only = TRUE)
source("../analysis/functions.R")
dat <- read.csv("../data/vertebrates/chromnum_clade_gnsz_rsq.csv")
tree <- read.tree("../data/vertebrates/chordates_species.nwk")



# prune and format tree
sp <- dat$species
spf <- sub("^([^_]*_[^_]*)_.*", "\\1", gsub(" ", "_", sp))
sp.intersect <- intersect(tree$tip.label, spf)
pruned.tree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% sp.intersect)])
spn <- sp[match(spf[match(pruned.tree$tip.label, spf)], spf)]
pruned.tree$tip.label <- spn
# subset and reorder dataframe based on pruned tree data
dat <- dat[dat$species %in% spn, ]
dat <- dat[match(spn, dat$species), ]
# PGLS using nlme
nlme.pgls.model <- gls(rsq ~ gnsz_Gbp, 
                  data = dat, 
                  correlation = corBrownian(phy = pruned.tree, form = ~species),
                  method = "ML")
nlme.pgls.summary <- summary(nlme.pgls.model)
# graph
summary <- nlme.pgls.summary
intercept <- signif(summary$tTable[1, 1], 3)
slope <- signif(summary$tTable[2, 1], 3)
pval <- signif(summary$tTable[2, 4], 3)
dat$clade <- factor(dat$clade, levels = c("Mammalia", "Actinopterygii", "Sauria", "Others"))
ggplot(dat, aes(x = gnsz_Gbp, y = rsq, color = clade)) +
  geom_point(shape = 16, alpha = 0.4, size = 2.3) +
  scale_color_manual(labels = c(
    paste0("Mammals\n(n = ", sum(dat$clade == "Mammalia"), ")"),
    paste0("Ray-finned fish\n(n = ", sum(dat$clade == "Actinopterygii"), ")"), 
    paste0("Reptiles\n(n = ", sum(dat$clade == "Sauria"), ")"),
    paste0("Amphibians\n(n = ", sum(dat$clade == "Others"), ")")
  ), values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3"))+
  theme(plot.title = element_text(hjust = 0.475),
        plot.subtitle = element_text(hjust = 0.475), 
        axis.line = element_line(color = "black"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "#f2f2f2", color = "black", linewidth = 0.5),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "black", linetype = "dotted", size = 0.25),
        legend.position = c(0.87, 0.75),
        legend.key.size = unit(21, "points"))+
  geom_abline(intercept = intercept, slope = slope, color = "black", linetype = "dashed", linewidth = 0.5) +
  xlim(c(0.3, 6.3)) +
  ylim(c(0, 1)) +
  labs(title = bquote(italic(r)^2 ~ "vs Estimated Genome Size"), 
       subtitle = bquote(italic(β) * "-coefficient" == .(slope) * "," ~~ italic(p) * "-value" == .(pval)),
       x = "Estimated Genome Size (Gbp)", 
       y = bquote(italic(r)^2))
ggsave(filename = "vert_rsq_vs_gnsz_corrected.jpg", 
       plot = last_plot(), 
       width = 7680, 
       height = 4320, 
       units = "px", 
       dpi = 1100)

