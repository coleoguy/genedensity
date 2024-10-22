

vert.invert <- "vertebrates"

packages <- c("ape", "ggplot2", "nlme")
lapply(packages, library, character.only = TRUE)
final.results <- read.csv(paste0("../results/", vert.invert, "/final_results.csv"))
pruned.tree <- read.tree(paste0("../data/", vert.invert, "/pruned_tree.nwk"))
pruned.tree$tip.label <- gsub("_", " ", pruned.tree$tip.label)
# subset relevant data for graphing
dat <- na.omit(final.results[, c("species", "rsq", "k2p.mean", "clade")])
# subset data for pgls line
dat.sbset <- dat[dat$species %in% pruned.tree$tip.label, ]
# PGLS using nlme
pgls.model <- gls(rsq ~ k2p.mean, 
                  data = dat.sbset, 
                  correlation = corBrownian(phy = pruned.tree, form = ~species),
                  method = "ML")
# graph
summary <- summary(pgls.model)
intercept <- signif(summary$tTable[1, 1], 3)
slope <- signif(summary$tTable[2, 1], 3)
pval <- signif(summary$tTable[2, 4], 3)
dat$clade <- factor(dat$clade, levels = c("Mammalia", "Actinopterygii", "Sauria", "Others"))
ggplot(dat, aes(x = k2p.mean, y = rsq, color = clade)) +
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
        legend.background = element_rect(fill = "#f2f2f2", color = "black", linewidth = 0.5),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "black", linetype = "dotted", size = 0.25),
        legend.position = c(0.86, 0.69),
        legend.key.size = unit(21, "points"))+
  xlim(c(10, 26)) +
  ylim(c(0, 1))+
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5, fullrange = TRUE)+
  labs(# title = bquote(italic(r)^2~"vs Estimated Genome Size"), 
       # subtitle = bquote(italic(β) * "-coefficient" == .(slope) * "," ~~ italic(p) * "-value" == .(pval)),
       x = "Mean K2P Distance (Subs per Site)", 
       y = bquote(italic(r)^2))
ggsave(filename = paste0("rsq_k2pmean_scatter_pgls_", vert.invert, ".jpg"), 
       plot = last_plot(), 
       width = 7680, 
       height = 4320, 
       units = "px", 
       dpi = 1100)

