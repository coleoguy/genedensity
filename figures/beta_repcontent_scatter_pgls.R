

vert.invert <- "vertebrates"

# load stuff in
packages <- c("ape", "ggplot2", "nlme")
lapply(packages, library, character.only = TRUE)
source("../analysis/functions.R")
final.results <- read.csv(paste0("../results/", vert.invert, "/final_results.csv"))
tree <- read.tree(paste0("../data/", vert.invert, "/pruned_tree.nwk"))
tree$tip.label <- gsub("_", " ", tree$tip.label)

# gather and subset relevant results
dat <- na.omit(final.results[, c("species", "beta", "rep.content_percent.of.assembly", "clade")])
dat$beta <- dat$beta * 1000000
sp.intersect <- intersect(dat$species, tree$tip.label)
dat <- dat[dat$species %in% sp.intersect, ]

# prune tree
pruned.tree <- keep.tip(tree, sp.intersect)

# create PGLS object for trendline
pgls.model <- gls(beta ~ rep.content_percent.of.assembly, 
                  data = dat, 
                  correlation = corBrownian(phy = pruned.tree, form = ~species),
                  method = "ML")
summary <- summary(pgls.model)
intercept <- signif(summary$tTable[1, 1], 3)
slope <- signif(summary$tTable[2, 1], 3)

# calculate PICs for permutation test of pearson correlation coefficient
y <- pic(setNames(dat$beta, dat$species), pruned.tree)
x <- pic(setNames(dat$rep.content_percent.of.assembly, dat$species), pruned.tree)
pval <- signif(permTest(x, y, 1000000, "pearson"), 3)

# set factors for figure legend
dat$clade <- factor(dat$clade, levels = c("Mammalia", "Actinopterygii", "Sauria", "Others"))

# graph
ggplot(dat, aes(x = rep.content_percent.of.assembly, y = beta, color = clade)) +
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
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5, fullrange = TRUE)+
  labs(title = bquote(italic(β)~"vs Repeat Content"), 
       subtitle = bquote(italic(β) * "-coefficient" == .(slope) * "," ~~ "permutation test" ~ italic(p) * "-value" == .(pval)),
       x = "Repeat Content (% coverage)", 
       y = bquote(italic(β) ~ "(genes per Mbp)"))
ggsave(filename = paste0("beta_repcontent_scatter_pgls_", vert.invert, ".jpg"), 
       plot = last_plot(), 
       width = 7680, 
       height = 4320, 
       units = "px", 
       dpi = 1100)


