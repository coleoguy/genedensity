
# load packages
library(data.table)
library(ggplot2)

dat <- fread("../results/vertebrates/parsed_results.csv")
species <- unique(dat$species)

for (i in species) {
  gene.count <- dat$contig.gene.count[dat$species == i]
  contig.size <- dat$contig.size_bp[dat$species == i]
  gene.count.vs.contig.size <- data.table(gene.count, contig.size)
  fit <- summary(lm(gene.count ~ contig.size))
  slope <- signif(fit$coefficients[2, 1], 3)
  intercept <- signif(fit$coefficients[1, 1], 3)
  rsq <- signif(fit$adj.r.squared, 3)
  ggplot(gene.count.vs.contig.size, aes(x = contig.size, y = gene.count)) +
    geom_point(shape = 16, alpha = 0.4, size = 2.3) +
    theme(plot.title = element_text(hjust = 0.475), 
          plot.subtitle = element_text(hjust = 0.475), 
          axis.line = element_line(color = "black"),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(color = "black", linetype = "dotted", size = 0.25)) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5, fullrange = TRUE)+
    labs(
      title = bquote("Gene Count vs Contig Size for"~italic(.(i))), 
      subtitle = bquote(y == .(slope) * x + .(intercept) * "," ~~ italic(r)^2 == .(rsq)),
      x = "Contig Size", 
      y = "Gene Count")
  ggsave(filename = paste0(i, ".jpg"), 
         path = "vertebrates/gene_count_vs_contig_size",
         plot = last_plot(), 
         width = 7680, 
         height = 4320, 
         units = "px", 
         dpi = 1100)
}