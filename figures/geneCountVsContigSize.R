

# constants
requiredPackages <- c("data.table", "ggplot2")

# paths
functionsPath <- "../../analysis/functions.R"
parsedResultsCsvPath <- "../../results/vertebrates/parsedResults.csv"

# load stuff in
source(functionsPath)
loadPackages(requiredPackages)
dat <- fread(parsedResultsCsvPath)


species <- unique(dat$species)

for (i in species) {
  geneCount <- dat$contigGeneCount[dat$species == i]
  contigSize <- dat$contigSize.bp[dat$species == i]
  geneCountVsContigSize <- data.table(geneCount, contigSize)
  fit <- summary(lm(geneCount ~ contigSize))
  slope <- signif(fit$coefficients[2, 1], 3)
  intercept <- signif(fit$coefficients[1, 1], 3)
  rsq <- signif(fit$adj.r.squared, 3)
  ggplot(geneCountVsContigSize, aes(x = contigSize, y = geneCount)) +
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
         path = "vertebrates",
         plot = last_plot(), 
         width = 7680, 
         height = 4320, 
         units = "px", 
         dpi = 1100)
}




