
# load functions
library(ggplot2)




parsed.results <- read.csv("../results/vertebrates/parsed_results.csv")
species <- unique(parsed.results$species)
rsq <- sapply(species, function(species) unique(parsed.results$species.rsquared[parsed.results$species == species]))
rsq <- data.frame(rsq)

# r-squared histogram
ggplot(rsq, aes(x = rsq)) +
  geom_histogram(bins = 10)+
  ggtitle(bquote("Histogram of"~italic(r)^2~"Across All Species")) +
  theme(plot.title = element_text(hjust = 0.475),
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "black", linetype = "dotted", size = 0.25)) +
  labs(x = bquote(italic(r)^2), y = "Frequency")
ggsave(filename = "rsq_histogram.jpg", 
       plot = last_plot(), 
       width = 7680, 
       height = 4320, 
       units = "px", 
       dpi = 1100)


