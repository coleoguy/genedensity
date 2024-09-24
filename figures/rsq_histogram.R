
# "vertebrates" or "invertebrates"? 
vert.invert <- "vertebrates"

library(ggplot2)
parsed.results <- read.csv(paste0("../results/", vert.invert, "/parsed_results.csv"))
species <- unique(parsed.results$species)
rsq <- sapply(species, function(species) unique(parsed.results$species.rsquared[parsed.results$species == species]))
rsq <- data.frame(rsq)
vert.invert <- paste0(toupper(substr(vert.invert, 1, 1)), 
                       substr(vert.invert, 2, nchar(vert.invert) - 1))
ggplot(rsq, aes(x = rsq)) +
  geom_histogram(bins = 10)+
  ggtitle(bquote("Histogram of"~italic(r)^2~"Across All"~.(vert.invert)~"Species")) +
  theme(plot.title = element_text(hjust = 0.475),
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "black", linetype = "dotted", size = 0.25)) +
  labs(x = bquote(italic(r)^2), y = "Frequency")
ggsave(filename = paste0(regmatches(vert.invert,
                                    regexpr(paste0(".*?", "vert"), vert.invert)), 
                         "_rsq_histogram.jpg"), 
       plot = last_plot(), 
       width = 7680, 
       height = 4320, 
       units = "px", 
       dpi = 1100)


