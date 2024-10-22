



vert.invert <- "vertebrates"


# list packages
packages <- c("ggplot2", "ggrepel")
# load packages
lapply(packages, library, character.only = TRUE)
# pic loadings
loadings <- read.csv(paste0("../results/", vert.invert, "/final_results_pic_loadings.csv"), row.names = 1)
# Create biplot
ggplot(loadings, aes(x = PC1, y = PC2)) +
  # vectors
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
               # arrows
               arrow = arrow(), color = "blue") +  
  # vector labels
  geom_text_repel(data = loadings, aes(x = PC1, y = PC2, label = rownames(loadings)),
                  size = 4, nudge_x = 0.1, nudge_y = 0.1,
                  max.overlaps = Inf, force = 180) +
  # axis labels
  labs(# title = "PCA Biplot", 
       x = "Principal Component 1", 
       y = "Principal Component 2") +
  theme(plot.title = element_text(hjust = 0.475),
        plot.subtitle = element_text(hjust = 0.475), 
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "black", linetype = "dotted", size = 0.25),
        legend.position = "none")
ggsave(filename = paste0("all_dimensions_pca_", vert.invert, ".jpg"), 
       plot = last_plot(), 
       width = 7680, 
       height = 4320, 
       units = "px", 
       dpi = 1100)





