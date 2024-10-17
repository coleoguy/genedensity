

# "vertebrates" or "invertebrates"?
vert.invert <- "vertebrates"
# verbose
verbose <- T

library(ggplot2)
files <- list.files(paste0("../results/", vert.invert, "/repeat_landscape_plotting"))
parsed.results <- read.csv(paste0("../results/", vert.invert, "/parsed_Results.csv"))
for (file in files) {
  dat <- read.csv(paste0("../results/", vert.invert, "/repeat_landscape_plotting/", file))
  species <- sub("\\..*$", "", file)
  species <- gsub("_", " ", species)
  if (verbose == TRUE) {
    start.time <- Sys.time()
    gc()
    print(noquote(paste0(species, " (", Sys.time(), ")")))
  }
  rsq <- unique(parsed.results$species.rsquared[parsed.results$species == species])
  rsq <- signif(rsq, 3)
  ggplot(dat, aes(y = percent.of.genome, x = divergence, fill = repeat.group)) + 
    geom_bar(position="stack", stat="identity") +
    ggtitle(bquote("Repeat Landscape of"~italic(.(species))~~(italic(r)^2==.(rsq)))) +
    theme(plot.title = element_text(hjust = 0.45),
          axis.line = element_line(color = "black"),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "#f2f2f2", 
                                           color = "black", 
                                           linewidth = 0.5),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(color = "black", 
                                          linetype = "dotted", 
                                          size = 0.25),
          legend.position = c(0.884, 0.557),
          legend.key.size = unit(23.5, "points")) +
    scale_fill_manual(
      breaks = c(unique(dat$repeat.group)[1:4], "Other", "Unknown"), 
      values = c("#e41a1c", "#ffff33", "#4daf4a", "#377eb8", "#ff7f00", "#984ea3")) +
    labs(x = "K2P Distance (Substitutions per Site)", 
         y = "% of Genome Size Occupied by Repeats")
  ggsave(filename = paste0(gsub(" ", "_", species), ".jpg"), 
         plot = last_plot(), 
         path = paste0("../figures/", vert.invert, "/repeat_landscape"), 
         width = 7680, 
         height = 4320, 
         units = "px", 
         dpi = 1100)
  if (verbose == TRUE) {
    print(noquote("   Saved!"))
  }
}

