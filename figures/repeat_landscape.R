
# verbose
verbose <- F

library(ggplot2)
files <- list.files(paste0("../results/vertebrates/repeat_landscape_plotting"))
parsed.results <- read.csv(paste0("../results/vertebrates/parsed.csv"))
parsed.results <- unique(parsed.results[1:30])
for (file in files) {
  dat <- read.csv(paste0("../results/vertebrates/repeat_landscape_plotting/", file))
  species <- sub("\\..*$", "", file)
  species <- gsub("_", " ", species)
  if (verbose == TRUE) {
    start.time <- Sys.time()
    gc()
    print(noquote(paste0(species, " (", Sys.time(), ")")))
  }
  wcv <- unique(parsed.results$weightcv[parsed.results$species == species])
  wcv <- signif(wcv, 3)
  ggplot(dat, aes(y = percent.cvrg, x = divergence, fill = repeat.class)) + 
    geom_bar(position="stack", stat="identity") +
    ggtitle(bquote("Repeat Landscape of"~italic(.(species))~~("WCV"==.(wcv)))) +
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
    ylim(c(0, 3)) +
    scale_fill_manual(
      breaks = c(unique(dat$repeat.class)[1:4], "Other", "Unknown"), 
      values = c("#e41a1c", "#ffff33", "#4daf4a", "#377eb8", "#ff7f00", "#984ea3")) +
    labs(x = "K2P Distance (Substitutions per Site)", 
         y = "Percent of Genome")
  ggsave(filename = paste0(gsub(" ", "_", species), ".jpg"), 
         plot = last_plot(), 
         path = paste0("repeat_landscape"), 
         width = 7680, 
         height = 4320, 
         units = "px", 
         dpi = 1100)
  if (verbose == TRUE) {
    print(noquote("   Saved!"))
  }
}

