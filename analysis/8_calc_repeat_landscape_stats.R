# Zhaobo Hu
# zhaobohu2002@gmail.com

# Description: calculates repeat content, k

# "vertebrates" or "invertebrates"
vert.invert <- "invertebrates"

source("functions.R")
files <- list.files(paste0("../results/", 
                           vert.invert, 
                           "/repeat_landscape_divsums"))
asmbly.sz <- read.csv(paste0("../results/", vert.invert, "/assembly_sizes.csv"))


rep.landsc.stats <- do.call(rbind, lapply(files, 
                                          calcRepLandscapeStats, 
                                          asmbly.sz = asmbly.sz, 
                                          vert.invert = vert.invert))
write.csv(rep.landsc.stats, 
          paste0("../results/", vert.invert, "/rep_landscape_stats.csv"), 
          row.names = FALSE)
