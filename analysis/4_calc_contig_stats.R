


# Zhaobo Hu
# zhaobohu2002@gmail.com

# Description: for each species, fit a linear model for parsed contigs and 
# calculate slope, p-value of slope, r-squared, and coefficient of varience

# "vertebrates" or "invertebrates"
vert.invert <- "invertebrates"

# source functions
source("functions.R")
# read results
all.contigs.results <- read.csv(paste0("../results/", 
                                       vert.invert, 
                                       "/all_contigs_results.csv"))
# list species in results
species <- unique(all.contigs.results$species)
# calculate stats
contig.stats <- do.call(rbind, lapply(species, 
                                      getContigStats, 
                                      results = all.contigs.results))
# add species to stats dataframe
contig.stats <- data.frame(species, contig.stats)
# write csv
write.csv(contig.stats, 
          paste0("../results/", vert.invert, "/contig_stats.csv"), 
          row.names = FALSE)



