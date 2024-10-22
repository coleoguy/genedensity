


## Zhaobo Hu
## zhaobohu2002@gmail.com

## Takes a species as input. drop contigs shorter than 10 million bp. drops
## assemblies with less than 3 contigs. 

# "vertebrates" or "invertebrates"?
vert.invert <- "invertebrates"
# all contigs shorter than this size will be dropped
min.contig.size <- 10000000

# source functions
source("functions.R")
# contig results for each species
combined.results <- read.csv(paste0("../results/", vert.invert, "/combined_results.csv"))
# unique species in results
species <- unique(combined.results$species)
# parse results
results <- do.call(rbind, lapply(species, 
                                 parseResults, 
                                 combined.results = combined.results, 
                                 min.contig.size = min.contig.size))
# write csv with parsed results
write.csv(results[1:5], 
          paste0("../results/", vert.invert, "/all_contigs_results.csv"), 
          row.names = FALSE)
# assembly sizes
asmbly.size <- unique(results[c(1, 6)])
# write csv with assembly sizes
write.csv(asmbly.size, 
          paste0("../results/", vert.invert, "/assembly_sizes.csv"), 
          row.names = FALSE)




