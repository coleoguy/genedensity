
## Zhaobo Hu
## zhaobohu2002@gmail.com

## Takes a species as input. drop contigs shorter than 10 million bp. drops
## assemblies that have less than 3 contigs or with p-values less than 0.05

# all contigs shorter than this size will be dropped
min.contig.size <- 10000000
# "vertebrates" or "invertebrates"?
vert.invert <- "vertebrates"

source("functions.R")
combined.results <- read.csv(paste0("../results/", vert.invert, "/combined_results.csv"))
species <- unique(combined.results$species)
results <- do.call(rbind, lapply(species, 
                                 parseResults, 
                                 combined.results = combined.results, 
                                 min.contig.size = min.contig.size))
write.csv(results, file = paste0("../results/", vert.invert, "/parsed_results.csv"), row.names = FALSE)
