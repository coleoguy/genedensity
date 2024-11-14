


## Zhaobo Hu
## zhaobohu2002@gmail.com

## Takes a species as input. drop contigs shorter than 10 million bp. drops
## assemblies with less than 3 contigs. 

# "vertebrates" or "invertebrates"?
vert.invert <- "vertebrates"

# source functions
source("functions.R")
# contig results for each species
results <- read.csv(paste0("../results/", vert.invert, "/combined_results.csv"))
# unique species in results
species <- unique(results$species)

parsed <- data.frame()
# parse results
for (i in species) {
  sub <- results[results$species == i, ]
  sub <- sub[sub$contig.size.Mbp > 10, ]
  if (nrow(sub) >= 3) {
    parsed <- rbind(parsed, sub)
  }
}

# write csv with parsed results
write.csv(parsed[1:5], 
          paste0("../results/", vert.invert, "/all_contigs_results.csv"), 
          row.names = FALSE)
# assembly sizes
asmbly.size <- unique(parsed[c(1, 6)])
# write csv with assembly sizes
write.csv(asmbly.size, 
          paste0("../results/", vert.invert, "/assembly_sizes.csv"), 
          row.names = FALSE)




