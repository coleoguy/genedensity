
## Zhaobo Hu
## zhaobohu2002@gmail.com

## Takes a species as input. drop contigs shorter than 10 million bp. drops
## assemblies that have less than 3 contigs or with p-values less than 0.05

# all contigs shorter than this size will be dropped
minContigSize <- 10000000

source("functions.R")
combined.results <- read.csv("../results/vertebrates/combined_results.csv")
species <- unique(combined.results$species)
results <- do.call(rbind, lapply(species, parseResults))
write.csv(results, file = "../results/vertebrates/parsed_results.csv", row.names = FALSE)
