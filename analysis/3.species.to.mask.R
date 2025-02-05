
## Zhaobo Hu
## zhaobohu2002@gmail.com

## description: creates a csv file for the repeat landscape shell script. this
## file contains species name and family

results <- read.csv("../results/parsed.csv")
results <- results[results$thrs == 0, ]
results <- results[!duplicated(results$species), ]

# drop species not in tree
# library(ape)
# tree <- read.tree("../data/formatted_tree.nwk")
# sp <- intersect(results$species, gsub("_", " ", tree$tip.label))
# results <- results[results$species %in% sp, ]

# drop species without chromosome number
results <- results[!is.na(results$chromnum.1n), ]

# format species
results <- results[, c("species", "family")]
results$species <- na.omit(gsub(" ", "_", results$species))

# drop analyzed species
sp.analyzed <- list.files("../results/divsums")
sp.analyzed <- gsub(".divsum$", "", sp.analyzed)
results <- results[!(results$species %in% sp.analyzed), ]

write.table(results, 
            paste0("../data/to.mask.csv"), 
            sep = ",", 
            col.names = FALSE, 
            row.names = FALSE,
            quote = FALSE)
