
## Zhaobo Hu
## zhaobohu2002@gmail.com

## description: creates a csv file for the repeat landscape shell script. this
## file contains species name and family

library(ape)
results <- read.csv("../results/vertebrates/parsed.csv")
tree <- read.tree("../data/vertebrates/formatted_tree.nwk")

# drop species not in tree
results <- unique(results[1:34])
sp <- intersect(results$species, gsub("_", " ", tree$tip.label))
results <- results[results$species %in% sp, ]

# drop species without chromosome number
results <- results[!is.na(results$chromnum.1n), ]

# format species
results <- results[, c("species", "family")]
results$species <- na.omit(gsub(" ", "_", results$species))

# drop analyzed species
sp.analyzed <- list.files("../results/vertebrates/repeat_landscape_divsums")
sp.analyzed <- gsub(".divsum$", "", sp.analyzed)
results <- results[!(results$species %in% sp.analyzed), ]

write.table(results, 
            paste0("../data/vertebrates/landsc_ref.csv"), 
            sep = ",", 
            col.names = FALSE, 
            row.names = FALSE,
            quote = FALSE)
