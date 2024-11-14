
## Zhaobo Hu
## zhaobohu2002@gmail.com

## description: creates a csv file for the repeat landscape shell script. this
## file contains species name, family, and assembly size

# "vertebrates" or "invertebrates"?
vert.invert <- "vertebrates"

results <- read.csv(paste0("../data/", vert.invert, "/taxo_gnsz.csv"))
# get species with genome size estimates
results <- results[!is.na(results$est.gnsz.Mbp), ]
# subset for species and family
results <- results[, c("species", "family")]
# format species name
results$species <- na.omit(gsub(" ", "_", results$species))
# look for analyzed species
sp.analyzed <- list.files(paste0("../results/", vert.invert, "/repeat_landscape_divsums"))
# remove file extenison
sp.analyzed <- gsub(".divsum$", "", sp.analyzed)
# drop analyzed species
results <- results[!(results$species %in% sp.analyzed), ]
write.table(results, 
            paste0("../data/", vert.invert, "/landsc_ref.csv"), 
            sep = ",", 
            col.names = FALSE, 
            row.names = FALSE,
            quote = FALSE)
