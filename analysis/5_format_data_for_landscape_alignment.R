
## Zhaobo Hu
## zhaobohu2002@gmail.com

## description: creates a csv file for the repeat landscape shell script. this
## file contains species name, family, and assembly size

source("functions.R")
parsed.results <- read.csv("../results/vertebrates/parsed_results.csv")
clades.gnsz <- read.csv("../data/vertebrates/clades_gnsz.csv")
species <- clades.gnsz$species[!is.na(clades.gnsz$genome.size.est_bp)]
family <- clades.gnsz$family[!is.na(clades.gnsz$genome.size.est_bp)]
asmbly.size_bp <- sapply(species, getAsmblysz)
# reformat species names to match Ensembl's URL
species2 <- tolower(gsub(" ", "_", species))
results <- na.omit(data.frame(species2, family, asmbly.size_bp))
write.table(results, 
            "../data/vertebrates/species_family_asmblysz.csv", 
            sep = ",", 
            col.names = FALSE, 
            row.names = FALSE)
