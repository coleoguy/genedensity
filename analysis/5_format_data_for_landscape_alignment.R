
## Zhaobo Hu
## zhaobohu2002@gmail.com

## description: creates a csv file for the repeat landscape shell script. this
## file contains species name, family, and assembly size

# "vertebrates" or "invertebrates"?
vert.invert <- "vertebrates"

source("functions.R")
parsed.results <- read.csv(
  paste0("../results/", vert.invert, "/parsed_results.csv"))
clades.gnsz <- read.csv("../data/", vert.invert, "/clades_gnsz.csv")
species <- clades.gnsz$species[!is.na(clades.gnsz$genome.size.est_bp)]
family <- clades.gnsz$family[!is.na(clades.gnsz$genome.size.est_bp)]
asmbly.size_bp <- sapply(species, getAsmblysz)
# reformat species names to match Ensembl's URL
species2 <- tolower(gsub(" ", "_", species))
results <- na.omit(data.frame(species2, family, asmbly.size_bp))
write.table(results, 
            paste0("../data/", vert.invert, "/species_family_asmblysz.csv"), 
            sep = ",", 
            col.names = FALSE, 
            row.names = FALSE,
            quote = FALSE)
