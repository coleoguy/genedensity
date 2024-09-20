

## Zhaobo Hu
## zhaobohu2002@gmail.com

## description: combines the files in the results folder into a 
## single .csv file. Then, creates a .csv file with species names
## for manual curation of clade data. If this file already exists, 
## the script will add new species to the file if possible.

# load package
library(data.table)

# source functions
source("functions.R")

# combine results
csvFullPaths <- getCsvFullPaths("../results/vertebrates/individual_species_results")
resultsList <- lapply(csvFullPaths, fread)
resultsDataTable <- na.omit(rbindlist(resultsList, fill = TRUE))
species <- unique(resultsDataTable$species)
fwrite(resultsDataTable, "../results/vertebrates/combined_results.csv", row.names = FALSE)










