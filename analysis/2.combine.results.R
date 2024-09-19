

## Zhaobo Hu
## zhaobohu2002@gmail.com

## description: combines the files in the results folder into a 
## single .csv file. Then, creates a .csv file with species names
## for manual curation of clade data. If this file already exists, 
## the script will add new species to the file if possible.

# constants
requiredPackages <- c("data.table", "ggplot2")

# paths
dataDirPath <- "../data/vertebrates"
resultsDirPath <- "../results/vertebrates"
individualResultsDirPath <- "../results/vertebrates/individualSpeciesResults"
cladeDataCsvPath <- "../data/vertebrates/cladeData.csv"

# load stuff in
source("functions.R")
loadPackages(requiredPackages)

# combine results
csvFullPaths <- getCsvFullPaths(individualResultsDirPath)
resultsList <- lapply(csvFullPaths, fread)
resultsDataTable <- na.omit(rbindlist(resultsList, fill = TRUE))
species <- unique(resultsDataTable$species)
fwrite(resultsDataTable, paste0(resultsDirPath, "/combinedResults.csv"), row.names = FALSE)

# if data file doesn't exist, create new file. else append to existing file
if (!file.exists(paste0(dataDirPath, "/cladeData.csv"))) {
  class <- NA
  order <- NA
  family <- NA
  cladeData <- data.table(species, class, order, family)
  fwrite(cladeData, cladeDataCsvPath, row.names = FALSE)
} else {
  cladeData <- fread(cladeDataCsvPath, header = TRUE)
  newSpecies <- species[which(!(species %in% cladeData$species))]
  if (length(newSpecies) > 0) {
    emptyMatrix <- matrix(nrow = length(newSpecies), ncol = ncol(cladeData) - 1)
    newData <- data.table(newSpecies, emptyMatrix)
    colnames(newData) <- colnames(cladeData)
    cladeData <- rbind(newData, cladeData)
    cladeData <- cladeData[order(cladeData$species), ]
    fwrite(cladeData, cladeDataCsvPath, row.names = FALSE, na = "NA")
  }
}







