



# constants
requiredPackages <- c("data.table", "ggplot2")

# paths
dataDirPath <- "../data"
parsedResults <- "../results/vertebrates/parsedResults.csv"
cladeDataCsvPath <- "../data/vertebrates/cladeData.csv"

# load stuff in
source("functions.R")
loadPackages(requiredPackages)
dat <- fread(parsedResults)


species <- unique(dat$species)

# if data file doesn't exist, create new file. else append to existing file
if (!file.exists(cladeDataCsvPath)) {
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


