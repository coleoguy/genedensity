

## Zhaobo Hu
## zhaobohu2002@gmail.com

## description: creates a csv file for the repeat landscape shell script. this
## file contains species name, family, and assembly size

# constants
requiredPackages <- c("data.table", "ggplot2")

# paths
resultsDirPath <- "../results/vertebrates"
parsedResultsCsvPath <- "../results/vertebrates/parsedResults.csv"
cladeDataCsvPath <- "../data/vertebrates/cladeData.csv"

# load stuff in
source("functions.R")
loadPackages(requiredPackages)
parsedResults <- fread(parsedResultsCsvPath)
cladeData <- fread(cladeDataCsvPath)

# make file
species <- cladeData$species[!is.na(cladeData$estimatedGenomeSize.bp)]
family <- cladeData$family[!is.na(cladeData$estimatedGenomeSize.bp)]
asmblySize.bp <- sapply(species, function(sp) 
  return(as.numeric(unique(
    parsedResults$asmblySize.bp[
      parsedResults$species == sp]))))
species2 <- tolower(gsub(" ", "_", species))
results <- na.omit(data.table(species2, family, asmblySize.bp))
fwrite(results, 
       paste0(resultsDirPath, "/speciesFamilyAsmblysz.csv"), 
       row.names = FALSE, 
       col.names = FALSE)


