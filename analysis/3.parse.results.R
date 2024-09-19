


# constants
requiredPackages <- c("data.table", "ggplot2")
minContigSize <- 10000000

# paths
functionsPath <- "functions.R"
combinedResults <- "../results/vertebrates/combinedResults.csv"

# load stuff in
source(functionsPath)
loadPackages(requiredPackages)
dat <- fread(combinedResults)


result <- data.table()
species <- unique(dat$species)
for (i in species) {
  # subset results for species
  table <- dat[dat$species == i]
  # filter results by contig size
  table2 <- table[table$contigSize.bp >= minContigSize]
  # if has more tahn 2 contigs
  if (nrow(table2) > 2) {
    # calculate p-value
    fit <- summary(lm(table2$contigGeneCount ~ table2$contigSize.bp))
    speciesPvalue <- fit$coefficients[2, 4]
    # if p-value is less than 0.05
    if (speciesPvalue < 0.05) {
      # add p-value and adjusted r-squared
      speciesRsquared <- fit$adj.r.squared
      table3 <- data.frame(table2, speciesRsquared, speciesPvalue)
      result <- rbind(result, table3)
    }
  }
}

fwrite(result, file = "../results/vertebrates/parsedResults.csv", row.names = FALSE)



