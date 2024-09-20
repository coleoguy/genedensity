
# load package
library(data.table)

# constants
minContigSize <- 10000000

combined.results <- fread("../results/vertebrates/combined_results.csv")
result <- data.table()
species <- unique(combined.results$species)
for (i in species) {
  # subset results for species
  table <- combined.results[combined.results$species == i]
  # filter results by contig size
  table2 <- table[table$contig.size_bp >= minContigSize]
  # if has more tahn 2 contigs
  if (nrow(table2) > 2) {
    # calculate p-value
    fit <- summary(lm(table2$contig.gene.count ~ table2$contig.size_bp))
    species.pvalue <- fit$coefficients[2, 4]
    # if p-value is less than 0.05
    if (species.pvalue < 0.05) {
      # add p-value and adjusted r-squared
      species.rsquared <- fit$adj.r.squared
      table3 <- data.frame(table2, species.rsquared, species.pvalue)
      result <- rbind(result, table3)
    }
  }
}

fwrite(result, file = "../results/vertebrates/parsed_results.csv", row.names = FALSE)
