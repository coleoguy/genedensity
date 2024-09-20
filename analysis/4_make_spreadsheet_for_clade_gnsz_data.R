




# paths
clades.gnsz.csv <- "../data/vertebrates/clades_gnsz.csv"

parsed.results <- read.csv("../results/vertebrates/parsed_results.csv")
species <- unique(parsed.results$species)

# if data file doesn't exist, create new file. else append to existing file
if (!file.exists(clades.gnsz.csv)) {
  class <- NA
  order <- NA
  family <- NA
  genome.size.est_bp <- NA
  gnsz.source <- NA
  data.sheet <- data.frame(species, 
                          class, 
                          order, 
                          family, 
                          genome.size.est_bp, 
                          gnsz.source)
  write.csv(data.sheet, clades.gnsz.csv, row.names = FALSE)
} else {
  data.sheet <- read.csv(clades.gnsz.csv, header = TRUE)
  newSpecies <- species[which(!(species %in% data.sheet$species))]
  if (length(newSpecies) > 0) {
    emptyMatrix <- matrix(nrow = length(newSpecies), ncol = ncol(data.sheet) - 1)
    newData <- data.frame(newSpecies, emptyMatrix)
    colnames(newData) <- colnames(data.sheet)
    data.sheet <- rbind(newData, data.sheet)
    data.sheet <- data.sheet[order(data.sheet$species), ]
    write.csv(data.sheet, clades.gnsz.csv, row.names = FALSE, na = "NA")
  }
}

