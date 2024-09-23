## Zhaobo Hu
## zhaobohu2002@gmail.com

## Description: makes a csv to record taxonomy data and estimated genome size
## for each species. if the csv already exists and new genomes have been
## analyzed, appends the new species to the existing csv.

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
  # existing data
  data.sheet <- read.csv(clades.gnsz.csv, header = TRUE)
  # list of new species
  new.species <- species[which(!(species %in% data.sheet$species))]
  if (length(new.species) > 0) {
    # creates an empty matrix with row number equal to number of new species and
    # column number 1 less than the width of available data
    empty.matrix <- matrix(nrow = length(new.species), ncol = ncol(data.sheet) - 1)
    # add a column of new species to the left of the empty matrix
    new.data <- data.frame(new.species, empty.matrix)
    # name the columns the same as the existing data to allow row binding
    colnames(new.data) <- colnames(data.sheet)
    # bind rows
    data.sheet <- rbind(new.data, data.sheet)
    # sort rows by alphabetical order of species
    data.sheet <- data.sheet[order(data.sheet$species), ]
    write.csv(data.sheet, clades.gnsz.csv, row.names = FALSE, na = "NA")
  }
}

