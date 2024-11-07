## Zhaobo Hu
## zhaobohu2002@gmail.com

## Description: makes a csv to record taxonomy data and estimated genome size
## for each species. if the csv already exists and new genomes have been
## analyzed, appends the new species to the existing csv without overwriting.

# "vertebrates" or "invertebrates"?
vert.invert <- "invertebrates"

taxo.gnsz.csv <- paste0("../data/", vert.invert, "/taxo_gnsz.csv")
input <- read.csv(paste0("../results/", vert.invert, "/all_contigs_results.csv"))
species <- unique(input$species)
# if data file doesn't exist, create new file. else append to existing file
if (!file.exists(taxo.gnsz.csv)) {
  class <- NA
  order <- NA
  family <- NA
  est.gnsz.Mbp <- NA
  est.gnsz.source <- NA
  est.gnsz.db <- NA
  data.sheet <- data.frame(species, 
                          class, 
                          order, 
                          family, 
                          est.gnsz.Mbp, 
                          est.gnsz.source,
                          est.gnsz.db)
  write.csv(data.sheet, taxo.gnsz.csv, row.names = FALSE)
} else {
  # existing data
  data.sheet <- read.csv(taxo.gnsz.csv, header = TRUE)
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
    write.csv(data.sheet, taxo.gnsz.csv, row.names = FALSE, na = "NA")
  }
}

