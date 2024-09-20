
## Zhaobo Hu
## zhaobohu2002@gmail.com

## Description: reads a reference CSV to get species names and
## chromosome numbers for each species. Then, for each species, pulls
## information from fasta/fa and gtf/gff3 files to calculate gene
## density and other relevant information for the longest contigs in
## the assembly. The number of contigs considered for each species is
## equal to the 2n chromosome number of the species. Then, calculates
## the r-squared values for genes~contigSize as well as the length of
## masked areas in each contig.

# load stuff in
library(data.table)
source("functions.R")
chromnums.csv <- read.csv("../data/vertebrates/chromnums.csv")

# constants
mito.keywords <- c("mt", "mito", "mitochondrial", "nonchromosomal")
verbose <- TRUE

genome.files <- list.files("../data/vertebrates/genomes")
all.species.underscore <- unique(gsub("\\..*$", "", genome.files))
all.species <- gsub("_", " ", all.species.underscore)

# begin loop
for (species in all.species) {
  if (verbose == TRUE) {
    start.time <- Sys.time()
    gc()
    print(noquote(paste0(species, " (", Sys.time(), ")")))
  }
  # proceed if result file is not found
  results.csv <- paste0("../results/vertebrates/individual_species_results/", 
                        gsub(" ", "_", species), 
                        ".csv")
  if (!file.exists(results.csv)) {
    # get chromosome number
    chromnum.1n <- chromnums.csv$chromnum.1n[chromnums.csv$species == species]
    # proceed if chromosme number is available
    if (!is.na(chromnum.1n)) {
      fasta.file.path <- paste0("../data/vertebrates/genomes/", 
                              gsub(" ", "_", species), 
                              ".fa")
      gtf.file.path <- paste0("../data/vertebrates/genomes/", 
                              gsub(" ", "_", species), 
                              ".gtf")
      # read fasta
      fasta.data <- dataFromFasta(fasta.file.path, chromnum.1n, mito.keywords)
      contig.name <- fasta.data$contig.name
      contig.size_bp <- fasta.data$contig.size_bp
      asmbly.size_bp <- fasta.data$asmbly.size_bp
      rm(fasta.data)
      gc()
      # read gtf
      gtf.data <- dataFromGtf(gtf.file.path, contig.name, mito.keywords)
      contig.gene.count <- gtf.data$contig.gene.count
      asmbly.gene.count <- gtf.data$asmbly.gene.count
      rm(gtf.data)
      gc()
      # assemble dataframe
      dat <- data.table(species,
                        contig.name, 
                        contig.size_bp, 
                        contig.gene.count,
                        asmbly.size_bp,
                        asmbly.gene.count)
      # drop contigs with less than 2 genes
      dat <- dat[dat$contig.gene.count >= 2, ]
      # proceed if there are rows in the dataframe
      if (nrow(dat) != 0) {
        fit <- summary(lm(dat$contig.gene.count ~ dat$contig.size_bp))
        contig.gene.dens_genes.per.bp <- dat$contig.gene.count/dat$contig.size_bp
        dat2 <- data.table(
          dat[, 1:4],
          contig.gene.dens_genes.per.bp,
          dat[, 5:6]
        )
        fwrite(dat2, file = results.csv, row.names = FALSE)
        if (verbose == TRUE) {
          print(noquote("   Successfully written!"))
        }
      } else {
        fwrite(data.table(), file = results.csv, row.names = FALSE)
        if (verbose == TRUE) {
          print(noquote("   Aborted (no genes that match contig names)"))
        }
      }
    } else {
      fwrite(data.table(), file = results.csv, row.names = FALSE)
      if (verbose == TRUE) {
        print(noquote("   Aborted (no chromosome number estimate)"))
      }
    }
  } else {
    if (verbose == TRUE) {
      print(noquote("   Aborted (file already exists)"))
    }
  }
  if (verbose == TRUE) {
    exec.time <- round(as.numeric(difftime(Sys.time(), start.time, units = "mins")), 2)
    print(noquote(paste0("   ", exec.time, " minutes")))
  }
}

