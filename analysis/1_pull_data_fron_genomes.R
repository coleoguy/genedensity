
## Zhaobo Hu
## zhaobohu2002@gmail.com

## Description: reads a directory with fasta and gtf files to get the names of
## each species to be analyzed. Skips the species if it has already been
## analyzed (as per the existence of a result file). Then, finds the chromosome
## number of the species from a csv file (skips the species if chromosome 
## number is not found). Reads the fasta file of the species, and for each of
## the longest 2N contigs, outputs the contig name, contig size, and assembly
## size (assembly size is the same for all contigs of a species). Next, reads
## the gtf file of the species and searches for the number of genes for each
## contig name (also gets the number of genes in the assembly). Contigs with
## less than 2 genes are then dropped, and if the species has any contigs
## remaining, calculates the gene densities of each contig and writes the 
## results to a csv file

# "vertebrates" or "invertebrates"?
vert.invert <- "vertebrates"
# verbose
verbose <- FALSE

library(data.table)
source("functions.R")
# load chromosome numbers
paste0("../data/", vert.invert, "/chromnums.csv")
chromnums.csv <- read.csv(paste0("../data/", vert.invert, "/chromnums.csv"))
# keywords to match and exclude mitochondrial contigs for assembly size calcs
mito.keywords <- c("mt", "mito", "mitochondrial", "nonchromosomal")
# list genome files
genome.files <- list.files(paste0("../data/", vert.invert, "/genomes"))
# remove file extensions from genome file list
all.species.underscore <- unique(gsub("\\..*$", "", genome.files))
# remove underscores from genome file list to match for species chromnum csv file
all.species <- gsub("_", " ", all.species.underscore)
# begin loop
for (species in all.species) {
  if (verbose == TRUE) {
    start.time <- Sys.time()
    gc()
    print(noquote(paste0(species, " (", Sys.time(), ")")))
  }
  # proceed if result file is not found, else go to next species
  results.csv <- paste0("../results/", 
                        vert.invert, 
                        "/individual_species_results/", 
                        gsub(" ", "_", species), 
                        ".csv")
  if (!file.exists(results.csv)) {
    # get chromosome number
    chromnum.1n <- chromnums.csv$chromnum.1n[chromnums.csv$species == species]
    # proceed if chromosme number is available, else go to next species
    if (!is.na(chromnum.1n)) {
      fasta.file.path <- paste0("../data/", 
                                vert.invert, 
                                "/genomes/", 
                                gsub(" ", "_", species), 
                                ".fa")
      gtf.file.path <- paste0("../data/", 
                              vert.invert, 
                              "/genomes/", 
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
      # assemble datatable
      dat <- data.table(species,
                        contig.name, 
                        contig.size_bp, 
                        contig.gene.count,
                        asmbly.size_bp,
                        asmbly.gene.count)
      # drop contigs with less than 2 genes
      dat <- dat[dat$contig.gene.count >= 2, ]
      # proceed if data.table is not empty
      if (nrow(dat) != 0) {
        # calculate gene density
        contig.gene.dens_genes.per.bp <- dat$contig.gene.count/dat$contig.size_bp
        # add gene density to datatable
        dat2 <- data.table(
          dat[, 1:4],
          contig.gene.dens_genes.per.bp,
          dat[, 5:6]
        )
        # write results
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

