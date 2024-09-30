
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
chromnums.table <- read.csv(paste0("../data/", vert.invert, "/chromnums.csv"))
# keywords to match and exclude mitochondrial contigs for assembly size calcs
mito.keywords <- c("mt", "mito", "mitochondrial", "nonchromosomal")
# list genome files
genome.files <- list.files(paste0("../data/", vert.invert, "/genomes"))
# pull species names from list of genome files
all.species <- gsub("_", " ", unique(gsub("\\..*$", "", genome.files)))
# list result files
result.files <- list.files(
  paste0("../results/", vert.invert, "/individual_species_results"))
# pull species names from list of result files
result.species <- gsub("_", " ", unique(gsub("\\..*$", "", result.files)))
# begin loop
for (species in all.species) {
  if (verbose == TRUE) {
    start.time <- Sys.time()
    gc()
    print(noquote(paste0(species, " (", Sys.time(), ")")))
  }
  # proceed if result file is not found, else go to next species
  if (!(species %in% result.species)) {
    # get chromosome number
    chromnum.1n <- chromnums.table$chromnum.1n[chromnums.table$species == species]
    # proceed if chromosme number is available, else go to next species
    if (!is.na(chromnum.1n)) {
      list.files(paste0("../data/", vert.invert, "/genomes"), gsub(" ", "_", species))
      fa.gtf.file.path <- sort(paste0("../data/", 
                                vert.invert, 
                                "/genomes/", 
                                list.files(paste0("../data/", 
                                                  vert.invert, 
                                                  "/genomes"), 
                                           gsub(" ", "_", species))))
      fasta.file.path <- fa.gtf.file.path[1]
      gtf.file.path <- fa.gtf.file.path[2]
      # read fasta
      fasta.data <- dataFromFasta(
        fasta.file.path, chromnum.1n, mito.keywords, verbose)
      contig.name <- fasta.data$contig.name
      contig.size_bp <- fasta.data$contig.size_bp
      asmbly.size_bp <- fasta.data$asmbly.size_bp
      rm(fasta.data)
      gc()
      # read gtf
      gtf.data <- dataFromGtf(gtf.file.path, contig.name, mito.keywords, verbose)
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
        dat <- data.table(
          dat[, 1:4],
          contig.gene.dens_genes.per.bp,
          dat[, 5:6]
        )
        # write results
        result.file <- paste0("../results/", 
                              vert.invert, 
                              "/individual_species_results/", 
                              gsub(" ", "_", species), 
                              ".csv")
        fwrite(dat, file = result.file, row.names = FALSE)
        if (verbose == TRUE) {
          print(noquote("   Successfully written!"))
        }
      } else {
        if (verbose == TRUE) {
          print(noquote("   Aborted (sequence names in gtf and fasta do not match)"))
        }
      }
    } else {
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

