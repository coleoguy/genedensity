


# Zhaobo Hu
# zhaobohu2002@gmail.com

# Description: reads files from a directory containing fasta and gtf/gff3 files
# of all species to be analyzed. The script will alphabetically sort all species
# in the directly and loop through each of them. Haploid chromosome numbers for 
# each species are recored in gene_density/data/chromnums.csv. If chromosome 
# number is recorded for a species, the script will read the fasta of the 
# species and extract the assembly size as well as the contig sizes of the 2N 
# longest contigs, where N is the haploid chromosome number. Next, gene counts 
# for each of the 2N longest contigs are gathered from the gtf/gff3 file. Gene
# density for each contig is calculated by dividing the gene count of each
# contig by the size of each contig. The results for each species is saved in
# its separate csv file


# verbose
verbose <- F
# is ensembl genome? setting to TRUE will run much faster
ensembl <- F
# keywords to match and exclude mitochondrial contigs for assembly size calcs
mito.keywords <- c("mt", "mito", "mitochondrial", "nonchromosomal")

# load library
library(data.table)
library(seqinr)
# source function
source("functions.R")
# load chromosome numbers
chromnums <- read.csv(paste0("../data/data.csv"))
# list genome files
genome.files <- list.files(paste0("../data/genomes"))
# pull species names from list of genome files
all.species <- gsub("_", " ", unique(gsub("\\..*$", "", genome.files)))
# list result files
results <- list.files(
  paste0("../results/individual_species_results"))
# pull species names from list of result files
result.species <- gsub("_", " ", unique(gsub("\\..*$", "", results)))
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
    chromnum.1n <- chromnums$chromnum.1n[chromnums$species == species]
    # proceed if chromosme number is available, else go to next species
    if (is.na(chromnum.1n)) {
      chromnum.1n <- chromnums$chromnum.est[chromnums$species == species]
    } 
    if (is.na(chromnum.1n)) {
      chromnum.1n <- 50
    } 
    # create full paths for files. sort alphabetically to place
    # the fasta file in front of the gff3/gtf file
    genome.files <- sort(paste0("../data/genomes/", 
                                list.files(paste0("../data/genomes"), 
                                           gsub(" ", "_", species))))
    # assume first file is fasta
    fasta.path <- genome.files[1]
    # assume second file is gff3/gtf/gbff
    annot.path <- genome.files[2]
    # read fasta
    fasta.data <- dataFromFasta(
      fasta.path, chromnum.1n, mito.keywords, verbose, ensembl)
    name <- fasta.data$name
    size.Mbp <- fasta.data$size.Mbp
    asmblysize.Gbp <- unique(fasta.data$asmblysize.Gbp)
    rm(fasta.data)
    gc()
    if (length(name) != 0) {
      # read gtf/gff3
      genecount <- dataFromGtf(annot.path, name, mito.keywords, verbose, annot)
      # proceed if gene count is available for more than one contig
      if (sum(is.na(genecount)) < length(name)) {
        # assemble datatable
        dat <- data.table(species,
                          name, 
                          size.Mbp, 
                          genecount)
        dat <- na.omit(dat)
        # calculate gene density
        genedens <- dat$genecount/dat$size.Mbp
        # add gene density to datatable
        dat <- data.table(
          dat,
          genedens,
          asmblysize.Gbp
        )
        # write results
        result.path <- paste0("../results/individual_species_results/", 
                              gsub(" ", "_", species), 
                              ".csv")
        fwrite(dat, file = result.path, row.names = FALSE)
        if (verbose == TRUE) {
          print(noquote("   Successfully written!"))
        }
      } else {
        # prompt to be displayed when sequence names in fasta and annotation 
        # files do not match
        if (verbose == TRUE) {
          print(noquote("   Aborted (sequence names in gtf and fasta do not match)"))
        }
      }
    } else {
      # prompt to be displayed when no contigs are found in the fasta file
      if (verbose == TRUE) {
        print(noquote("   Aborted (no contigs found in FASTA)"))
      }
    }
  } else {
    # prompt to be displayed when the current genome has already been analyzed
    if (verbose == TRUE) {
      print(noquote("   Aborted (file already exists)"))
    }
  }
  if (verbose == TRUE) {
    exec.time <- round(as.numeric(difftime(Sys.time(), start.time, units = "mins")), 2)
    print(noquote(paste0("   ", exec.time, " minutes")))
  }
}


