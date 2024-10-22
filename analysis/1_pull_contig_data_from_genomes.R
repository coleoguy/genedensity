


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

# "vertebrates" or "invertebrates"?
vert.invert <- "invertebrates"
# verbose
verbose <- F
# keywords to match and exclude mitochondrial contigs for assembly size calcs
mito.keywords <- c("mt", "mito", "mitochondrial", "nonchromosomal")

# load library
library(data.table)
# source function
source("functions.R")
# load chromosome numbers
chromnums.table <- read.csv(paste0("../data/", vert.invert, "/chromnums.csv"))
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
      # list files corresponding to the current species name (there shoule be two)
      list.files(paste0("../data/", vert.invert, "/genomes"), gsub(" ", "_", species))
      # create full paths for the previous files. sort alphabetically to place
      # the fasta file in front of the gff3/gtf file
      fa.gtf.file.path <- sort(paste0("../data/", 
                                vert.invert, 
                                "/genomes/", 
                                list.files(paste0("../data/", 
                                                  vert.invert, 
                                                  "/genomes"), 
                                           gsub(" ", "_", species))))
      # assume the first file to be the fasta
      fasta.file.path <- fa.gtf.file.path[1]
      # assume the second file to be the gtf/gff3
      gtf.file.path <- fa.gtf.file.path[2]
      # read fasta
      fasta.data <- dataFromFasta(
        fasta.file.path, chromnum.1n, mito.keywords, verbose)
      contig.name <- fasta.data$contig.name
      contig.size_bp <- fasta.data$contig.size_bp
      asmbly.size_bp <- fasta.data$asmbly.size_bp
      rm(fasta.data)
      gc()
      if (length(contig.name) != 0) {
        # read gtf/gff3
        contig.gene.count <- dataFromGtf(gtf.file.path, contig.name, mito.keywords, verbose)
        # proceed if gene count is available for more than one contig
        if (sum(is.na(contig.gene.count)) < length(contig.name)) {
          # assemble datatable
          dat <- data.table(species,
                            contig.name, 
                            contig.size_bp, 
                            contig.gene.count,
                            asmbly.size_bp)
          dat <- na.omit(dat)
          # calculate gene density
          contig.gene.dens_genes.per.bp <- dat$contig.gene.count/dat$contig.size_bp
          # add gene density to datatable
          dat <- data.table(
            dat[, 1:4],
            contig.gene.dens_genes.per.bp,
            dat[, 5]
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
          # prompt to be displayed when sequence names in fasta and annotation 
          # files do not match
          if (verbose == TRUE) {
            print(noquote("   Aborted (sequence names in gtf/gff3 and fasta do not match)"))
          }
        }
      } else {
        # prompt to be displayed when no contigs are found in the fasta file
        if (verbose == TRUE) {
          print(noquote("   Aborted (no contigs found in FASTA)"))
        }
      }
    } else {
      # prompt to be displayed when chromosome numbers are unavailable
      if (verbose == TRUE) {
        print(noquote("   Aborted (no chromosome number)"))
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
write.csv(chromnums.table, 
          paste0("../results/", vert.invert, "/final_results.csv"),
          row.names = FALSE)


