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

# load library
library(data.table)
source("functions.R")

# verbose
verbose <- T

# list genome files
genome.files <- list.files(paste0("../data/genomes"))
genome.files <- genome.files[genome.files != "readme.txt"]

# pull species names from list of genome files
all.species <- unique(gsub("\\..*$", "", genome.files))
max.contig <- 60

# make a results object
results <- as.data.frame(matrix(NA, 0, 3))
colnames(results) <- c("species","rsq","assem.sz")
# begin loop
for (i in c(1:length(all.species))){ # TODO change range
  print(paste("Working on", all.species[i]))
  # assume first file is fasta
  fasta.path <- paste0("../data/genomes/",all.species[i], ".fa")
  # assume second file is gff3/gtf/gbff
  annot.path <- paste0("../data/genomes/",all.species[i], ".gtf")
  # read fasta
  fasta.data <- dataFromFasta(fasta.path = fasta.path, 
                              max.contig = max.contig,
                              verbose = TRUE)
  fasta.data <- fasta.data[fasta.data$size.Mb >= 10, ]
  
  # skip to next species if we have less than three retained contigs
  if(nrow(fasta.data) < 3){
    if (verbose == TRUE) {
      print(noquote("   Aborted (less than 3 retained contigs)"))
    }
    next
  }
  
  # skip to next species if sum of captured size < 0.8 of assembly size
  if(unique(sum(fasta.data$size) < 0.8 * fasta.data$asmblysize.Mb)) {
    if (verbose == TRUE) {
      print(noquote("   Aborted (low assembly contiguity)"))
    }
    next
  }
  
  # get names and sizes of retained contigs
  name <- fasta.data$name
  size.Mb <- fasta.data$size.Mb
  asmblysize.Mb <- fasta.data$asmblysize.Mb[1]
  rm(fasta.data)
  gc()
  
  # read gtf
  genecount <- dataFromGtf(annot.path, name, verbose)
  
  # skip to next species if gene count is unavailable for at least 3 contigs
  if (sum(!is.na(genecount)) < 3) {
    if (verbose == TRUE) {
      print(noquote("   Aborted (insufficient contigs with data)"))
    }
    next
  }
  
  # assemble datatable
  dat <- data.table(all.species[i],
                    size.Mb, 
                    genecount)
  dat <- na.omit(dat)
  # calculate and add gene density
  dat$genedens <- dat$genecount/dat$size.Mb
  rsq <- summary(lm(dat$genecount~dat$size.Mb))$r.squared
  j <- nrow(results) + 1
  results[j, 1] <- all.species[i]
  results[j, 2:3] <- c(rsq, round(asmblysize.Mb))
}

tax <- read.csv("../data/taxonomy.csv")
tax$clade <- tax$class
tax[tax$clade %in% "Aves", ]$clade <- "Sauria"
tax[tax$clade %in% "Reptilia", ]$clade <- "Sauria"
tax[!(tax$clade %in% c("Actinopterygii", "Mammalia", "Sauria")), ]$clade <- "Others"
results <- merge(results, tax, by = "species", all.x = T)

write.csv(results, row.names=F, file = "../results/rsq.08.csv") # TODO change name
