# Zhaobo Hu
# zhaobohu2002@gmail.com

# Runs on Ensembl genome assemblies. If an assembly has enough large contigs, 
# fit a linear model using these contigs where contig size predicts the number of
# genes on each contig. Record the R2 value as a measure of gene density homogeneity.
# Also computes Gini coefficient and coefficient of variation (CV) of gene density
# as direct dispersion measures.

# load library
library(data.table)
library(DescTools)
library(future.apply)
source("functions.R")

# verbose
verbose <- T

# list genome files
genome.files <- list.files(paste0("../data/genomes"))
genome.files <- genome.files[genome.files != "readme.txt"]

# pull species names from list of genome files
all.species <- unique(gsub("\\..*$", "", genome.files))
max.contig <- 60

# set up parallel workers (set to number of available cores)
plan(multisession, workers = 12)

# begin parallel loop
results.list <- future_lapply(seq_along(all.species), function(i) {

  # assume first file is fasta
  # fasta.path <- paste0("../data/genomes/", all.species[i], ".fa")
  fai.path <- paste0("../data/genomes/", all.species[i], ".fa.fai")
  # assume second file is gff3/gtf/gbff
  annot.path <- paste0("../data/genomes/", all.species[i], ".gtf")
  # read fasta
  #fasta.data <- dataFromFasta(fasta.path = fasta.path,
  #                            max.contig = max.contig,
  #                            verbose = FALSE)
  fasta.data <- dataFromFai(fai.path = fai.path, max.contig = max.contig)
  fasta.data <- fasta.data[fasta.data$size.Mb >= 10, ]

  # skip to next species if we have less than three retained contigs
  if (nrow(fasta.data) < 3) return(NULL)

  # skip to next species if sum of captured size < 0.8 of assembly size
  if (unique(sum(fasta.data$size) < 0.8 * fasta.data$asmblysize.Mb)) return(NULL)

  # get names and sizes of retained contigs
  name <- fasta.data$name
  size.Mb <- fasta.data$size.Mb
  asmblysize.Mb <- fasta.data$asmblysize.Mb[1]
  rm(fasta.data)

  # read gtf
  genecount <- dataFromGtf(annot.path, name, verbose = FALSE)

  # skip to next species if gene count is unavailable for at least 3 contigs
  if (sum(!is.na(genecount)) < 3) return(NULL)

  # assemble datatable
  dat <- data.table(species = all.species[i],
                    size.Mb,
                    genecount)
  dat <- na.omit(dat)
  # calculate and add gene density
  dat$genedens <- dat$genecount / dat$size.Mb

  # R-squared
  rsq <- summary(lm(dat$genecount ~ dat$size.Mb))$r.squared

  # Gini coefficient (length-weighted by contig size)
  gini <- Gini(dat$genedens, weights = dat$size.Mb)

  # coefficient of variation
  cv <- sd(dat$genedens) / mean(dat$genedens)

  return(data.frame(species  = all.species[i],
                    rsq      = rsq,
                    gini     = gini,
                    cv       = cv,
                    assem.sz = round(asmblysize.Mb)))
})

# print progress (species name for any that returned NULL)
for (i in seq_along(results.list)) {
  if (is.null(results.list[[i]])) {
    print(paste("Aborted:", all.species[i]))
  } else {
    print(paste("Done:", all.species[i]))
  }
}

# collapse results
results <- do.call(rbind, results.list)

tax <- read.csv("../data/taxonomy.csv")
tax$clade <- tax$class
tax[tax$clade %in% "Aves", ]$clade <- "Sauropsida"
tax[tax$clade %in% "Reptilia", ]$clade <- "Sauropsida"
tax[!(tax$clade %in% c("Actinopterygii", "Mammalia", "Sauropsida")), ]$clade <- "Others"
tax$species <- tolower(tax$species)
results <- merge(results, tax, by = "species", all.x = T)

write.csv(results, row.names = F, file = "../results/rsq.csv")
