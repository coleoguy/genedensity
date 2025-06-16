
## Takes a fasta file path as an input and extracts the contig names and 
## contig sizes of a species. gets the assembly size by summing the sizes of 
## all contigs except those identified by mitochondrial keywords. drops all but 
## the 2N longest contigs and outputs the results as a datatable
dataFromFasta <- function(fasta.path, max.contig, verbose) {
  # load fasta file as a vector
  fasta <- fread(fasta.path, header = FALSE, showProgress = verbose)$V1
  # get indices for the fasta vector
  header.ind <- which(grepl("^>", fasta))
  # get headers from fasta vector
  header <- fasta[header.ind]
  # get contig sizes from each header
  size.Mb <- as.numeric(sub(".*:([0-9]+):1 REF", "\\1", header)) / 1000000
  # get contig names from each header
  name <- sub(">\\s*([^ ]+).*", "\\1", header)
  # get assembly size from all contigs
  asmblysize.Mb <- sum(size.Mb,  na.rm = T)
  # create table
  fasta.data <- data.table(name, size.Mb, asmblysize.Mb)
  # sort by size
  fasta.data <- fasta.data[order(fasta.data$size.Mb, decreasing = TRUE), ]
  # Keep large contigs
  fasta.data <- head(fasta.data, max.contig)
  return(fasta.data)
}

## Takes a species' gtf file path as an input and reads the file. then, filter
## for rows that are annotated as "gene". for each contig name input, get the
## number of rows with matching contig names. also gets the number of nuclear
## genes in the assembly by counting the number of rows in the gtf minus the
## rows with contig names that match the mitochondrial keywords. outputs the
## results as a datatable
dataFromGtf <- function(annot.path, name, verbose) {
  # read gtf
  gtf <- read.table(annot.path, header = FALSE, sep = "\t")
  # filter for genes only
  gtf <- gtf[which(gtf[, 3] == "gene"), ]
  # get the number of genes in each contig
  genecount <- sapply(
    name, function(name) table(gtf[[1]])[name])
  return(genecount)
}

# Shapiro-Wilk test
sw.test <- function(model) {
  res <- residuals(model)
  sw.p <- shapiro.test(res)$p.value
  return(sw.p)
}

# Pagel's lambda
lambda.test <- function(model) {
  res <- setNames(residuals(model), dat$species)
  lambda.p <- phylosig(tree, res, method = "lambda", test = TRUE, niter = 100)$P
  return(lambda.p)
}

# how many exons in a gene
library(data.table)
howmany <- function(genes, exons) {
  nrow(exons[exons$V4 >= as.numeric(genes[4]) 
             & exons$V5 <= as.numeric(genes[5]), ])
}

# assign alpha to color
co <- function(col, alpha = 1) {
  col.rgb <- c(col2rgb(col))
  new.col <- rgb(col.rgb[1], 
                 col.rgb[2], 
                 col.rgb[3], 
                 alpha = alpha * 255, 
                 maxColorValue = 255)
  return(new.col)
}