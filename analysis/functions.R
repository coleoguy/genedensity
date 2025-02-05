
## Takes a fasta file path as an input and extracts the contig names and 
## contig sizes of a species. gets the assembly size by summing the sizes of 
## all contigs except those identified by mitochondrial keywords. drops all but 
## the 2N longest contigs and outputs the results as a datatable
dataFromFasta <- function(fasta.path, chromnum.1n, mito.keywords, verbose, ensembl) {
  if (ensembl == TRUE) {
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
  } else {
    fasta <- read.fasta(fasta.path)
    size.bp <- c()
    for (i in 1:length(fasta)) {
      size.bp <- c(size.bp, length(fasta[[i]]))
    }
    size.Mb <- size.bp / 1000000
    name <- names(fasta)
  }
  # get assembly size from all contigs excluding mitochondrial
  mito.indices <- which(!tolower(name) %in% mito.keywords)
  asmblysize.Mb <- sum(na.omit(size.Mb[mito.indices]))
  # create table
  fasta.data <- data.table(name, size.Mb, asmblysize.Mb)
  # sort by size
  fasta.data <- fasta.data[order(fasta.data$size.Mb, decreasing = TRUE), ]
  # Keep large contigs
  fasta.data <- head(fasta.data, 2 * chromnum.1n)
  return(fasta.data)
}

## Takes a species' gtf file path as an input and reads the file. then, filter
## for rows that are annotated as "gene". for each contig name input, get the
## number of rows with matching contig names. also gets the number of nuclear
## genes in the assembly by counting the number of rows in the gtf minus the
## rows with contig names that match the mitochondrial keywords. outputs the
## results as a datatable
dataFromGtf <- function(annot.path, name, mito.keywords, verbose, annot) {
  # read gtf
  gtf <- read.table(annot.path, header = TRUE, sep = "\t")
  # filter for genes only
  gtf <- gtf[which(gtf[, 3] == "gene"), ]
  # get the number of genes in each contig
  genecount <- sapply(
    name, function(name) table(gtf[[1]])[name])
  return(genecount)
}


permTest <- function(x, y, reps, method) {
  permuted.y <- replicate(reps, sample(y))
  permuted.cor <- apply(permuted.y, MARGIN = 2, function(col) cor(x, col, method = method))
  pval <- mean(abs(permuted.cor) >= abs(cor(x, y, method = method)))
  return(pval)
}
