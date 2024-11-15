
## Takes a fasta file path as an input and extracts the contig names and 
## contig sizes of a species. gets the assembly size by summing the sizes of 
## all contigs except those identified by mitochondrial keywords. drops all but 
## the 2N longest contigs and outputs the results as a datatable
dataFromFasta <- function(fasta.path, chromnum.1n, mito.keywords, verbose) {
  # load fasta file as a vector
  fasta <- fread(fasta.path, header = FALSE, showProgress = verbose)$V1
  # get indices for the fasta vector
  header.ind <- which(grepl("^>", fasta))
  # get headers from fasta vector
  header <- fasta[header.ind]
  # get contig sizes from each header
  size.Mbp <- as.numeric(sub(".*:([0-9]+):1 REF", "\\1", header)) / 1000000
  # get contig names from each header
  name <- sub(">\\s*([^ ]+).*", "\\1", header)
  # get assembly size from all contigs excluding mitochondrial
  mito.indices <- which(!tolower(name) %in% mito.keywords)
  asmblysize.Gbp <- sum(na.omit(size.Mbp[mito.indices])) / 1000
  # create table
  fasta.data <- data.table(name, size.Mbp, asmblysize.Gbp)
  # sort by size
  fasta.data <- fasta.data[order(fasta.data$size.Mbp, decreasing = TRUE), ]
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
dataFromGtf <- function(gtf.file.path, name, mito.keywords, verbose) {
  # read gtf
  gtf <- read.table(gtf.file.path, header = TRUE, sep = "\t")
  # filter for genes only
  gtf <- gtf[which(gtf[, 3] == "gene"), ]
  # get the number of genes in each contig
  genecount <- sapply(
    name, function(name) table(gtf[[1]])[name])
  return(genecount)
}


## Takes a character vector of divsum files as an input. reads each file and
## calculates the mean K2P distance. outputs the mean K2P distance of each
## species as a dataframe
calcRepLandscStats <- function(species, file, asmblysz.Mbp) {
  # read text file into lines
  divsum.vector <- readLines(
    paste0("../results/vertebrates/repeat_landscape_divsums/", file))
  # look for the start of relevant information
  phrase <- "Coverage for each repeat class and divergence (Kimura)"
  start.index <- match(phrase, divsum.vector) + 1
  # condense relevant lines into a table
  divsum.vector <- divsum.vector[start.index:length(divsum.vector)]
  divsum.table <- read.table(textConnection(divsum.vector), 
                             sep = " ", 
                             header = TRUE)
  # drop NA columns
  divsum.table <- divsum.table[
    -c(which(sapply(divsum.table, function(col) all(is.na(col)))))]
  # vector of divergence scores
  divergence <- divsum.table$Div
  # vector of the frequencies of each divergence score
  frequency <- rowSums(divsum.table[, !names(divsum.table) == "Div"])
  # repeat content in Mbp
  repcontent.Mbp <- sum(frequency) / 1000000
  # the bin that contains the median
  median.bin <- which(cumsum(frequency) > sum(frequency)/2)[1]
  # frequency of the previous bin
  lower <- cumsum(frequency)[median.bin-1]
  # frequency of the next bin
  upper <- cumsum(frequency)[median.bin+1]
  # median frequency
  mid <- sum(frequency)/2
  # median bin
  k2p.median <- median.bin + (mid-lower)/(upper-lower)
  # calculate mean
  k2p.mean <- sum(divergence*frequency)/sum(frequency)
  # repeat content in percent coverage
  repcontent.pctcvrg <- (repcontent.Mbp / asmblysz.Mbp) * 100
  
  # build dataframe
  df <- data.frame(species, 
           repcontent.Mbp,
           repcontent.pctcvrg, 
           k2p.mean, 
           k2p.median)
  return(df)
}


calcPic <- function(col, species, tree) {
  df <- data.frame(species, col)
  df <- df[order(df$species, tree$tip.label), ]
  pic <- pic(df$col, phy = tree)
  return(pic)
}

permTest <- function(x, y, reps, method) {
  permuted.y <- replicate(reps, sample(y))
  permuted.cor <- apply(permuted.y, MARGIN = 2, function(col) cor(x, col, method = method))
  pval <- mean(abs(permuted.cor) >= abs(cor(x, y, method = method)))
  return(pval)
}
