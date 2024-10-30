
## Takes a fasta file path as an input and extracts the contig names and 
## contig sizes of a species. gets the assembly size by summing the sizes of 
## all contigs except those identified by mitochondrial keywords. drops all but 
## the 2N longest contigs and outputs the results as a datatable
dataFromFasta <- function(fasta.file.path, chromnum.1n, mito.keywords, verbose) {
  # load fasta file as a vector
  fasta <- fread(fasta.file.path, header = FALSE, showProgress = verbose)$V1
  # get indices for the fasta vector
  contig.header.index <- which(grepl("^>", fasta))
  # get headers from fasta vector
  contig.header <- fasta[contig.header.index]
  # get contig sizes from each header
  contig.size_bp <- as.numeric(sub(".*:([0-9]+):1 REF", "\\1", contig.header))
  # get contig names from each header
  contig.name <- sub(">\\s*([^ ]+).*", "\\1", contig.header)
  # get assembly size from all contigs excluding mitochondrial
  mito.indices <- which(!tolower(contig.name) %in% mito.keywords)
  asmbly.size_bp <- sum(na.omit(contig.size_bp[mito.indices]))
  # create table
  fasta.data <- data.table(contig.name, contig.size_bp, asmbly.size_bp)
  # sort by size
  fasta.data <- fasta.data[order(fasta.data$contig.size_bp, decreasing = TRUE), ]
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
dataFromGtf <- function(gtf.file.path, contig.name, mito.keywords, verbose) {
  # read gtf
  gtf <- read.table(gtf.file.path, header = TRUE, sep = "\t")
  # filter for genes only
  gtf <- gtf[which(gtf[, 3] == "gene"), ]
  # get the number of genes in each contig
  contig.gene.count <- sapply(
    contig.name, function(contig.name) table(gtf[[1]])[contig.name])
  return(contig.gene.count)
}

## Takes a species as input. drop contigs shorter than 10 million bp. drops
## assemblies that have less than 3 contigs or with p-values less than 0.05
parseResults <- function(species, combined.results, min.contig.size) {
  # subset results for species
  table <- combined.results[combined.results$species == species, ]
  # drop contigs shorter than minContigSize
  table <- table[table$contig.size_bp >= min.contig.size, ]
  # if the are more than 2 contigs, continue. otherwise drop the assembly
  if (nrow(table) > 2) {
    return(table)
  }
}

## Takes a character vector of divsum files as an input. reads each file and
## calculates the mean K2P distance. outputs the mean K2P distance of each
## species as a dataframe
calcRepLandscapeStats <- function(files, asmbly.sz, vert.invert) {
  # read text file into lines
  divsum.vector <- readLines(
    paste0("../results/", vert.invert, "/repeat_landscape_divsums/", files))
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
  # repeat content in bp
  rep.content_bp <- sum(frequency)
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
  # species name, all lowercase with underscore
  sp.under <- gsub("_summary\\.divsum$", "", files)
  # species name, uppercase genus with underscore
  sp.under <- sub("^(\\w)", "\\U\\1", sp.under, perl = TRUE)
  # species name, uppercase genus without underscore
  species <- gsub("_", " ", sp.under)
  # assembly size
  asmbly.sz <- asmbly.sz$asmbly.size_bp[asmbly.sz$species == species]
  # repeat content in percent coverage
  rep.content_percent.of.assembly <- (rep.content_bp / asmbly.sz) * 100
  # build dataframe
  dat <- data.frame(species, 
                    rep.content_bp,
                    rep.content_percent.of.assembly, 
                    k2p.mean, 
                    k2p.median)
  return(dat)
}



calcContigStats <- function(species, results) {
  cursp <- results[results$species == species, ]
  fit <- summary(glm(cursp$contig.gene.count ~ cursp$contig.size_bp))
  beta <- fit$coefficients[2, 1]
  pval.beta <- fit$coefficients[2, 4]
  rsq <- summary(lm(cursp$contig.gene.count ~ cursp$contig.size_bp))$r.squared
  coef.of.var <-sd(cursp$contig.gene.dens_genes.per.bp) / mean(cursp$contig.gene.dens_genes.per.bp)
  contig.stats <- data.frame(beta, pval.beta, rsq, coef.of.var)
  return(contig.stats)
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
