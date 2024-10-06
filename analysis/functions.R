
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
  # create table
  fasta.data <- data.table(contig.name, contig.size_bp)
  # create new column with assembly sizes
  mito.indices <- which(!tolower(fasta.data$contig.name) %in% mito.keywords)
  asmbly.size_bp <- sum(fasta.data$contig.size_bp[mito.indices])
  fasta.data <- data.table(fasta.data, asmbly.size_bp)
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
  asmbly.gene.count <- nrow(gtf[!(gtf[[1]] %in% mito.keywords), ])
  # make table
  geneFreqTable <- table(gtf[[1]])
  # get the number of genes in each contig
  contig.gene.count <- sapply(
    contig.name, function(contig.name) table(gtf[[1]])[contig.name])
  gtf.data <- data.table(contig.gene.count, asmbly.gene.count)
  return(gtf.data)
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
    # calculate p-value
    fit <- summary(lm(table$contig.gene.count ~ table$contig.size_bp))
    species.pvalue <- fit$coefficients[2, 4]
    # add p-value and adjusted r-squared to the results
    species.rsquared <- fit$adj.r.squared
    table <- data.frame(table, species.rsquared, species.pvalue)
    return(table)
  }
}

## Takes a vector of species, looks for each species in "parsed_results.csv",
## and outputs a vector of assembly sizes for each species
getAsmblysz <- function(species) {
  return(as.numeric(unique(
    parsed.results$asmbly.size_bp[parsed.results$species == species])))
}

## Takes a character vector of divsum files as an input. reads each file and
## calculates the mean K2P distance. outputs the mean K2P distance of each
## species as a dataframe
getK2pMean <- function(files) {
  # read text file into lines
  divsum.vector <- readLines(
    paste0("../results/", vert.invert, "/repeat_landscape/", files))
  # look for the start of useful information
  phrase <- "Coverage for each repeat class and divergence (Kimura)"
  start.index <- match(phrase, divsum.vector) + 1
  # condense the useful lines into a table
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
  # calculate mean
  k2p.mean <- sum(divergence*frequency)/sum(frequency)
  return(k2p.mean)
}

## Takes a vector of species, looks for each species in "parsed_results.csv",
## and outputs a vector containing r-squared values for each species
getRsq <- function(species) {
  filter.names <- parsed.results$species == species
  if (any(filter.names)) {
    rsq <- unique(parsed.results$species.rsquared[filter.names])
    return(rsq)
  } else {
    return(NA)
  }
}

## Takes a vector of species, looks for each species in "clades_gnsz.csv",
## and outputs a vector containing class (taxonomy) assignments for each 
## species
getClass <- function(species) {
  filter.names <- clades.gnsz$species == species
  if (any(filter.names)) {
    class <- clades.gnsz$class[filter.names]
    return(class)
  } else {
    return(NA)
  }
}


## Takes a vector of species, looks for each species in "clades_gnsz.csv",
## and outputs a vector containing estimated genome sizes for each species
getGnszEst <- function(species) {
  filter.names <- clades.gnsz$species == species
  if (any(filter.names)) {
    gnsz <- clades.gnsz$genome.size.est_bp[filter.names]
    return(gnsz)
  } else {
    return(NA)
  }
}




getPIC <- function(dataframe, tree){
  sp.intersect <- intersect(tree$tip.label, gsub(" ", "_", dataframe$species))
  dataframe <- dataframe[dataframe$species %in% gsub("_", " ", sp.intersect), ]
  pruned.tree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% sp.intersect)])
  tree.species <- gsub("_", " ", pruned.tree$tip.label[-length(pruned.tree$tip.label)])
  pic <- pic(dataframe[, 2], pruned.tree)
  names(pic) <- tree.species
  pic <- pic[order(names(pic))]
  return(pic)
}



