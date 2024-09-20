

dataFromFasta <- function(fasta.file.path, chromnum.1n, mito.keywords) {
  # load fasta file as a vector
  fasta <- fread(fasta.file.path, header = FALSE, showProgress = TRUE)$V1
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

dataFromGtf <- function(gtf.file.path, contig.name, mito.keywords) {
  # read gtf
  gtf <- fread(gtf.file.path, header = FALSE, showProgress = TRUE)
  # filter for genes only
  gtf <- gtf[which(gtf[, 3] == "gene"), ]
  asmbly.gene.count <- nrow(gtf[!(gtf[[1]] %in% mito.keywords), ])
  # make table
  geneFreqTable <- table(gtf[[1]])
  # get the number of genes in each contig
  contig.gene.count <- sapply(contig.name, 
                              function(contig.name) table(gtf[[1]])[contig.name])
  gtf.data <- data.table(contig.gene.count, asmbly.gene.count)
  return(gtf.data)
}

getCsvFullPaths <- function(csv.dir.path) {
  csv.files <- list.files(csv.dir.path)
  csv.full.paths <- as.character(sapply(csv.files, function(csv.files) {
    paste0(csv.dir.path, "/", csv.files)
    }))
  return(csv.full.paths)
}

getAsmblysz <- function(species) {
  return(as.numeric(unique(parsed.results$asmbly.size_bp[parsed.results$species == species])))
}


getK2pMean <- function(files) {
  # read text file into lines
  divsum.vector <- readLines(paste0("../results/vertebrates/repeatLandscape/", files))
  # look for the start of useful information
  phrase <- "Coverage for each repeat class and divergence (Kimura)"
  start.index <- match(phrase, divsum.vector) + 1
  # condense the useful lines into a table
  divsum.vector2 <- divsum.vector[start.index:length(divsum.vector)]
  divsum.table <- read.table(textConnection(divsum.vector2), sep = " ", header = TRUE)
  # drop NA columns
  divsum.table2 <- divsum.table[-c(which(sapply(divsum.table, 
                                              function(col) all(is.na(col)))))]
  divergence <- divsum.table$Div
  frequency <- rowSums(divsum.table2[, !names(divsum.table2) == "Div"])
  k2p.mean <- sum(divergence*frequency)/sum(frequency)
  return(k2p.mean)
}


getRsq <- function(species) {
  filter.names <- parsed.results$species == species
  if (any(filter.names)) {
    rsq <- unique(parsed.results$species.rsquared[filter.names])
    return(rsq)
  } else {
    return(NA)
  }
}


getClass <- function(species) {
  filter.names <- clades.gnsz$species == species
  if (any(filter.names)) {
    class <- clades.gnsz$class[filter.names]
    return(class)
  } else {
    return(NA)
  }
}


