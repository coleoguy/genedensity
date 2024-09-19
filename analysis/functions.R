

dataFromFasta <- function(fastaFilePath, chromNum.1n, mitoContigKeywords) {
  # load fasta file as a vector
  fasta <- fread(fastaFilePath, header = FALSE, showProgress = TRUE)$V1
  # get indices for the fasta vector
  contigHeaderIndex <- which(grepl("^>", fasta))
  # get headers from fasta vector
  contigHeader <- fasta[contigHeaderIndex]
  # get contig sizes from each header
  contigSize.bp <- as.numeric(sub(".*:([0-9]+):1 REF", "\\1", contigHeader))
  # get contig names from each header
  contigName <- sub(">\\s*([^ ]+).*", "\\1", contigHeader)
  # create table
  fastaData <- data.table(contigName, contigSize.bp)
  # create new column with assembly sizes
  asmblySize.bp <- sum(fastaData$contigSize.bp[which(!tolower(fastaData$contigName) %in% mitoContigKeywords)])
  fastaData <- data.table(fastaData, asmblySize.bp)
  # sort by size
  fastaData <- fastaData[order(fastaData$contigSize.bp, decreasing = TRUE), ]
  # Keep large contigs
  fastaData <- head(fastaData, 2 * chromNum.1n)
  return(fastaData)
}

dataFromGtf <- function(gtfFilePath, contigName, mitoContigKeywords) {
  # read gtf
  gtf <- fread(gtfFilePath, header = FALSE, showProgress = TRUE)
  # filter for genes only
  gtf <- gtf[which(gtf[, 3] == "gene"), ]
  asmblyGeneCount <- nrow(gtf[!(gtf[[1]] %in% mitoContigKeywords), ])
  # make table
  geneFreqTable <- table(gtf[[1]])
  # get the number of genes in each contig
  contigGeneCount <- sapply(contigName, function(contigName) table(gtf[[1]])[contigName])
  gtfData <- data.table(contigGeneCount, asmblyGeneCount)
  return(gtfData)
}

getCsvFullPaths <- function(csvDirPath) {
  csvFiles <- list.files(csvDirPath)
  csvFullPaths <- as.character(sapply(csvFiles, function(csvFiles) paste0(csvDirPath, "/", csvFiles)))
  return(csvFullPaths)
}

loadPackages <- function(requiredPackages) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  installedPackages <- requiredPackages %in% rownames(installed.packages())
  if (any(installedPackages == FALSE)) {
    BiocManager::install(requiredPackages[!installedPackages])
  }
  invisible(lapply(requiredPackages, library, character.only = TRUE))
}



