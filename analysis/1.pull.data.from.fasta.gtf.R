
## Zhaobo Hu
## zhaobohu2002@gmail.com

## Description: reads a reference CSV to get species names and
## chromosome numbers for each species. Then, for each species, pulls
## information from fasta/fa and gtf/gff3 files to calculate gene
## density and other relevant information for the longest contigs in
## the assembly. The number of contigs considered for each species is
## equal to the 2n chromosome number of the species. Then, calculates
## the r-squared values for genes~contigSize as well as the length of
## masked areas in each contig.

# constants
requiredPackages <- c("data.table", "ggplot2")
mitoContigKeywords <- c("mt", "mito", "mitochondrial", "nonchromosomal")

# paths
speciesChromnumCsvPath <- "../data/vertebrates/speciesChromnum.csv"
fastaDirPath <- "../data/vertebrates/genomes"
gtfDirPath <- "../data/vertebrates/genomes"
individualResultsDirPath <- "../results/vertebrates/individualSpeciesResults"

# load stuff in
source("functions.R")
loadPackages(requiredPackages)
speciesChromnumCsv <- fread(speciesChromnumCsvPath)

# begin loop
for (species in speciesChromnumCsv$species) {
  startTime <- Sys.time()
  gc()
  print(noquote(paste0(species, " (", Sys.time(), ")")))
  # proceed if result file is not found
  resultsCsv <- paste0(individualResultsDirPath, "/", gsub(" ", "_", species), ".csv")
  if (!file.exists(resultsCsv)) {
    # get chromosome number
    chromNum.1n <- speciesChromnumCsv$chromNum.1n[speciesChromnumCsv$species == species]
    # proceed if chromosme number is available
    if (!is.na(chromNum.1n)) {
      fastaFilePath <- paste0(fastaDirPath, "/", gsub(" ", "_", species), ".fa")
      gtfFilePath <- paste0(gtfDirPath, "/", gsub(" ", "_", species), ".gtf")
      # read fasta
      fastaData <- dataFromFasta(fastaFilePath, chromNum.1n, mitoContigKeywords)
      contigName <- fastaData$contigName
      contigSize.bp <- fastaData$contigSize.bp
      asmblySize.bp <- fastaData$asmblySize.bp
      rm(fastaData)
      gc()
      # read gtf
      gtfData <- dataFromGtf(gtfFilePath, contigName, mitoContigKeywords)
      contigGeneCount <- gtfData$contigGeneCount
      asmblyGeneCount <- gtfData$asmblyGeneCount
      rm(gtfData)
      gc()
      # assemble dataframe
      dat <- data.table(species,
                        contigName, 
                        contigSize.bp, 
                        contigGeneCount,
                        asmblySize.bp,
                        asmblyGeneCount)
      # drop contigs with less than 2 genes
      dat <- dat[dat$contigGeneCount >= 2, ]
      # proceed if there are rows in the dataframe
      if (nrow(dat) != 0) {
        fit <- summary(lm(dat$contigGeneCount~dat$contigSize.bp))
        contigGeneDens.genesPerBp <- dat$contigGeneCount/dat$contigSize.bp
        dat2 <- data.table(
          dat[, 1:4],
          contigGeneDens.genesPerBp,
          dat[, 5:6]
        )
        fwrite(dat2, file = resultsCsv, row.names = FALSE)
        print(noquote("   Successfully written!"))
      } else {
        fwrite(data.table(), file = resultsCsv, row.names = FALSE)
        print(noquote("   Aborted (no genes that match contig names)"))
      }
    } else {
      fwrite(data.table(), file = resultsCsv, row.names = FALSE)
      print(noquote("   Aborted (no chromosome number estimate)"))
    }
  } else {
    print(noquote("   Aborted (file already exists)"))
  }
  execTime <- round(as.numeric(difftime(Sys.time(), startTime, units = "mins")), 2)
  print(noquote(paste0("   ", execTime, " minutes")))
}

