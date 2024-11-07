# Zhaobo Hu
# zhaobohu2002@gmail.com

# Description: for each species in the divergence summary folder, pulls the 
# coverage for each repeat class and divergence for use in plotting

# "vertebrates" or "invertebrates"?
vert.invert <- "invertebrates"
# verbose
verbose <- F

asmblysz <- read.csv(paste0("../results/", vert.invert, "/assembly_sizes.csv"))
files <- list.files(paste0("../results/", vert.invert, "/repeat_landscape_divsums"))
for (file in files) {
  # species name
  species <- gsub("_", " ", gsub(".divsum$", "", file))
  if (verbose == TRUE) {
    start.time <- Sys.time()
    gc()
    print(noquote(paste0(species, " (", Sys.time(), ")")))
  }
  # read text file into lines
  divsum.vector <- readLines(paste0("../results/", 
                                    vert.invert, 
                                    "/repeat_landscape_divsums/", 
                                    file))
  # look for the start of relevant information
  phrase <- "Coverage for each repeat class and divergence (Kimura)"
  # get index of the starting line
  start.index <- match(phrase, divsum.vector) + 1
  # condense lines into a table
  divsum.vector <- divsum.vector[start.index:length(divsum.vector)]
  divsum.table <- read.table(textConnection(divsum.vector), sep = " ", header = TRUE)
  # drop NA columns
  divsum.table <- divsum.table[-c(which(sapply(divsum.table, 
                                               function(col) all(is.na(col)))))]
  # get repeat categories from column names
  repeat.cats <- sapply(colnames(divsum.table), 
                        function(names) strsplit(names, "\\.")[[1]][1])
  # find non-unique repeat categories among column names
  non.unique.cats <- names(table(repeat.cats)[table(repeat.cats) > 1])
  # get the indices of the non-unique colnames
  non.unique.indices <- sapply(non.unique.cats, 
                               function(str) which(as.character(repeat.cats) == str))
  # get the row sums of columns with non-unique names
  non.unique.sum <- lapply(non.unique.indices, 
                           function(index) rowSums(divsum.table[c(index)]))
  # remove non-unique columns from table
  results.table <- divsum.table[-c(as.numeric(unlist(non.unique.indices)))]
  # rename table
  colnames(results.table) <- repeat.cats[-c(as.numeric(unlist(non.unique.indices)))]
  # add summed columns to table
  results.table <- data.frame(results.table, non.unique.sum)
  # sort repeat categories by length
  cats.sort.length <- sort(colSums(results.table), decreasing = TRUE)
  # find categories that have been classified
  known.cats <- names(sort(colSums(results.table), decreasing = TRUE)) != "Unknown"
  # find the 4 largest repeat groups
  top.4.cat.names <- names(head(cats.sort.length[known.cats], 4))
  # find indices of the 4 largest repeat groups
  top.4.cat.ind <- which(colnames(results.table) %in% top.4.cat.names)
  # find index of uncharacterized column
  unknown.ind <- which(colnames(results.table) %in% "Unknown")
  # find index of divergence column
  div.ind <- which(colnames(results.table) %in% "Div")
  # find and sum columns that don't belong to any aforementioned categories
  other <- rowSums(results.table[-c(div.ind, top.4.cat.ind, unknown.ind)])
  # subset table with divergence, top 4 classes, and uncharacterized
  results.table <- results.table[c(div.ind, top.4.cat.ind, unknown.ind)]
  # add "Other" to the table
  results.table$Other <- other
  
  # ------------------------format results for plotting------------------------
  
  # format repeat divergence for plotting
  divergence <- rep(results.table$Div, length(colnames(results.table))-1)
  # get sequence lengths from the table as a vector
  sequence.lengths.Mbp <- as.numeric(unlist(results.table[-c(1)])) / 1000000
  # find the assembly size of the species
  gnsz.Mbp <- asmblysz$asmbly.size.Mbp[asmblysz$species == species]
  # calculate the percentage of the genome for repeats in each repeat class and divergence
  percent.of.genome <- sequence.lengths.Mbp / gnsz.Mbp * 100
  # format repeat classes for plotting
  repeat.group <- as.character(
    sapply(colnames(results.table)[colnames(results.table) != "Div"], 
           function(rowname) rep(rowname, length(results.table$Div))))
  # assemble dataframe for plotting
  final <- data.frame(species, divergence, sequence.lengths.Mbp, percent.of.genome, repeat.group)
  write.csv(final, 
            paste0("../results/vertebrates/repeat_landscape_plotting/", 
                   gsub(" ", "_", species), 
                   ".csv"),
            row.names = FALSE)
  if (verbose == TRUE) {
    print(noquote("   Successfully written!"))
  }
}





