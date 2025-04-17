




# Snippets
# list result files
results <- list.files(
  paste0("../results/individual_species_results"))
# pull species names from list of result files
result.species <- gsub("_", " ", unique(gsub("\\..*$", "", results)))


if (verbose == TRUE) {
  start.time <- Sys.time()
  gc()
  print(noquote(paste0(species, " (", Sys.time(), ")")))
}
# proceed if result file is not found, else go to next species
if (!(species %in% result.species)) {
  # fasta.path, chromnum.1n, mito.keywords, verbose, ensembl)
  
  # create full paths for files. sort alphabetically to place
  # the fasta file in front of the gff3/gtf file
  genome.files <- sort(paste0("../data/genomes/", 
                              list.files(paste0("../data/genomes"), 
                                         gsub(" ", "_", species))))
  # write results
  result.path <- paste0("../results/individual_species_results/", 
                        gsub(" ", "_", species), 
                        ".csv")
  fwrite(dat, file = result.path, row.names = FALSE)
  if (verbose == TRUE) {
    print(noquote("   Successfully written!"))
  }
  # keywords to match and exclude mitochondrial contigs for assembly size calcs
  mito.keywords <- c("mt", "mito", "mitochondrial", "nonchromosomal")
  
  
  
  
  
  
  # reorganize and save results TODO should produce file where each row is a 
  # species and each column is either a proportion or a age for a kind of repeat
  dat <- merge(dat, repstats, by= "species", all = T)
  dat <- dat[order(dat$size.Mb, decreasing = TRUE), ]
  dat <- dat[order(dat$species), ]
  dat <- dat[order(dat$thrs), ]
  
  
  
  permTest <- function(x, y, reps, method) {
    permuted.y <- replicate(reps, sample(y))
    permuted.cor <- apply(permuted.y, MARGIN = 2, function(col) cor(x, col, method = method))
    pval <- mean(abs(permuted.cor) >= abs(cor(x, y, method = method)))
    return(pval)
  }
  