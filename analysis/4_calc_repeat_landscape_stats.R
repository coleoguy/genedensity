# Zhaobo Hu
# zhaobohu2002@gmail.com

# Description: Calculates statistics to summarize repeat landscape
# characteristics for each species. saves one file with unparsed
# contigs and another file for parsed contigs

unparsed <- read.csv("../results/vertebrates/unparsed.csv")

# calculate stats
files <- list.files("../results/vertebrates/repeat_landscape_divsums")
species <- gsub("_", " ", gsub(".divsum$", "", files))
asmblysz <- unique(unparsed[, c(1, 13)])
asmblysz <- asmblysz[asmblysz$species %in% species, ]
asmblysz <- asmblysz[order(asmblysz$species == species), ]
dat <- data.frame(asmblysz, files)
repstats <- data.frame()
for (i in dat$species) {
  file <- dat[dat$species == i, ]$files
  asmblysz.Mbp <- dat[dat$species == i, ]$asmblysize.Mbp
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
  perdivrep.bp <- rowSums(divsum.table[, !names(divsum.table) == "Div"])
  
  # repeat content in Mbp
  totalrep.Mbp <- sum(perdivrep.bp) / 1000000
  
  # repeat content in percent coverage
  perdivrep.pct <- 0.0001 * (perdivrep.bp / asmblysz.Mbp)
  totalrep.pct <- sum(perdivrep.pct)
  
  # median
  median.bin <- which(cumsum(perdivrep.pct) > sum(perdivrep.pct)/2)[1]
  lower <- cumsum(perdivrep.pct)[median.bin-1]
  upper <- cumsum(perdivrep.pct)[median.bin+1]
  mid <- sum(perdivrep.pct)/2
  mediandvg <- median.bin + (mid-lower)/(upper-lower)
  
  # mean
  meandvg <- sum(divergence*perdivrep.pct)/sum(perdivrep.pct)
  
  # build dataframe
  df <- data.frame(i, 
                   totalrep.Mbp,
                   totalrep.pct, 
                   meandvg, 
                   mediandvg)
  repstats <- rbind(repstats, df)
}
colnames(repstats)[1] <- "species"
df <- merge(unparsed, repstats, by = "species", all.x = TRUE)

# reorganize and save results
df <- df[, c(1:22, 27:30, 23:26)]
write.csv(df,
          "../results/vertebrates/unparsed.csv", 
          row.names = FALSE)

# parse contigs
df <- df[df$size.Mbp >= 10, ]
sp.lessthanthree <- names(which(table(df$species) < 3))
df <- df[!(df$species %in% sp.lessthanthree), ]
write.csv(df,
          "../results/vertebrates/parsed.csv", 
          row.names = FALSE)

