# Zhaobo Hu
# zhaobohu2002@gmail.com

# Description: Calculates statistics to summarize repeat landscape
# characteristics for each species. saves one file with parsed
# contigs and another file for parsed contigs

parsed <- read.csv("../results/vertebrates/parsed.csv")

# calculate stats
files <- list.files("../results/vertebrates/repeat_landscape_divsums")
species <- gsub("_", " ", gsub(".divsum$", "", files))
asmblysz <- unique(parsed[, c(1, 13)])
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
  # vector of the repeat content of each divergence score
  rep.bp <- rowSums(divsum.table[, !names(divsum.table) == "Div"])
  
  # repeat content in Mbp
  totalrep.Mbp <- sum(rep.bp) / 1000000
  
  # repeat content in percent coverage
  rep.pct <- 0.0001 * (rep.bp / asmblysz.Mbp) # 0.0001% = (bp / Mbp) * (1 Mbp / 1000000 bp) * (100%) 
  totalrep.pct <- sum(rep.pct)
  
  # divergence bin with median repeat
  median <- which(cumsum(rep.pct) > sum(rep.pct)/2)[1]
  
  # unused
  # mean <- sum(divergence*perdivrep.pct)/sum(perdivrep.pct)
  # k <- kurtosis(perdivrep.pct)
  # s <- skewness(perdivrep.pct)
  # max <- max(perdivrep.pct)
  # which <- which.max(perdivrep.pct)
  
  # build dataframe
  df <- data.frame(i, 
                   totalrep.Mbp,
                   totalrep.pct, 
                   median)
  repstats <- rbind(repstats, df)
}
colnames(repstats)[1] <- "species"
df <- merge(parsed, repstats, by = "species", all.x = TRUE)

# reorganize and save results
df <- df[, c(1:22, 27:29, 23:26)]
write.csv(df,
          "../results/vertebrates/parsed.csv", 
          row.names = FALSE)


