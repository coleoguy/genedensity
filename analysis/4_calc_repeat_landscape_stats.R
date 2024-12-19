# Zhaobo Hu
# zhaobohu2002@gmail.com

# Description: Calculates statistics to summarize repeat landscape
# characteristics for each species. saves one file with parsed
# contigs and another file for parsed contigs

# calculate stats
# library(e1071)
files <- list.files("../results/divsums")
sp <- gsub("_", " ", gsub(".divsum$", "", files))

dat <- read.csv("../results/parsed.csv")
asmbsz <- dat[!duplicated(dat$species), ]
asmbsz <- asmbsz[asmbsz$species %in% sp, ]
asmbsz <- asmbsz[order(asmbsz$species == sp), ]
asmbsz <- asmbsz$asmblysize.Mbp*1000000

repstats <- data.frame()
for (i in 1:length(sp)) {
  species <- sp[i]
  # read text file into lines
  lines <- readLines(paste0("../results/divsums/", files[i]))
  # look for the start of relevant information
  phrase <- "Coverage for each repeat class and divergence (Kimura)"
  start.index <- match(phrase, lines) + 1
  # condense relevant lines into a table
  lines <- lines[start.index:length(lines)]
  table <- read.table(textConnection(lines), 
                             sep = " ", 
                             header = TRUE)
  # drop NA columns
  table <- table[-c(which(sapply(table, function(col) all(is.na(col)))))]
  # vector of divergence scores
  div <- table$Div
  # vector of the repeat content of each divergence score
  rep.bp <- rowSums(table[, !names(table) == "Div"])
  # repeat content in Mbp
  totalrep.Mbp <- sum(rep.bp) / 1000000
  # repeat content in percent coverage
  rep.pct <- 0.0001 * (rep.bp / asmbsz[i]) # 0.0001% = (bp / Mbp) * (1 Mbp / 1000000 bp) * (100%) 
  totalrep.pct <- sum(rep.pct)
  # divergence bin with median repeat
  median <- which(cumsum(rep.pct) > sum(rep.pct)/2)[1]
  
  # unused
  # mean <- sum(divergence*rep.pct)/sum(rep.pct)
  # k <- kurtosis(rep.pct)
  # s <- skewness(rep.pct)
  # max <- max(rep.pct)
  # which <- which.max(rep.pct)
  # s <- skewness(rep.pct)
  # k <- kurtosis(rep.pct)
  
  # build dataframe
  df <- data.frame(species, 
                   totalrep.Mbp,
                   totalrep.pct, 
                   median
                   )
  repstats <- rbind(repstats, df)
}
df <- merge(dat, repstats, by = "species", all.x = TRUE)

# reorganize and save results
df <- df[, c(1:23, 28:30, 24:27)]
write.csv(df,
          "../results/parsed.csv", 
          row.names = FALSE)


