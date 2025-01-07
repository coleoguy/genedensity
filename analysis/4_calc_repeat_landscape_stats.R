# Zhaobo Hu
# zhaobohu2002@gmail.com

# Description: Calculates statistics to summarize repeat landscape
# characteristics for each species. saves one file with parsed
# contigs and another file for parsed contigs

# calculate stats
files <- list.files("../results/divsums")
sp <- gsub("_", " ", gsub(".divsum$", "", files))

dat <- read.csv("../results/parsed.csv")
asmbsz <- dat[!duplicated(dat$species), ]
asmbsz <- asmbsz[asmbsz$species %in% sp, ]
asmbsz <- setNames(asmbsz$asmblysize.Mbp*1000000, asmbsz$species)

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
  
  # condense table
  classes <- c("LINE", "SINE", "LTR", "DNA", "RC", "Div", "Unknown")
  for (j in classes) {
    pat <- paste0("^", j, "(\\.|$)")
    headers <- grep(pat, names(table), value = TRUE)
    sub <- table[, headers]
    sums <- rowSums(as.matrix(sub))
    table <- table[, !names(table) %in% headers]
    assign(j, sums)
  }
  Others <- rowSums(as.matrix(table))
  table <- data.frame(Div, LINE, SINE, LTR, DNA, RC, Others, Unknown)
  
  # all repeat total and median
  rep.bp <- rowSums(table[, !names(table) == "Div"])
  total.rep.pct <- sum((rep.bp / asmbsz[i]) * 100)
  total.rep.median <- which(cumsum(rep.bp) > sum(rep.bp)/2)[1]
  
  for (k in classes) {
    assign(paste0(tolower(k), ".rep.pct"), sum(table[k] / asmbsz[i] * 100))
    assign(paste0(tolower(k), ".rep.median"), which(cumsum(table[k]) > sum(table[k])/2)[1])
  }
  
  # build dataframe
  df <- data.frame(species, 
                   total.rep.pct, 
                   total.rep.median, 
                   line.rep.pct, 
                   line.rep.median, 
                   sine.rep.pct, 
                   sine.rep.median, 
                   ltr.rep.pct, 
                   ltr.rep.median, 
                   dna.rep.pct, 
                   dna.rep.median, 
                   rc.rep.pct, 
                   rc.rep.median
                   )
  repstats <- rbind(repstats, df)
}
df <- merge(dat, repstats, by = "species", all.x = TRUE)

# reorganize and save results
df <- df[, c(1:20, 25:36, 21:24)]
write.csv(df,
          "../results/parsed.csv", 
          row.names = FALSE)


