
dat <- read.csv("../results/parsed.csv")
files <- list.files("../results/divsums")
sp <- gsub("_", " ", gsub(".divsum$", "", files))
asmbsz <- dat[!duplicated(dat$species), ]
asmbsz <- asmbsz[asmbsz$species %in% sp, ]
asmbsz <- setNames(asmbsz$asmblysize.Mbp*1000000, asmbsz$species)
repstats <- data.frame()

for (i in 1:length(sp)) {
  species <- sp[i]
  # read text file into lines
  divsum <- readLines(paste0("../results/divsums/", files[i]))
  # look for the start of relevant information
  phrase <- "Coverage for each repeat class and divergence (Kimura)"
  start.index <- match(phrase, divsum) + 1
  # condense relevant lines into a table
  divsum <- divsum[start.index:length(divsum)]
  divsum <- read.table(textConnection(divsum), 
                      sep = " ", 
                      header = TRUE)
  # drop columns with all NA
  divsum <- divsum[-c(which(sapply(divsum, function(col) all(is.na(col)))))]
  
  # condense table
  classes <- c("LINE", "SINE", "LTR", "DNA", "RC", "Div", "Unknown")
  for (j in classes) {
    pat <- paste0("^", j, "(\\.|$)")
    headers <- grep(pat, names(divsum), value = TRUE)
    sub <- divsum[, headers]
    sums <- rowSums(as.matrix(sub))
    divsum <- divsum[, !names(divsum) %in% headers]
    assign(j, sums)
  }
  Others <- rowSums(as.matrix(divsum))
  divsum <- data.frame(Div, LINE, SINE, LTR, DNA, RC, Others, Unknown)
  
  # all repeat total and median
  rep.bp <- rowSums(divsum[, !names(divsum) == "Div"])
  total.rep.pct <- sum((rep.bp / asmbsz[sp[i]]) * 100)
  total.rep.median <- which(cumsum(rep.bp) > sum(rep.bp)/2)[1]
  
  for (k in classes) {
    assign(paste0(tolower(k), ".rep.pct"), sum(divsum[k] /  asmbsz[sp[i]] * 100))
    assign(paste0(tolower(k), ".rep.median"), which(cumsum(divsum[k]) > sum(divsum[k])/2)[1])
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

# reorganize and save results
dat <- merge(dat, repstats, by= "species", all = T)
dat <- dat[, c(1:15, 20:31, 16:19)]
dat <- dat[order(dat$size.Mbp, decreasing = TRUE), ]
dat <- dat[order(dat$species), ]
dat <- dat[order(dat$thrs), ]
write.csv(dat,
          "../results/parsed.csv", 
          row.names = FALSE)


