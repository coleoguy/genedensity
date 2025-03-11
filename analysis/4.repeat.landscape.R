
dat <- read.csv("../results/parsed.csv")
files <- list.files("../results/divsums")
sp <- gsub("_", " ", gsub(".divsum$", "", files))
asmbsz <- dat[!duplicated(dat$species), ]
asmbsz <- asmbsz[asmbsz$species %in% sp, ]
asmbsz <- setNames(asmbsz$asmblysize.Mb*1000000, asmbsz$species)
repstats <- data.frame()

lis <- vector("list", length(sp))
for (h in 1:length(sp)) {
  species <- sp[h]
  # read text file into lines
  divsum <- readLines(paste0("../results/divsums/", files[h]))
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
  lis[[h]] <- colnames(divsum)
}

to.collapse <- sort(unique(
  sub("\\..*", "", unique(grep("\\.", unlist(lis), value = TRUE)))
  ))
others <- unique(grep("\\.", unlist(lis), value = TRUE, invert = TRUE))
others <- others[!others %in% to.collapse]
others <- others[!others %in% "Div"]

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
  
  # collapse
  for (j in to.collapse) {
    pat <- paste0("^", j, "(\\.|$)")
    headers <- grep(pat, names(divsum), value = TRUE)
      sub <- divsum[, headers]
      sums <- rowSums(as.matrix(sub))
      #divsum <- divsum[, !names(divsum) %in% headers]
      assign(j, sums)
  }
  
  # record others
  for (k in others) {
    pat <- paste0("^", k, "(\\.|$)")
    headers <- grep(pat, names(divsum), value = TRUE)
    sub <- divsum[, headers]
    sums <- rowSums(as.matrix(sub))
    assign(k, sums)
  }
  
  div <- divsum$Div
  divsum <- data.frame(div, mget(to.collapse), mget(others))
  
  # collapse further
  to.keep <- c("div", "DNA", "LINE", "LTR", "SINE", "Unknown")
  to.sum <- divsum[, colnames(divsum)[!colnames(divsum) %in% to.keep]]
  divsum <- divsum[, to.keep]
  divsum$others <- rowSums(to.sum)
  
  df <- data.frame(species)
  for (l in colnames(divsum)[-1]) {
    # proportion
    df[[paste0("prop.", tolower(l))]] <- sum(divsum[l] / asmbsz[sp[i]])

    # age
    df[[paste0("age.", tolower(l))]] <- which(cumsum(divsum[l]) > sum(divsum[l])/2)[1]
  }
  repstats <- rbind(repstats, df)
}

# reorganize and save results
dat <- merge(dat, repstats, by= "species", all = T)
dat <- dat[order(dat$size.Mb, decreasing = TRUE), ]
dat <- dat[order(dat$species), ]
dat <- dat[order(dat$thrs), ]
write.csv(dat,
          "../results/parsed.csv", 
          row.names = FALSE)




