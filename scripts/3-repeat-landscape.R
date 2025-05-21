
dat <- read.csv("../results/rsq.csv")
files <- list.files("../results/divsums")
sp <- gsub(".divsum$", "", files)
repstats <- data.frame()

lis <- vector("list", length(sp))
for (h in 1:length(sp)) {
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
  asmbsz <- dat[dat$species == species, ]$assem.sz * 1000000
  if (length(asmbsz) == 0L) {
    asmbsz <- NA
  }
  
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
    vec <- divsum[[l]]
    
    # proportion
    df[[paste0("prop.", tolower(l))]] <- sum(vec / asmbsz)
    
    # age
    idx <- seq(vec)-1
    # vec[vec <= 0] <- 0
    if(any(vec[vec <= 0])) {
      df[[paste0("age.", tolower(l))]] <- NA
      next
    }
    values <- rep(idx, vec)
    df[[paste0("age.", tolower(l))]] <- median(values)
    rm(values)
    gc()
  }
  repstats <- rbind(repstats, df)
  print(species)
}

write.csv(repstats, "../results/repeat-results.csv", row.names = FALSE)


