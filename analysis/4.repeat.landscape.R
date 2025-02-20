
dat <- read.csv("../results/parsed.csv")
files <- list.files("../results/divsums")
sp <- gsub("_", " ", gsub(".divsum$", "", files))
asmbsz <- dat[!duplicated(dat$species), ]
asmbsz <- asmbsz[asmbsz$species %in% sp, ]
asmbsz <- setNames(asmbsz$asmblysize.Mb*1000000, asmbsz$species)
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
  col.to.keep <- c("DNA", "LINE", "LTR", "SINE", "Div", "Unknown")
  for (j in col.to.keep) {
    pat <- paste0("^", j, "(\\.|$)")
    headers <- grep(pat, names(divsum), value = TRUE)
    sub <- divsum[, headers]
    sums <- rowSums(as.matrix(sub))
    divsum <- divsum[, !names(divsum) %in% headers]
    assign(j, sums)
  }
  Others <- rowSums(as.matrix(divsum))
  divsum <- data.frame(Div, DNA, LINE, LTR, SINE, Others, Unknown)
  rep.bp <- rowSums(divsum[, !names(divsum) == "Div"])
  
  # set which repeats to record
  rep <- colnames(divsum)[!colnames(divsum) %in% c("Div")]
  
  # find repeat proportion
  df <- data.frame(species)
  for (k in rep) {
    name <- paste0("rep.prop.", tolower(k))
    assign(name, sum(divsum[k] /  asmbsz[sp[i]]))
    df[[name]] <- get(name)
  }
  
  # find repeat age
  rep <- unlist(lapply(1:length(rep), function(x) {
    apply(combn(rep, x), 2, paste, collapse = ".")
  }))
  for (l in rep) {
    sub <- rowSums(divsum[unlist(strsplit(l, "\\."))])
    name <- paste0("rep.age.", tolower(l))
    assign(name, which(cumsum(sub) > sum(sub)/2)[1])
    df[[name]] <- get(name)
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




