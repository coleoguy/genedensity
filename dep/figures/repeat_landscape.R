
# repeat landscape
par(mar = c(3, 4, 1, 7)+0.1)
dat <- read.csv("../results/parsed.csv")
dat <- dat[dat$thrs == 0, ]
dat <- dat[!duplicated(dat$species), ]
dat <- na.omit(dat[, c("species", "rsq", "clade", "asmblysize.Mbp")])
files <- list.files("../results/divsums")
sp <- gsub("_", " ", gsub(".divsum$", "", files))
asmbsz <- dat[dat$species %in% sp, ]
asmbsz <- setNames(asmbsz$asmblysize.Mbp*1000000, asmbsz$species)
cols <- c("#e31a1c", "#fb9a99", "#33a02c", "#b2df8a", "#1f78b4", "#a6cee3")
for (i in 1:length(sp)) {
  species <- sp[i]
  divsum <- readLines(paste0("../results/divsums/", files[i]))
  phrase <- "Coverage for each repeat class and divergence (Kimura)"
  start.index <- match(phrase, divsum) + 1
  divsum <- divsum[start.index:length(divsum)]
  divsum <- read.table(textConnection(divsum), 
                       sep = " ", 
                       header = TRUE)
  divsum <- divsum[-c(which(sapply(divsum, function(col) all(is.na(col)))))]
  classes <- c("LINE", "SINE", "LTR", "DNA", "Div", "Unknown")
  for (j in classes) {
    pat <- paste0("^", j, "(\\.|$)")
    headers <- grep(pat, names(divsum), value = TRUE)
    sub <- divsum[, headers]
    sums <- rowSums(as.matrix(sub)) / asmbsz[i] * 100
    divsum <- divsum[, !names(divsum) %in% headers]
    assign(j, sums)
  }
  Others <- rowSums(as.matrix(divsum)) / asmbsz[i] * 100
  divsum <- data.frame(LINE, SINE, LTR, DNA, Others, Unknown)
  divsum <- t(as.matrix(divsum))
  divsum <- divsum[, 1:51]
  divsum <- divsum[nrow(divsum):1, ]
  svg(paste0("repeat_landscape/", gsub(" ", "_", species), ".svg"), 
      width = 8.6, 
      height = 6.6)
  barplot(divsum, 
          col = cols, 
          space = 0, 
          border = NA, 
          xlab = NA, 
          ylab = NA, 
          ylim = c(0, 1.15*max(colSums(divsum))))
  axis(1)
  pos <- seq(0.2 * max(colSums(divsum)), 0.85 * max(colSums(divsum)), length.out = 6)
  text(55, pos[6], "LINE", xpd = NA, adj = c(0, 0.5))
  text(55, pos[5], "SINE", xpd = NA, adj = c(0, 0.5))
  text(55, pos[4], "LTR", xpd = NA, adj = c(0, 0.5))
  text(55, pos[3], "DNA", xpd = NA, adj = c(0, 0.5))
  text(55, pos[2], "Others", xpd = NA, adj = c(0, 0.5))
  text(55, pos[1], "Unclassified", xpd = NA, adj = c(0, 0.5))
  points(53, pos[6], pch = 15, col = cols[6], xpd = NA, cex = 1.5)
  points(53, pos[5], pch = 15, col = cols[5], xpd = NA, cex = 1.5)
  points(53, pos[4], pch = 15, col = cols[4], xpd = NA, cex = 1.5)
  points(53, pos[3], pch = 15, col = cols[3], xpd = NA, cex = 1.5)
  points(53, pos[2], pch = 15, col = cols[2], xpd = NA, cex = 1.5)
  points(53, pos[1], pch = 15, col = cols[1], xpd = NA, cex = 1.5)
  dev.off()
}
par(mar = c(5.1, 4.1, 4.1, 2.1))
