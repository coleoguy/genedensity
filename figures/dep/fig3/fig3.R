



# filter data
library(data.table)
library(beeswarm)
library(viridis)
dat <- read.csv("../../results/rsq.csv")
dat <- na.omit(dat[, c("species", "clade", "rsq", "assem.sz")])
dat$assem.sz <- dat$assem.sz * 1000000
files <- list.files("../../results/divsums")
masked.sp <- gsub(".divsum$", "", files)
int <- intersect(masked.sp, dat$species)
dat <- dat[dat$species %in% int, ]
dat <- dat[dat$species %in% c("Bos_taurus", "Peromyscus_maniculatus_bairdii"), ]

# repeat landscape
cols <- c("#e31a1c", "#fb9a99", "#33a02c", "#b2df8a", "#1f78b4", "#a6cee3")
# par(mar = c(3, 4, 1, 7)+0.1)
par(mfrow = c(1, 2))
for (i in dat$species) {
  if (i == dat$species[1]) {
    par(mar = c(5.1, 4.1, 4.1, 0.8))
  } else {
    par(mar = c(5.1, 0.7, 4.1, 4.2))
  }
  species <- i
  divsum <- readLines(paste0("../../results/divsums/", i, ".divsum"))
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
    sums <- rowSums(as.matrix(sub)) / dat[dat$species == i, ]$assem.sz * 100
    divsum <- divsum[, !names(divsum) %in% headers]
    assign(j, sums)
  }
  Others <- rowSums(as.matrix(divsum)) / dat[dat$species == i, ]$assem.sz * 100
  
  divsum <- data.frame(LINE, SINE, LTR, DNA, Others, Unknown)
  divsum <- t(as.matrix(divsum))
  divsum <- divsum[, 1:51]
  divsum <- divsum[nrow(divsum):1, ]
  if (i == dat$species[1]) {
    barplot(divsum, 
            col = cols, 
            space = 0, 
            border = NA, 
            xlab = NA, 
            ylab = NA, 
            ylim = c(0, 1.1*3.049599))
    axis(1)
  } else {
    barplot(divsum, 
            col = cols, 
            space = 0, 
            border = NA, 
            axes = FALSE, 
            xlab = NA, 
            ylab = NA, 
            ylim = c(0, 1.1*3.049599))
    axis(1)
    
    pos <- seq(0.3 * 3.049599, 0.95 * 3.049599, length.out = 6)
    text(47, pos[6], "LINE", xpd = NA, adj = c(0, 0.5))
    text(47, pos[5], "SINE", xpd = NA, adj = c(0, 0.5))
    text(47, pos[4], "LTR", xpd = NA, adj = c(0, 0.5))
    text(47, pos[3], "DNA", xpd = NA, adj = c(0, 0.5))
    text(47, pos[2], "Others", xpd = NA, adj = c(0, 0.5))
    text(47, pos[1], "Unclassified", xpd = NA, adj = c(0, 0.5))
    
    points(45, pos[6], pch = 15, col = cols[6], xpd = NA, cex = 1.5)
    points(45, pos[5], pch = 15, col = cols[5], xpd = NA, cex = 1.5)
    points(45, pos[4], pch = 15, col = cols[4], xpd = NA, cex = 1.5)
    points(45, pos[3], pch = 15, col = cols[3], xpd = NA, cex = 1.5)
    points(45, pos[2], pch = 15, col = cols[2], xpd = NA, cex = 1.5)
    points(45, pos[1], pch = 15, col = cols[1], xpd = NA, cex = 1.5)
  }

}
par(mfrow = c(1, 1))
par(mar = c(5.1, 4.1, 4.1, 2.1))

