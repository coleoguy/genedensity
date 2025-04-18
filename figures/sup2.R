

# largest repeats in each clade
library(data.table)
library(beeswarm)
library(viridis)
dat <- read.csv("../results/parsed.csv")
dat <- dat[dat$thrs == 0.8, ]
dat <- dat[!duplicated(dat$species), ]
dat <- dat[!is.na(dat$chromnum.1n), ]
dat <- na.omit(dat[, c("species", "clade", "rsq")])

files <- list.files("../results/divsums")
masked.sp <- gsub("_", " ", gsub(".divsum$", "", files))

int <- intersect(masked.sp, dat$species)
dat <- dat[dat$species %in% int, ]

largest <- c()
for (i in dat$species) {
  species <- i
  divsum <- readLines(paste0("../results/divsums/", gsub(" ", "_", i), ".divsum"))
  phrase <- "Coverage for each repeat class and divergence (Kimura)"
  start.index <- match(phrase, divsum) + 1
  divsum <- divsum[start.index:length(divsum)]
  divsum <- read.table(textConnection(divsum), 
                       sep = " ", 
                       header = TRUE)
  divsum <- divsum[-c(which(sapply(divsum, function(col) all(is.na(col)))))]
  divsum <- sort(colSums(divsum))
  largest <- c(largest, names(tail(divsum, 1)))
}
dat$largest <- largest
dat <- dat[dat$clade %in% c("Mammalia", "Actinopterygii", "Sauropsida"), ]
dat$clade <- factor(dat$clade, levels = c("Mammalia", "Actinopterygii", "Sauropsida"))
dat$largest[dat$largest %in% c("DNA.TcMar.Tc1", "LINE.L2", "SINE.tRNA.Core")] <- "Others"
map <- c("Mammalia" = "#d95f02", "Actinopterygii" = "#7570b3", "Sauropsida" = "#1b9e77")
cols <- map[dat$clade]
par(mar = c(5, 4, 4, 7)+0.1)
beeswarm(rsq ~ largest, 
         data = dat, 
         xlab = NA, 
         ylab = "R2", 
         pch = 16, 
         pwcol = cols, 
         spacing = 1.1, 
         at = c(1, 1.92, 3, 4.08, 5), 
         labels = c("LINE/CR1", "LINE/L1", "LINE/RTE-BovB", "Others", "Unclassified"))
par(mar = c(5.1, 4.1, 4.1, 2.1))



