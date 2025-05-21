



# filter data
library(data.table)
library(beeswarm)
library(viridis)
dat <- read.csv("../results/parsed.csv")
dat <- dat[dat$thrs == 0.8, ]
dat <- dat[!duplicated(dat$species), ]
dat <- dat[!is.na(dat$chromnum.1n), ]
dat <- na.omit(dat[, c("species", "clade", "rsq", "asmblysize.Mb")])
names(dat)[names(dat) == "asmblysize.Mb"] <- "asmblysize.bp"
dat$asmblysize.bp <- dat$asmblysize.bp * 1000000
files <- list.files("../results/divsums")
masked.sp <- gsub("_", " ", gsub(".divsum$", "", files))
int <- intersect(masked.sp, dat$species)
dat <- dat[dat$species %in% int, ]
dat <- dat[dat$species %in% c("Bos taurus", "Peromyscus maniculatus bairdii"), ]

# repeat landscape
cols <- c("#e31a1c", "#fb9a99", "#33a02c", "#b2df8a", "#1f78b4", "#a6cee3")
# par(mar = c(3, 4, 1, 7)+0.1)
# par(mfrow = c(1, 2))
i <- dat$species[1]
if (i == dat$species[1]) {
  par(mar = c(5.1, 4.1, 4.1, 0.8))
} else {
  par(mar = c(5.1, 0.7, 4.1, 4.2))
}
species <- i
divsum <- readLines(paste0("../results/divsums/", gsub(" ", "_", i), ".divsum"))
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
  sums <- rowSums(as.matrix(sub)) / dat[dat$species == i, ]$asmblysize.bp * 100
  divsum <- divsum[, !names(divsum) %in% headers]
  assign(j, sums)
}
Others <- rowSums(as.matrix(divsum)) / dat[dat$species == i, ]$asmblysize.bp * 100
# divsum <- data.frame(LINE, SINE, LTR, DNA, Others, Unknown)
divsum <- data.frame(LINE)
divsum <- t(as.matrix(divsum))
divsum <- divsum[, 1:51]
# divsum <- divsum[nrow(divsum):1, ]
barplot(divsum, 
        # col = cols, 
        col = "#a6cee3",
        space = 0, 
        border = NA, 
        xlab = NA, 
        ylab = NA, 
        #yaxt = "n",
        ylim = c(0, 1.1*max(divsum)))
#axis(2, at = round(c(0, max(divsum)), 2))
axis(1)


par(mfrow = c(1, 1))
par(mar = c(5.1, 4.1, 4.1, 2.1))

