# Zhaobo Hu
# zhaobohu2002@gmail.com

# Description: calculates repeat content, k

# "vertebrates" or "invertebrates"
vert.invert <- "vertebrates"

source("functions.R")
files <- list.files(paste0("../results/", 
                           vert.invert, 
                           "/repeat_landscape_divsums"))

# species name
species <- gsub("_", " ", gsub(".divsum$", "", files))

asmblysz <- read.csv(paste0("../results/", vert.invert, "/assembly_sizes.csv"))
asmblysz <- asmblysz[asmblysz$species %in% species, ]
asmblysz <- asmblysz[order(asmblysz$species == species), ]

dat <- data.frame(asmblysz, files)

replandsc.stats <- data.frame()

for (i in dat$species) {
  file <- dat[dat$species == i, ]$files
  asmblysz.Mbp <- dat[dat$species == i, ]$asmbly.size.Mbp
  result <- lapply(i, 
                  calcRepLandscStats, 
                  file = file,
                  asmblysz.Mbp = asmblysz.Mbp,
                  vert.invert = vert.invert)[[1]]
  replandsc.stats <- rbind(replandsc.stats, result)
}

write.csv(replandsc.stats, 
          paste0("../results/", vert.invert, "/rep_landscape_stats.csv"), 
          row.names = FALSE)









######### ignore stuff below
library(pracma)

all_peaks <- findpeaks(frequency)
all_mins <- findpeaks(-frequency)
all_mins[, 1] <- -all_mins[, 1]
plot(frequency, type = 'l', main = "Frequency Data with Detected Peaks/Mins", 
     xlab = "Index", ylab = "Amplitude")
points(all_peaks[, 2], frequency[all_peaks[, 2]], col = "red", pch = 19)
points(all_mins[, 2], all_mins[, 1], col = "blue", pch = 19)
