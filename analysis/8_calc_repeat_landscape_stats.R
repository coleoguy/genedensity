# Zhaobo Hu
# zhaobohu2002@gmail.com

# Description: calculates repeat content, k

# "vertebrates" or "invertebrates"
vert.invert <- "invertebrates"

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
