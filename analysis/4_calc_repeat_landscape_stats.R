# Zhaobo Hu
# zhaobohu2002@gmail.com

# Description: Calculates statistics to summarize repeat landscape
# characteristics for each species. saves one file with unparsed
# contigs and another file for parsed contigs

unparsed <- read.csv("../results/vertebrates/unparsed.csv")

# calculate stats
source("functions.R")
files <- list.files("../results/vertebrates/repeat_landscape_divsums")
species <- gsub("_", " ", gsub(".divsum$", "", files))
asmblysz <- unique(unparsed[, c(1, 17)])
asmblysz <- asmblysz[asmblysz$species %in% species, ]
asmblysz <- asmblysz[order(asmblysz$species == species), ]
dat <- data.frame(asmblysz, files)
replandsc.stats <- data.frame()
for (i in dat$species) {
  file <- dat[dat$species == i, ]$files
  asmblysz.Mbp <- dat[dat$species == i, ]$asmblysize.Mbp
  result <- lapply(i, 
                  calcRepLandscStats, 
                  file = file,
                  asmblysz.Mbp = asmblysz.Mbp)[[1]]
  replandsc.stats <- rbind(replandsc.stats, result)
}
df <- merge(unparsed, replandsc.stats, by = "species", all.x = TRUE)

# reorganize and save results
dg <- df[, -c(13:16)]
dg <- data.frame(dg, df[13:16])
write.csv(dg,
          "../results/vertebrates/unparsed.csv", 
          row.names = FALSE)

# parse contigs
dg <- dg[dg$size.Mbp >= 10, ]
sp.lessthanthree <- names(which(table(dg$species) < 3))
dg <- dg[!(dg$species %in% sp.lessthanthree), ]
write.csv(dg,
          "../results/vertebrates/parsed.csv", 
          row.names = FALSE)

