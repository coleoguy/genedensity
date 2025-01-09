
# Zhaobo Hu
# zhaobohu2002@gmail.com

# Description: Parses results and calculates additional statistics
# to summarize contigs for each species

# contig sum/assembly size ratio threshold
thrs <- 0

dat <- read.csv("../data/data.csv")

# combine raw contig results
library(data.table)
dir <- "../results/individual_species_results"
files <- paste0(dir, "/",  list.files(dir))
contigs <- lapply(files, fread)
contigs <- as.data.frame(rbindlist((contigs), fill = TRUE))

# parse by contig size
contigs <- contigs[contigs$size.Mbp >= 10, ]

# remove species with less than 2 contigs
rm <- names(table(contigs$species)[table(contigs$species) < 2])
contigs <- contigs[!(contigs$species %in% rm), ]

# test new method
parsed <- data.frame()
for (z in unique(contigs$species)) {
  sub <- contigs[contigs$species == z, ]
  cont <- sum(sub$size.Mbp)
  total <- contigs[contigs$species == z, ]$asmblysize[1]
  if (cont/total >= thrs) {
    parsed <- rbind(parsed, sub)
  }
}

# calculate stats based on parsed results
sp <- unique(parsed$species)
final <- data.frame()
for (species in sp) {
  i <- species
  sub <- parsed[which(parsed$species == i), ]
  if (nrow(sub) > 0){
    fit <- summary(glm(sub$genecount ~ sub$size.Mbp))
    beta <- fit$coefficients[2, 1]
    pval.beta <- fit$coefficients[2, 4]
    rsq <- summary(lm(sub$genecount ~ sub$size.Mbp))$r.squared
    weightmean <- sum(sub$genedens * sub$size.Mbp) / sum(sub$size.Mbp)
    weightsd <- sqrt(sum(sub$size.Mbp * (sub$genedens - weightmean)^2) / sum(sub$size.Mbp))
    weightcv <- weightsd / weightmean
    contig.stats <- data.frame(species, beta, pval.beta, rsq, weightmean, weightsd, weightcv)
  } else {
    beta <- pval.beta <- rsq <- weightmean <- weightsd <- weightcv <- NA
    contig.stats <- data.frame(species, beta, pval.beta, rsq, weightmean, weightsd, weightcv)
  }
  final <- rbind(final, merge(merge(dat[dat$species == species, ], contig.stats, by = "species"), sub, by = "species", all = TRUE))
}

#assign clades
final$clade <- final$class
final[final$clade %in% "Aves", ]$clade <- "Sauria"
final[final$clade %in% "Reptilia", ]$clade <- "Sauria"
final[!(final$clade %in% c("Actinopterygii", "Mammalia", "Sauria")), ]$clade <- "Others"

final$cont.asmb.rat.cutoff <- thrs

# reorder columns
final <- final[, c(1, 23, 2:8, 12:17, 22, 24, 9:11, 18:21)]

# write csv
write.csv(final, "../results/parsed.csv", row.names = FALSE)






