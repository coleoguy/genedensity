
# Zhaobo Hu
# zhaobohu2002@gmail.com

# Description: Parses results and calculates additional statistics
# to summarize contigs for each species

dat <- read.csv("../data/data.csv")

# combine raw contig results
library(data.table)
dir <- "../results/individual_species_results"
files <- paste0(dir, "/",  list.files(dir))
raw <- lapply(files, fread)
raw <- as.data.frame(rbindlist((raw), fill = TRUE))

# parse by contig size
raw <- raw[raw$size.Mbp >= 10, ]

# test new method
parsed <- data.frame()
for (z in unique(dat$species)) {
  sub <- raw[raw$species == z, ]
  total <- sum(sub$size.Mbp)
  est <- dat[dat$species == z, ]$est.gnsz.Mbp
  if (is.na(est) == TRUE) {
    next
  } else {
    if (total/est >= 0.5) {
      parsed <- rbind(parsed, sub)
    }
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
    sdgd <- sd(sub$genedens)
    meangd <- mean(sub$genedens)
    cor <- cor(sub$size.Mbp, sub$genecount)
    weightmean <- sum(sub$genedens * sub$size.Mbp) / sum(sub$size.Mbp)
    weightsd <- sqrt(sum(sub$size.Mbp * (sub$genedens - weightmean)^2) / sum(sub$size.Mbp))
    weightcv <- weightsd / weightmean
    chromsd <- sd(sub$size.Mbp)/mean(sub$size.Mbp)
    contig.stats <- data.frame(species, beta, meangd, sdgd, pval.beta, rsq, cor, weightmean, weightsd, weightcv, chromsd)
  } else {
    beta <- meangd <- sdgd <- pval.beta <- rsq <- cor <- weightmean <- weightsd <- weightcv <- chromsd <- NA
    contig.stats <- data.frame(species, beta, meangd, sdgd, pval.beta, rsq, cor, weightmean, weightsd, weightcv, chromsd)
  }
  final <- rbind(final, merge(merge(dat[dat$species == species, ], contig.stats, by = "species"), sub, by = "species", all = TRUE))
}

#assign clades
final$clade <- final$class
final[final$clade %in% "Aves", ]$clade <- "Sauria"
final[final$clade %in% "Reptilia", ]$clade <- "Sauria"
final[!(final$clade %in% c("Actinopterygii", "Mammalia", "Sauria")), ]$clade <- "Others"

# reorder columns
final <- final[, c(1, 27, 2:11, 26, 12:25)]

# write csv
write.csv(final, "../results/parsed.csv", row.names = FALSE)






