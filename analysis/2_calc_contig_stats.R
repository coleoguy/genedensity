
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

# parse results
parsed <- raw[raw$size.Mbp >= 10, ]
sp.lessthanthree <- names(which(table(parsed$species) < 3))
parsed <- parsed[!(parsed$species %in% sp.lessthanthree), ]

# calculate stats based on parsed results and record the unparsed contigs
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
    contig.stats <- data.frame(species, beta, meangd, sdgd, pval.beta, rsq, cor, weightmean, weightsd, weightcv)
  } else {
    beta <- meangd <- sdgd <- pval.beta <- rsq <- cor <- weightmean <- weightsd <- weightcv <- NA
    contig.stats <- data.frame(species, beta, meangd, sdgd, pval.beta, rsq, cor, weightmean, weightsd, weightcv)
  }
  final <- rbind(final, merge(merge(dat[dat$species == species, ], contig.stats, by = "species"), sub, by = "species", all = TRUE))
}

#assign clades
final$clade <- final$class
final[final$clade %in% "Aves", ]$clade <- "Sauria"
final[final$clade %in% "Reptilia", ]$clade <- "Sauria"
final[!(final$clade %in% c("Actinopterygii", "Mammalia", "Sauria")), ]$clade <- "Others"

# reorder columns
final <- final[, c(1, 26, 2:11, 25, 12:24)]

# write csv
write.csv(final, "../results/parsed.csv", row.names = FALSE)






