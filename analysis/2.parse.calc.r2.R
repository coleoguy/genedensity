
# Zhaobo Hu
# zhaobohu2002@gmail.com

# Description: Parses results and calculates additional statistics
# to summarize contigs for each species

dat <- read.csv("../data/data.csv")

# combine raw contig results
library(data.table)
# dir <- "../results/individual_species_results"
dir <- "../results/indiv.contigs"
files <- paste0(dir, "/",  list.files(dir))
contigs <- lapply(files, fread)
contigs <- as.data.frame(rbindlist((contigs), fill = TRUE))

# parse by contig size
contigs <- contigs[contigs$size.Mb >= 10, ]

# remove species with less than 3 contigs
rm <- names(table(contigs$species)[table(contigs$species) < 3])
contigs <- contigs[!(contigs$species %in% rm), ]

df <- data.frame()
# contig sum/assembly size ratio threshold
# for (thrs in seq(from = 0, to = 1, by = 0.01)) {
for (thrs in c(0.8, 0.9)) {
# for (thrs in c(0.95)) {
  # remove species if sum of contig sizes is not within some multiple of assembly size
  parsed <- data.frame()
  for (z in unique(contigs$species)) {
    sub <- contigs[contigs$species == z, ]
    cont <- sum(sub$size.Mb)
    total <- contigs[contigs$species == z, ]$asmblysize[1]
    if (cont/total >= thrs) {
      parsed <- rbind(parsed, sub)
    }
  }
  
  # calculate stats based on parsed results
  sp <- unique(parsed$species)
  for (species in sp) {
    sub <- parsed[which(parsed$species == species), ]
    # fit <- summary(glm(sub$genecount ~ sub$size.Mb))
    # beta <- fit$coefficients[2, 1]
    # pval.beta <- fit$coefficients[2, 4]
    rsq <- summary(lm(sub$genecount ~ sub$size.Mb))$r.squared
    # weightmean <- sum(sub$genedens * sub$size.Mb) / sum(sub$size.Mb)
    # weightsd <- sqrt(sum(sub$size.Mb * (sub$genedens - weightmean)^2) / sum(sub$size.Mb))
    # weightcv <- weightsd / weightmean
    stats <- data.frame(species, rsq, thrs)
    sub <- merge(sub, stats, by = "species", all = TRUE)
    sub <- merge(dat[dat$species == species, ], sub, by = "species", all = TRUE)
    df <- rbind(df, sub)
  }
}









#assign clades
df$clade <- df$class
df[df$clade %in% "Aves", ]$clade <- "Sauria"
df[df$clade %in% "Reptilia", ]$clade <- "Sauria"
df[!(df$clade %in% c("Actinopterygii", "Mammalia", "Sauria")), ]$clade <- "Others"

# reorder columns
df <- df[, c(18, 1, 19, 2:11, 16:17, 12:15)]

# write csv
write.csv(df, "../results/parsed.csv", row.names = FALSE)






