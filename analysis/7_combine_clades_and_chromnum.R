

chromnums <- read.csv("../data/vertebrates/chromnums.csv")
clades.gnsz <- read.csv("../data/vertebrates/clades_gnsz.csv")
parsed.results <- read.csv("../results/vertebrates/parsed_results.csv")
species <- unique(clades.gnsz$species)
rsq <- unique(parsed.results$species.rsquared[parsed.results$species %in% species])
dat <- data.frame(
  chromnums[chromnums$species %in% species, ],
  rsq,
  clades.gnsz[-1]
)
dat$genome.size.est_bp <- dat$genome.size.est_bp / 1000000000
colnames(dat)[colnames(dat) == "genome.size.est_bp"] <- "gnsz_Gbp"

write.csv(dat, "../data/vertebrates/chromnum_clade_gnsz_rsq.csv", row.names = FALSE)
