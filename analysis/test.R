# hypothesis: fusion reduces gene density variation by combining many
# chromosomes with various gene densities into a few chromosomes with
# similar gene densities. when fission splits a chromosome, it's more
# likely that the new chromosomes are unequal in gene density rather 
# than equal in gene density. If this is true, then smaller chromosomes
# should have more gene density variation than larger chromosomes. 

# For each genome, this plot makes a glm of gene count ~ contig size 
# and uses the absolute values of the model residuals in a glm against 
# contig size. If the hypothesis holds, smaller chromosomes will have
# higher residuals and the slope will be negative.

# Currently, only a few species have models with significant p values 
# and those models all have positive slope. Since many of my genomes are
# not chromosome level, I will look for higher-quality genomes

dat <- read.csv("../results/parsed.csv")
dat <- dat[dat$thrs == 0.8, ]
dat <- dat[dat$clade %in% c("Mammalia", "Sauria", "Actinopterygii"), ]
dat <- dat[!is.na(dat$chromnum.1n), ]
sp <- unique(dat$species)
b <- p <- c()
for (i in sp) {
  sub <- dat[dat$species == i, ]
  model <- glm(genecount ~ size.Mbp, data = sub)
  res <- abs(resid(model))
  b <- c(b, summary(glm(res ~ sub$size.Mbp))$coefficients[2, 1])
  p <- c(p, summary(glm(res ~ sub$size.Mbp))$coefficients[2, 4])
}
df <- data.frame(sp, b, p)