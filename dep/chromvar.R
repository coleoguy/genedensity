
# does gene density variation across chromosomes in a genome vary less
# as chromosomes get larger?

# for each species:
# 1. make gls(gene number ~ chromosome size)
# 2. take absolute residuals of model
# 3. make gls(abs(resid) ~ chromosome size)
# 4. record beta and p

# only 10 species show significant relationships between residuals and 
# chromosome size. in all 10 species, gene density variation increases 
# as chromosomes get bigger

dat <- read.csv("../results/parsed.csv")
dat <- dat[!is.na(dat$chromnum.1n), ]
df <- data.frame()
for (i in unique(dat$species)) {
  species <- i
  sub <- dat[dat$species == i, ]
  clade <- unique(sub$clade)
  res <- resid(glm(sub$genecount ~ sub$size.Mbp))
  b <- summary(glm(abs(res) ~ sub$size.Mbp))$coefficients[2, 1]
  p <- summary(glm(abs(res) ~ sub$size.Mbp))$coefficients[2, 4]
  df <- rbind(df, data.frame(species, clade, b, p))
}
df <- df[df$p < 0.05, ]
