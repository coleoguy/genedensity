
# why do reptiles have low gene density variation if microchromosomes are 
# more gene dense than normal chromosomes?

# a plot of gene density vs chromosome size for all reptile chromosomes shows a 
# decreasing monotonic pattern with increasing slope. the smallest chromosomes have 
# 100-150 genes per Mb and this gene density decreases quickly with small increases
# in chromosome size. When chromosomes are large, even a large change in chromosome
# size would only produce a small change in gene density. Although microchromosomes
# can have up to 10x gene density of main chromosomes, linear models of 
# gene count ~ chromosome size treat this gene density as noise because the actual number
# of genes found in very small chromosomes are very small regardless of how gene
# dense the chromosomes are. 

dat <- read.csv("../results/parsed.csv")
dat <- dat[dat$thrs == 0.9, ]
dat <- dat[!is.na(dat$species), ]
dat <- dat[!is.na(dat$chromnum.1n), ]
rep <- dat[dat$clade == "Sauria", ]
mam <- dat[dat$clade == "Mammalia", ]
plot(dat$size.Mbp, dat$genedens)
plot(rep$size.Mbp, rep$genedens, xlim = c(0, 800), ylim = c(0, 160))
plot(mam$size.Mbp, mam$genedens, xlim = c(0, 800), ylim = c(0, 160))

# larger genomes have more variation in chromosome size
i <- "Numida meleagris"
sub <- dat[dat$species == i, ]
plot(sub$size.Mbp, sub$genecount, xlim = c(0, 800))
abline(glm(sub$genecount ~ sub$size.Mbp))


lis <- list()
for (i in unique(dat$species)) {
  sub <- dat[dat$species == i, ]
  cl <- unique(sub$clade)
  b <- summary(glm(sub$genecount ~ sub$size.Mbp))$coefficients[2, 1]
  p <- summary(glm(sub$genecount ~ sub$size.Mbp))$coefficients[2, 4]
  sizerange.Mb <- diff(range(sub$size.Mbp))
  gnsz.Mb <- unique(sub$asmblysize.Mbp)
  lis <- c(lis, list(c(i, cl, b, p, gnsz.Mb, sizerange.Mb)))
}
df <- as.data.frame(do.call(rbind, lis))
names(df) <- c("species", "cl", "beta", "p", "gnsz.Mb", "sizerange.Mb")
num <- c("beta", "p", "gnsz.Mb", "sizerange.Mb")
df[, num] <- lapply(df[, num], as.numeric)
df <- df[df$p < 0.05, ]
plot(df$gnsz.Mb, df$sizerange.Mb)


library(beeswarm)
beeswarm(sizerange.Mb ~ cl, data = df)
# phylogenetic anova
library(phytools)
tree <- read.tree("../data/formatted_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
int <- intersect(df$species, tree$tip.label)
pruned.tree <- keep.tip(tree, int)
sig <- phylosig(pruned.tree,
                setNames(df$sizerange.Mb, df$species),
                method = "lambda",
                test = TRUE,
                nsim = 10000)[[4]]
if (sig < 0.05) {
  x <- setNames(df$cl, df$species)
  y <- setNames(df$sizerange.Mb, df$species)
  phylANOVA(pruned.tree, x, y, nsim = 10000, posthoc = TRUE)
}


# repeat expansion should be more useful for explaining gene density variation in 
# mammals while chromosome size variation should be more useful for explaining 
# gene density variation in reptiles






