
dat <- read.csv("../results/parsed.csv")
dat <- dat[dat$thrs == 0, ]
dens <- c()
for (i in  unique(dat$species)) {
  sub <- dat[dat$species == i, ]
  dens <- c(dens, setNames(mean(sub$genedens), unique(sub$clade)))
}
for (j in unique(names(dens))) {
  x <- mean(dens[which(names(dens) == j)])
  print(j)
  print(x)
}


# what happens to gene density variation if we break some chromosomes?


# zero gene density variation
size <- c(900, 800, 700, 600, 500, 400, 300, 200)
genes <- c(90, 80, 70, 60, 50, 40, 30, 20)
df <- data.frame(size, genes)
df$density <- df$genes/df$size
summary(lm(genes ~ size, data = df))$r.squared
# break the smallest chromosome in the middle, allocate genes nearly equally
size <- c(900, 800, 700, 600, 500, 400, 300, 100, 100)
genes <- c(90, 80, 70, 60, 50, 40, 30, 9, 11)
df <- data.frame(size, genes)
df$density <- df$genes/df$size
summary(lm(genes ~ size, data = df))$r.squared
# break the smallest chromosome in the middle, allocate genes unequally
size <- c(900, 800, 700, 600, 500, 400, 300, 100, 100)
genes <- c(90, 80, 70, 60, 50, 40, 30, 1, 11)
df <- data.frame(size, genes)
df$density <- df$genes/df$size
summary(lm(genes ~ size, data = df))$r.squared
# break the smallest chromosome at the tips, allocate genes equally
size <- c(900, 800, 700, 600, 500, 400, 300, 20, 180)
genes <- c(90, 80, 70, 60, 50, 40, 30, 2, 18)
df <- data.frame(size, genes)
df$density <- df$genes/df$size
summary(lm(genes ~ size, data = df))$r.squared
# break the smallest chromosome at the tips, give smallest fragnemt no genes
size <- c(900, 800, 700, 600, 500, 400, 300, 20, 180)
genes <- c(90, 80, 70, 60, 50, 40, 30, 0, 20)
df <- data.frame(size, genes)
df$density <- df$genes/df$size
summary(lm(genes ~ size, data = df))$r.squared
# break the smallest chromosome at the tips, give smallest fragnemt every gene
size <- c(900, 800, 700, 600, 500, 400, 300, 20, 180)
genes <- c(90, 80, 70, 60, 50, 40, 30, 20, 0)
df <- data.frame(size, genes)
df$density <- df$genes/df$size
summary(lm(genes ~ size, data = df))$r.squared
# break the largest chromosome in the middle, allocate genes nearly equally
size <- c(450, 450, 800, 700, 600, 500, 400, 300, 200)
genes <- c(44, 46, 80, 70, 60, 50, 40, 30, 20)
df <- data.frame(size, genes)
df$density <- df$genes/df$size
summary(lm(genes ~ size, data = df))$r.squared
# break the largest chromosome in the middle, allocate genes unequally
size <- c(450, 450, 800, 700, 600, 500, 400, 300, 200)
genes <- c(1, 89, 80, 70, 60, 50, 40, 30, 20)
df <- data.frame(size, genes)
df$density <- df$genes/df$size
summary(lm(genes ~ size, data = df))$r.squared
# break the largest chromosome at the tip, allocate genes equally
size <- c(100, 800, 800, 700, 600, 500, 400, 300, 200)
genes <- c(10, 80, 80, 70, 60, 50, 40, 30, 20)
df <- data.frame(size, genes)
df$density <- df$genes/df$size
summary(lm(genes ~ size, data = df))$r.squared
# break the largest chromosome at the tip, give larger fragment more genes
size <- c(100, 800, 800, 700, 600, 500, 400, 300, 200)
genes <- c(5, 85, 80, 70, 60, 50, 40, 30, 20)
df <- data.frame(size, genes)
df$density <- df$genes/df$size
summary(lm(genes ~ size, data = df))$r.squared
# break the largest chromosome at the tip, give smaller fragment more genes
size <- c(100, 800, 800, 700, 600, 500, 400, 300, 200)
genes <- c(80, 10, 80, 70, 60, 50, 40, 30, 20)
df <- data.frame(size, genes)
df$density <- df$genes/df$size
summary(lm(genes ~ size, data = df))$r.squared
# break the largest chromosome at the tip, give smaller fragment every gene
size <- c(100, 800, 800, 700, 600, 500, 400, 300, 200)
genes <- c(90, 0, 80, 70, 60, 50, 40, 30, 20)
df <- data.frame(size, genes)
df$density <- df$genes/df$size
summary(lm(genes ~ size, data = df))$r.squared
# break the largest chromosome at the tip, give smaller fragment no genes
size <- c(100, 800, 800, 700, 600, 500, 400, 300, 200)
genes <- c(0, 90, 80, 70, 60, 50, 40, 30, 20)
df <- data.frame(size, genes)
df$density <- df$genes/df$size
summary(lm(genes ~ size, data = df))$r.squared
# break the largest chromosome at the tip, give smaller fragment more genes
size <- c(100, 800, 800, 700, 600, 500, 400, 300, 200)
genes <- c(50, 40, 80, 70, 60, 50, 40, 30, 20)
df <- data.frame(size, genes)
df$density <- df$genes/df$size
summary(lm(genes ~ size, data = df))$r.squared




