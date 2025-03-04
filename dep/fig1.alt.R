


# rsq scatter
dat <- read.csv("../results/parsed.csv")
dat <- dat[dat$thrs == 0.8, ]
dat <- dat[dat$clade %in% c("Mammalia", "Sauria", "Actinopterygii"), ]
dat <- dat[!is.na(dat$chromnum.1n), ]
sub <- dat[dat$species %in% c("Bos taurus", "Peromyscus maniculatus bairdii"), ]
sub$genecount <- sub$genecount / 1000
cols <- c("#f4a582", "#0571b0")[as.factor(sub$species)]
par(mar = c(6, 4, 2, 2)+0.1)
plot(x = sub$size.Mb,
     y = sub$genecount,
     xlab = NA,
     ylab = "Gene count (thousands)",
     ylim = c(0.3, 3), 
     xlim = c(30, 200), 
     col = cols,
     cex = 1.5, 
     pch = 16)
text(91, -0.5, "Chromosome size (Mb)", xpd = NA, adj = 0)
x1 <- sub[sub$species == "Bos taurus", ]$size.Mb
y1 <- sub[sub$species == "Bos taurus", ]$genecount
model1 <- lm(y1 ~ x1)
abline(model1,
col = "#f4a582",
lwd = 1.5)
x2 <- sub[sub$species == "Peromyscus maniculatus bairdii", ]$size.Mb
y2 <- sub[sub$species == "Peromyscus maniculatus bairdii", ]$genecount
model2 <- lm(y2 ~ x2)
abline(model2,
col = "#0571b0",
lwd = 1.5)
par(mar = c(5.1, 4.1, 4.1, 2.1))

