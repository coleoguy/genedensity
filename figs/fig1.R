


# rsq scatter
dat <- read.csv("../results/parsed.csv")
dat <- dat[dat$thrs == 0.8, ]
dat <- dat[dat$clade %in% c("Mammalia", "Sauria", "Actinopterygii"), ]
dat <- dat[!is.na(dat$chromnum.1n), ]
dat$genecount <- dat$genecount / 1000
horse <- dat[dat$species == "Equus caballus", ]
meerkat <- dat[dat$species == "Suricata suricatta", ]
par(mfrow = c(1, 2))
par(mar = c(5.1, 4.1, 4.1, 0.9))
color <- adjustcolor("#405070", alpha.f = 0.7)
plot(horse$size.Mb, 
     horse$genecount, 
     pch = 16, 
     cex = 0.8, 
     col = color, 
     xlab = NA, 
     ylab = "Gene count (thousands)")
abline(glm(genecount ~ size.Mb, data = horse), lwd = 1.5, col = color)
par(mar = c(5.1, 2.9, 4.1, 2.1))
color <- adjustcolor("#405070", alpha.f = 0.7)
plot(meerkat$size.Mb, 
     meerkat$genecount, 
     pch = 16, 
     cex = 0.8, 
     col = color, 
     xlab = NA, 
     ylab = NA)
abline(glm(genecount ~ size.Mb, data = meerkat), lwd = 1.5, col = color)

text(-40, -0.01, "Chromosome size (Mb)", xpd = NA, adj = 0)
par(mfrow = c(1, 1))
par(mar = c(5.1, 4.1, 4.1, 2.1))
