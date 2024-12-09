

dat <- read.csv("../results/vertebrates/parsed.csv")
dat <- dat[!duplicated(dat$species), ]

d <- na.omit(dat[, c("rsq", "clade", "median", "totalrep.pct")])
d$median <- 70 - d$median # subtract median by max divergence
summary(step(glm(d$rsq ~ d$median * d$totalrep.pct)))

m <- d[d$clade == "Mammalia", ]
summary(step(glm(m$rsq ~ m$median * m$totalrep.pct)))

f <- d[d$clade == "Actinopterygii", ]
summary(step(glm(f$rsq ~ f$median * f$totalrep.pct)))

r <- d[d$clade == "Sauria", ]
summary(step(glm(r$rsq ~ r$median * r$totalrep.pct)))

















permTest <- function(x, y, reps, method) {
  perm.y <- replicate(reps, sample(y))
  perm.cor <- apply(perm.y, MARGIN = 2, function(col) cor(x, col, method = method))
  pval <- mean(abs(perm.cor) >= abs(cor(x, y, method = method)))
  return(pval)
}

y <- d$rsq
x <- d$median * d$totalrep.pct
r <- signif(cor(x, y, method = "pearson"), 3)
perm.pval <- signif(permTest(x, y, 100000, "pearson"), 3)






