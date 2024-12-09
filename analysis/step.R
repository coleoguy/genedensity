

dat <- read.csv("../results/vertebrates/parsed.csv")
dat <- dat[!duplicated(dat$species), ]

d <- na.omit(dat[, c("rsq", "clade", "median", "totalrep.pct")])
summary(step(glm(d$rsq ~ d$median : d$totalrep.pct)))
int <- d$median * d$totalrep.pct
plot(d$rsq ~ int)
abs(diff(range(int))*summary(glm(d$rsq ~ d$median : d$totalrep.pct))$coefficients[2, 1]/diff(range(d$rsq)))

m <- d[d$clade == "Mammalia", ]
summary(glm(m$rsq ~ m$median : m$totalrep.pct))
int <- m$median * m$totalrep.pct
plot(m$rsq ~ int)
abs(diff(range(int))*summary(glm(m$rsq ~ m$median : m$totalrep.pct))$coefficients[2, 1]/diff(range(m$rsq)))

f <- d[d$clade == "Actinopterygii", ]
summary(glm(f$rsq ~ f$median : f$totalrep.pct))
int <- f$median * f$totalrep.pct
plot(f$rsq ~ int)
abs(diff(range(int))*summary(glm(f$rsq ~ f$median : f$totalrep.pct))$coefficients[2, 1]/diff(range(f$rsq)))

r <- d[d$clade == "Sauria", ]
summary(glm(r$rsq ~ r$median : r$totalrep.pct))
int <- r$median * r$totalrep.pct
plot(r$rsq ~ int)
abs(diff(range(int))*summary(glm(r$rsq ~ r$median : r$totalrep.pct))$coefficients[2, 1]/diff(range(r$rsq)))



















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


