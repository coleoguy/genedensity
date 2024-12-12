


# transform results
dat <- read.csv("../results/vertebrates/parsed.csv")
dat <- dat[!duplicated(dat$species), ]
dat <- dat[!is.na(dat$chromnum.1n), ]
dat$median.trans <- 1 - (dat$median/70)
dat$totalrep.prop <- dat$totalrep.pct * 0.01

# subset data
dat <- na.omit(dat[, c("species", "rsq", "clade", "median.trans", "totalrep.prop", "asmblysize.Mbp", "est.gnsz.Mbp")])
dat <- dat[dat$clade == "Mammalia", ]

# remove assembly with bloated assembly size
dat <- dat[dat$species != "Callithrix jacchus", ]

# assign weights
dat$w <- 1 - (abs(dat$asmblysize.Mbp - dat$est.gnsz.Mbp) / dat$est.gnsz.Mbp)

# visualize weights
hist(dat$w)

# raw data regression for proportion of repeats
plot(dat$totalrep.prop, dat$rsq)
abline(glm(dat$rsq ~ dat$totalrep.prop, weights = dat$w), col = "blue")
abline(glm(dat$rsq ~ dat$totalrep.prop))

# regression for bin with median repeat content
plot(dat$median.trans, dat$rsq)
abline(glm(dat$rsq ~ dat$median.trans, weights = dat$w), col = "blue")
abline(glm(dat$rsq ~ dat$median.trans))

# regression for the product of previous predictors
plot(dat$median.trans * dat$totalrep.prop, dat$rsq)
term <- dat$median.tran * dat$totalrep.prop
abline(glm(dat$rsq ~ term, weights = dat$w), col = "blue")
abline(glm(dat$rsq ~ term))




