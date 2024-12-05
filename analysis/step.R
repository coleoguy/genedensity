


dat <- read.csv("../results/vertebrates/parsed.csv")
dat <- dat[!duplicated(dat$species), ]
cols <- colnames(dat)

# add mean chromosome size to dat
meansize <- c()
for (i in unique(dat$species)) {
  sub <- subset(dat, species == i)
  meansize <- c(meansize, mean(sub$size.Mbp))
}
dat$meansize <- meansize



# all variables tested: 
# mediandvg -> bin containing median divergence
# beta -> beta coefficient of glm(gene number ~ contig size)
# chromnum.1n
# meansize -> mean chromosome size
# totalrep.pct -> repeat percent

################### stepwise selection for 2 variables ######################

# testing combinations
d <- na.omit(dat[, c("rsq", "clade", "mediandvg", "beta", "chromnum.1n", "meansize", "totalrep.pct")])
fit <- glm(d$rsq ~ d$beta * d$meansize)
step(fit) # -22.56
fit <- glm(d$rsq ~ d$chromnum.1n * d$meansize)
step(fit) # -21.65
fit <- glm(d$rsq ~ d$mediandvg * d$meansize)
step(fit) # -21.4
fit <- glm(d$rsq ~ d$meansize * d$totalrep.pct)
step(fit) # -20.83
fit <- glm(d$rsq ~ d$mediandvg * d$totalrep.pct)
step(fit) # -16.67
fit <- glm(d$rsq ~ d$mediandvg * d$beta)
step(fit) # -15.22
fit <- glm(d$rsq ~ d$mediandvg * d$chromnum.1n)
step(fit) # -15.22
fit <- glm(d$rsq ~ d$beta * d$chromnum.1n)
step(fit) # -13.6
fit <- glm(d$rsq ~ d$beta * d$totalrep.pct)
step(fit) # -13.6
fit <- glm(d$rsq ~ d$chromnum.1n * d$totalrep.pct)
step(fit) # -13.33

# lowest AIC
# glm(rsq ~ beta * meansize)
d <- na.omit(dat[, c("rsq", "clade", "beta", "meansize")])
summary(glm(d$rsq ~ d$beta * d$meansize))
m <- d[d$clade == "Mammalia", ]
summary(glm(m$rsq ~ m$beta * m$meansize))
f <- d[d$clade == "Actinopterygii", ]
summary(glm(d$rsq ~ d$beta * d$meansize))
r <- d[d$clade == "Sauria", ]
summary(glm(r$rsq ~ r$beta * r$meansize))

# second-lowest AIC
# significant in all clades
# glm(rsq ~ chromnum.1n * meansize)
d <- na.omit(dat[, c("rsq", "clade", "chromnum.1n", "meansize")])
summary(glm(d$rsq ~ d$chromnum.1n + d$meansize))
m <- d[d$clade == "Mammalia", ]
summary(glm(m$rsq ~ m$chromnum.1n + m$meansize))
f <- d[d$clade == "Actinopterygii", ]
summary(glm(d$rsq ~ d$chromnum.1n + d$meansize))
r <- d[d$clade == "Sauria", ]
summary(glm(r$rsq ~ r$chromnum.1n + r$meansize))

# third-lowest AIC
# glm(rsq ~ mediandvg + meansize)
d <- na.omit(dat[, c("rsq", "clade", "mediandvg", "meansize")])
summary(glm(d$rsq ~ d$mediandvg + d$meansize))
m <- d[d$clade == "Mammalia", ]
summary(glm(m$rsq ~ m$mediandvg + m$meansize))
f <- d[d$clade == "Actinopterygii", ]
summary(glm(f$rsq ~ f$mediandvg + f$meansize))
m <- d[d$clade == "Mammalia", ]
summary(glm(r$rsq ~ r$mediandvg + r$meansize))







#################### stepwise selection for 3 variables ######################

# testing combinations
d <- na.omit(dat[, c("rsq", "clade", "mediandvg", "beta", "chromnum.1n", "meansize", "totalrep.pct")])
fit <- glm(d$rsq ~ d$mediandvg * d$meansize * d$totalrep.pct)
step(fit) # -29.18
fit <- glm(d$rsq ~ d$mediandvg * d$chromnum.1n * d$meansize)
step(fit) # -27.92
fit <- glm(d$rsq ~ d$beta * d$chromnum.1n * d$meansize)
step(fit) # -26.6
fit <- glm(d$rsq ~ d$beta * d$meansize * d$totalrep.pct)
step(fit) # -24.99
fit <- glm(d$rsq ~ d$mediandvg * d$beta * d$meansize)
step(fit) # -22.56
fit <- glm(d$rsq ~ d$chromnum.1n * d$meansize * d$totalrep.pct)
step(fit) # -20.59
fit <- glm(d$rsq ~ d$mediandvg * d$beta * d$chromnum.1n)
step(fit) # -16.89
fit <- glm(d$rsq ~ d$mediandvg * d$chromnum.1n * d$totalrep.pct)
step(fit) # -16.67
fit <- glm(d$rsq ~ d$mediandvg * d$beta * d$totalrep.pct)
step(fit) # -16.67
fit <- glm(d$rsq ~ d$beta * d$chromnum.1n * d$totalrep.pct)
step(fit) # -13.6

# lowest AIC
# glm(rsq ~ mediandvg + meansize + totalrep.pct + mediandvg:totalrep.pct)
d <- na.omit(dat[, c("rsq", "clade", "mediandvg", "meansize", "totalrep.pct")])
summary(glm(d$rsq ~ d$mediandvg + d$meansize + d$totalrep.pct + d$mediandvg:d$totalrep.pct))
m <- d[d$clade == "Mammalia", ]
summary(glm(m$rsq ~ m$mediandvg + m$meansize + m$totalrep.pct + m$mediandvg:m$totalrep.pct))
f <- d[d$clade == "Actinopterygii", ]
summary(glm(f$rsq ~ f$mediandvg + f$meansize + f$totalrep.pct + f$mediandvg:f$totalrep.pct))
r <- d[d$clade == "Sauria", ]
summary(glm(r$rsq ~ r$mediandvg + r$meansize + r$totalrep.pct + r$mediandvg:r$totalrep.pct))

# second-lowest AIC
# significant for mammals
# glm(rsq ~ mediandvg + chromnum.1n + meansize + mediandvg:chromnum.1n + mediandvg:meansize)  
d <- na.omit(dat[, c("rsq", "clade", "mediandvg", "chromnum.1n", "meansize")])
summary(glm(d$rsq ~ d$mediandvg + d$chromnum.1n + d$meansize + d$mediandvg:d$chromnum.1n + d$mediandvg:d$meansize))
m <- d[d$clade == "Mammalia", ]
summary(glm(m$rsq ~ m$mediandvg + m$chromnum.1n + m$meansize + m$mediandvg:m$chromnum.1n + m$mediandvg:m$meansize))
f <- d[d$clade == "Actinopterygii", ]
summary(glm(f$rsq ~ f$mediandvg + f$chromnum.1n + f$meansize + f$mediandvg:f$chromnum.1n + f$mediandvg:f$meansize))
r <- d[d$clade == "Sauria", ]
summary(glm(r$rsq ~ r$mediandvg + r$chromnum.1n + r$meansize + r$mediandvg:r$chromnum.1n + r$mediandvg:r$meansize))

# third-lowest AIC
# glm(rsq ~ beta + chromnum.1n + meansize + beta:meansize)
d <- na.omit(dat[, c("rsq", "clade", "beta", "chromnum.1n", "meansize")])
summary(glm(d$rsq ~ d$beta + d$chromnum.1n + d$meansize + d$beta:d$meansize))
m <- d[d$clade == "Mammalia", ]
summary(glm(m$rsq ~ m$beta + m$chromnum.1n + m$meansize + m$beta:m$meansize))
f <- d[d$clade == "Actinopterygii", ]
summary(glm(f$rsq ~ f$beta + f$chromnum.1n + f$meansize + f$beta:f$meansize))
r <- d[d$clade == "Sauria", ]
summary(glm(r$rsq ~ r$beta + r$chromnum.1n + r$meansize + r$beta:r$meansize))



