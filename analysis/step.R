


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



# variables tested: 
# beta -> beta coefficient of glm(gene number ~ contig size)
# chromnum.1n
# est.gnsz.Mbp -> estimated geome size
# meansize -> mean chromosome size
# mediandvg -> bin containing median divergence
# totalrep.pct -> repeat percent

################### stepwise model selection #####################

d <- na.omit(dat[, c("rsq", "clade", "beta", "chromnum.1n", "est.gnsz.Mbp", "meansize", "mediandvg", "totalrep.pct")])

# testing 2-variable combinations
fit <- glm(d$rsq ~ d$est.gnsz.Mbp * d$meansize)
step(fit) # -31.86
fit <- glm(d$rsq ~ d$beta * d$est.gnsz.Mbp)
step(fit) # -20.26
fit <- glm(d$rsq ~ d$beta * d$meansize)
step(fit) # -19.75
fit <- glm(d$rsq ~ d$chromnum.1n * d$meansize)
step(fit) # -19.73
fit <- glm(d$rsq ~ d$meansize * d$mediandvg)
step(fit) # -19.59
fit <- glm(d$rsq ~ d$est.gnsz.Mbp * d$mediandvg)
step(fit) # -19.14
fit <- glm(d$rsq ~ d$chromnum.1n * d$est.gnsz.Mbp)
step(fit) # -18.85
fit <- glm(d$rsq ~ d$est.gnsz.Mbp * d$totalrep.pct)
step(fit) # -18.85
fit <- glm(d$rsq ~ d$meansize * d$totalrep.pct) 
step(fit) # -18.51
fit <- glm(d$rsq ~ d$mediandvg * d$totalrep.pct)  
step(fit) # -14.64
fit <- glm(d$rsq ~ d$beta * d$mediandvg) 
step(fit) # -13.27
fit <- glm(d$rsq ~ d$chromnum.1n * d$mediandvg)
step(fit) # -13.27
fit <- glm(d$rsq ~ d$beta * d$chromnum.1n)
step(fit) # -11.55
fit <- glm(d$rsq ~ d$beta * d$totalrep.pct)
step(fit) # -11.28
fit <- glm(d$rsq ~ d$chromnum.1n * d$totalrep.pct)
step(fit) # -11.18


# testing 3-variable combinations
fit <- glm(d$rsq ~ d$beta * d$est.gnsz.Mbp * d$meansize)
step(fit) # -43.22
fit <- glm(d$rsq ~ d$beta * d$est.gnsz.Mbp * d$mediandvg)
step(fit) # -37.46
fit <- glm(d$rsq ~ d$beta * d$est.gnsz.Mbp * d$totalrep.pct)
step(fit) # -36.21
fit <- glm(d$rsq ~ d$chromnum.1n * d$est.gnsz.Mbp * d$meansize)
step(fit) # -34.18
fit <- glm(d$rsq ~ d$est.gnsz.Mbp * d$meansize * d$totalrep.pct)
step(fit) # -31.86
fit <- glm(d$rsq ~ d$est.gnsz.Mbp * d$meansize * d$mediandvg)
step(fit) # -31.39
fit <- glm(d$rsq ~ d$beta * d$chromnum.1n * d$est.gnsz.Mbp) 
step(fit) # -29.26
fit <- glm(d$rsq ~ d$meansize * d$mediandvg * d$totalrep.pct)
step(fit) # -27.11
fit <- glm(d$rsq ~ d$chromnum.1n * d$meansize * d$mediandvg)
step(fit) # -25.43
fit <- glm(d$rsq ~ d$beta * d$chromnum.1n * d$meansize)  
step(fit) # -23.62
fit <- glm(d$rsq ~ d$beta * d$meansize * d$totalrep.pct)
step(fit) # -21.55
fit <- glm(d$rsq ~ d$beta * d$meansize * d$mediandvg)
step(fit) # -19.75
fit <- glm(d$rsq ~ d$chromnum.1n * d$meansize * d$totalrep.pct)
step(fit) # -19.2
fit <- glm(d$rsq ~ d$chromnum.1n * d$est.gnsz.Mbp * d$mediandvg)
step(fit) # -19.14
fit <- glm(d$rsq ~ d$chromnum.1n * d$est.gnsz.Mbp * d$totalrep.pct)
step(fit) # -18.85
fit <- glm(d$rsq ~ d$est.gnsz.Mbp * d$mediandvg * d$totalrep.pct)
step(fit) # -18.11
fit <- glm(d$rsq ~ d$beta * d$chromnum.1n * d$mediandvg)
step(fit) # -14.74
fit <- glm(d$rsq ~ d$beta * d$mediandvg * d$totalrep.pct)
step(fit) # -14.64
fit <- glm(d$rsq ~ d$chromnum.1n * d$mediandvg * d$totalrep.pct)
step(fit) # -14.64
fit <- glm(d$rsq ~ d$beta * d$chromnum.1n * d$totalrep.pct)
step(fit) # -11.55




############# model summaries ##############


# glm(formula = rsq ~ est.gnsz.Mbp + meansize)
# lowest AIC for two-variable models
d <- na.omit(dat[, c("rsq", "clade", "est.gnsz.Mbp", "meansize")])
summary(glm(d$rsq ~ d$est.gnsz.Mbp * d$meansize))
m <- d[d$clade == "Mammalia", ]
summary(glm(m$rsq ~ m$est.gnsz.Mbp * m$meansize))
f <- d[d$clade == "Actinopterygii", ]
summary(glm(f$rsq ~ f$est.gnsz.Mbp * f$meansize))
r <- d[d$clade == "Sauria", ]
summary(glm(r$rsq ~ r$est.gnsz.Mbp * r$meansize))

# glm(formula = rsq ~ beta * est.gnsz.Mbp)
# second lowest AIC for two-variable models
# highly significant overall but not in each clade
d <- na.omit(dat[, c("rsq", "clade", "beta", "est.gnsz.Mbp")])
summary(glm(d$rsq ~ d$beta * d$est.gnsz.Mbp))
m <- d[d$clade == "Mammalia", ]
summary(glm(m$rsq ~ m$beta * m$est.gnsz.Mbp))
f <- d[d$clade == "Actinopterygii", ]
summary(glm(f$rsq ~ f$beta * f$est.gnsz.Mbp))
r <- d[d$clade == "Sauria", ]
summary(glm(r$rsq ~ r$beta * r$est.gnsz.Mbp))

# glm(formula = rsq ~ beta * meansize)
# third lowest AIC for two-variable models
d <- na.omit(dat[, c("rsq", "clade", "beta", "meansize")])
summary(glm(d$rsq ~ d$beta * d$meansize))
m <- d[d$clade == "Mammalia", ]
summary(glm(m$rsq ~ m$beta * m$meansize))
f <- d[d$clade == "Actinopterygii", ]
summary(glm(f$rsq ~ f$beta * f$meansize))
r <- d[d$clade == "Sauria", ]
summary(glm(r$rsq ~ r$beta * r$meansize))



# glm(formula = rsq ~ beta + est.gnsz.Mbp + meansize + beta:est.gnsz.Mbp)
# lowest AIC for three-variable models
# highly significant overall but not in each clade
d <- na.omit(dat[, c("rsq", "clade", "beta", "est.gnsz.Mbp", "meansize")])
summary(glm(d$rsq ~ d$beta + d$est.gnsz.Mbp + d$meansize + d$beta:d$est.gnsz.Mbp))
m <- d[d$clade == "Mammalia", ]
summary(glm(m$rsq ~ m$beta + m$est.gnsz.Mbp + m$meansize + m$beta:m$est.gnsz.Mbp))
f <- d[d$clade == "Actinopterygii", ]
summary(glm(f$rsq ~ f$beta + f$est.gnsz.Mbp + f$meansize + f$beta:f$est.gnsz.Mbp))
r <- d[d$clade == "Sauria", ]
summary(glm(r$rsq ~ r$beta + r$est.gnsz.Mbp + r$meansize + r$beta:r$est.gnsz.Mbp))


# glm(formula = rsq ~ beta + est.gnsz.Mbp + mediandvg + beta:est.gnsz.Mbp)
# second lowest AIC for three-variable models
# significant overall but not in each clade
d <- na.omit(dat[, c("rsq", "clade", "beta", "est.gnsz.Mbp", "mediandvg")])
summary(glm(d$rsq ~ d$beta * d$est.gnsz.Mbp + d$mediandvg))
m <- d[d$clade == "Mammalia", ]
summary(glm(m$rsq ~ m$beta * m$est.gnsz.Mbp + m$mediandvg))
f <- d[d$clade == "Actinopterygii", ]
summary(glm(f$rsq ~ f$beta * f$est.gnsz.Mbp + f$mediandvg))
r <- d[d$clade == "Sauria", ]
summary(glm(r$rsq ~ r$beta * r$est.gnsz.Mbp + r$mediandvg))

# glm(formula = rsq ~ beta + est.gnsz.Mbp + totalrep.pct + beta:est.gnsz.Mbp + beta:totalrep.pct)
# third lowest AIC for three- variable models
# significant overall but not in each clade
d <- na.omit(dat[, c("rsq", "clade", "beta", "est.gnsz.Mbp", "totalrep.pct")])
summary(glm(d$rsq ~ d$beta + d$est.gnsz.Mbp + d$totalrep.pct + d$beta:d$est.gnsz.Mbp + d$beta:d$totalrep.pct))
m <- d[d$clade == "Mammalia", ]
summary(glm(m$rsq ~ m$beta + m$est.gnsz.Mbp + m$totalrep.pct + m$beta:m$est.gnsz.Mbp + m$beta:m$totalrep.pct))
f <- d[d$clade == "Actinopterygii", ]
summary(glm(f$rsq ~ f$beta + f$est.gnsz.Mbp + f$totalrep.pct + f$beta:f$est.gnsz.Mbp + f$beta:f$totalrep.pct))
r <- d[d$clade == "Sauria", ]
summary(glm(r$rsq ~ r$beta + r$est.gnsz.Mbp + r$totalrep.pct + r$beta:r$est.gnsz.Mbp + r$beta:r$totalrep.pct))

