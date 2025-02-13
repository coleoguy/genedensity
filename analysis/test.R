

library(data.table)
dat <- fread("../results/parsed.csv")
dat <- dat[dat$thrs == 0.8, ]
dat <- dat[!is.na(dat$chromnum.1n), ]
dat <- dat[!duplicated(dat$species), ]
dat <- dat[dat$clade == "Actinopterygii", ]
dat <- na.omit(dat[, c("species","clade", "rsq", "sine.rep.pct")])
dat$sine.rep.prop <- dat$sine.rep.pct / 100
plot(dat$sine.rep.prop, 
     dat$rsq, 
     xlab = "repeat proportion", 
     ylab = "r2", 
     pch = 16, 
     cex = 0.5, 
     main = "fish SINE proportion", 
     xlim = c(0, 1), 
     ylim = c(0, 1))
abline(glm(rsq ~ sine.rep.prop, data = dat))


















library(data.table)
dat <- fread("../results/parsed.csv")
dat <- dat[dat$thrs == 0.97, ]
dat <- dat[!is.na(dat$chromnum.1n), ]
dat <- dat[!duplicated(dat$species), ]
dat <- dat[dat$clade == "Mammalia", ]
dat <- na.omit(dat[, c("species", "clade", "rsq", "ltr.rep.pct", "ltr.rep.median")])
  
median <- 1 - (dat$ltr.rep.median/70)
median.range <- range(na.omit(median))
dat$median.trans <- (median - median.range[1]) / diff(median.range)
vol <- dat$ltr.rep.pct / 100
vol.range <- range(na.omit(vol))
dat$rep.prop <- (vol - vol.range[1]) / diff(vol.range)
dat <- na.omit(dat[, c("species", "clade", "rsq", "median.trans", "rep.prop")])
dat <- dat[order(dat$median.trans), ]

nsplit <- 4
base <- floor(nrow(dat)/nsplit)
nrow(dat) - 4*base


plot(dat$ltr.rep.prop, 
     dat$rsq, 
     xlab = "repeat proportion", 
     ylab = "r2", 
     pch = 16, 
     cex = 0.5, 
     main = "fish ltr proportion", 
     xlim = c(0, 1), 
     ylim = c(0, 1))
abline(glm(rsq ~ ltr.rep.prop, data = dat))

