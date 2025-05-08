
dat <- read.csv("../results/maskedsz.csv", 
                colClasses = c("date" = "character"))
dat$lower <- dat$lower / 1000000
dat$N <- dat$lower / 1000000
dat$assembly_size <- dat$assembly_size / 1000000
dat$date <- as.Date(dat$date, format = "%m%d%Y")

op <- par(ask = TRUE)
for (i in unique(dat$sp)) {
  sub <- dat[dat$sp == i, ]
  plot(sub$date, sub$lower, main = i, xlab = "date", ylab = "size (Mb)")
  if (nrow(sub) > 2) {
    fit <- glm(sub$lower ~ sub$date)
    abline(fit$coefficients[1], fit$coefficients[2])
  }
}
par(op)
