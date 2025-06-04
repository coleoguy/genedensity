


dat <- read.csv("../../results/major-assem.csv")
dat$date <- as.Date(dat$date, format = "%d-%b-%Y")

for (i in unique(dat$sp)) {
  da <- dat[dat$sp == i, ]
  y <- da$lowercase / da$assembly_size
  plot(da$date, y, 
       xlab = "date", ylab = "y", 
       main = gsub("_", " ", i))
  readline("Press [Enter] to continue...")
}


dat <- read.csv("../../results/major-assem.csv")
dat$date <- as.Date(dat$date, format = "%d-%b-%Y")
dat <- dat[dat$sp == "Mus_musculus", ]
dat <- dat[dat$assem != "GRCm39", ]
y <- dat$lowercase / dat$assembly_size
plot(dat$date, y, 
     pch = 16, cex = 1, 
     xlim = c(11300, 15600), 
     ylim = c(0.25, 0.47), 
     xlab = "Release date", ylab = "Proportion of masked nucleotides", 
     main = "Mus musculus")
fit <- glm(y ~ date, data = dat)
abline(fit, lwd = 1.5)

