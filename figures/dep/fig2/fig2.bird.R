
par(mar = c(5, 4, 4, 5) + 0.1, 
    xpd = NA)

library(beeswarm)
dat <- read.csv("../../results/rsq.csv")

dat <- dat[, c("species", "rsq", "clade", "class")]
dat <- na.omit(dat)
dat <- dat[dat$clade != "Others", ]
clades <- c("Mammalia", "Actinopterygii", "Sauropsida")
num.sp <- sapply(clades, function(cl) {
  nrow(dat[dat$clade == cl, ])})
map <- setNames(c("#d95f02", "#7570b3", "#1b9e77"), clades)
cols <- map[dat$clade]
cols[which(dat$class == "Reptilia")] <- "#e7298a"
beeswarm(rsq ~ clade, 
         xlab = NA, 
         ylab = NA, 
         ylim = c(0.16, 1.03), 
         data = dat,
         pwcol = cols,
         pch = 16,
         spacing = 1,
         xaxt = "n", yaxt = "n")

# xaxt
axis(side   = 1,
     at     = seq(num.sp),
     labels = c("Mammals", "Ray-finned fish", "Reptiles"), 
     mgp = c(3, 0.9, 0))
axis(side   = 1,
     at     = seq(num.sp),
     labels = paste0("n=", num.sp), 
     mgp = c(3, 2.1, 0))
# yaxt
axis(side = 2, mgp = c(3, 0.9, 0))
mtext("R-squared", side = 2, line = 2.4)



# legend
points(3.8, 0.82, pch = 16, col = "#1b9e77")
points(3.8, 0.72, pch = 16, col = "#e7298a")
text(3.9, 0.82, 
     adj = c(0, 0.5), cex = 0.9, 
     labels = "Aves")
text(3.9, 0.72, 
     adj = c(0, 0.5), cex = 0.9, 
     labels = "Reptilia")
#rect(2.75, 0.18, 3.44, 0.37, 
#     border = "black", lwd = 1)
