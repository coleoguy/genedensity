
library(beeswarm)
taxo <- read.csv("../../results/rsq.csv")
repeats <- read.csv("../../results/repeat-results.csv")


dat <- merge(taxo, repeats, by.x = "species", by.y = "species", all = T)
dat$repeats <- rowSums(dat[, c("prop.dna", "prop.line", "prop.ltr", "prop.sine", "prop.others", "prop.unknown")])
dat <- dat[, c("species", "repeats", "clade", "class")]
dat <- na.omit(dat)
dat <- dat[dat$clade != "Others", ]
clades <- c("Mammalia", "Actinopterygii", "Sauropsida")
num.sp <- sapply(clades, function(cl) {
  nrow(dat[dat$clade == cl, ])})
map <- setNames(c("#d95f02", "#7570b3", "#1b9e77"), clades)
cols <- map[dat$clade]
cols[which(dat$class == "Reptilia")] <- "#e7298a"
beeswarm(repeats ~ clade, 
         xlab = NA, 
         ylab = NA, 
         ylim = c(0, 0.8), 
         data = dat,
         pwcol = cols,
         pch = 16,
         spacing = 1.4,
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
mtext("proportion of masked nucleotides", side = 2, line = 2.4)
# legend
points(2.9, 0.71, pch = 16, col = "#1b9e77")
points(2.9, 0.64, pch = 16, col = "#e7298a")
text(3, 0.71, 
     adj = c(0, 0.5), cex = 0.9, 
     labels = "Aves")
text(3, 0.64, 
     adj = c(0, 0.5), cex = 0.9, 
     labels = "Reptilia")
rect(2.75, 0.58, 3.44, 0.77, 
     border = "black", lwd = 1)

