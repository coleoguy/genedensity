
# before pgls: new parsing method kept 21 mammals and rsq~transformed median has slope -1.6988 and p value 0.016207; no need to exponentiate weights
# after pgls: 20 species, beta = 1.747489, p = 0.0183, r2 = 0.2731145, predicted rsq diff between highest and lowest medians: 0.37446190
# what about other clades?


# new model
library(viridis)
dat <- read.csv("../results/parsed.csv")
d <- dat
dat <- dat[!duplicated(dat$species), ]
dat <- dat[!is.na(dat$chromnum.1n), ]
dat$median.trans <- 1 - (dat$median/70)
dat$totalrep.prop <- dat$totalrep.pct * 0.01
dat <- na.omit(dat[, c("species", "rsq", "class", "order", "family", "clade", "median.trans", "totalrep.prop", "chromnum.1n", "est.gnsz.Mbp", "asmblysize.Mbp")])
dat <- dat[dat$clade == "Mammalia", ]
dat <- dat[dat$species != "Callithrix jacchus", ]
library(phytools)
tree <- read.tree("../data/formatted_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
int <- intersect(tree$tip.label, dat$species)
pruned.tree <- keep.tip(tree, int)
dat1 <- dat[dat$species %in% int, ]
dat1 <- dat1[match(pruned.tree$tip.label, dat1$species), ]
model <- glm(rsq ~ median.trans, data = dat1)
res <- setNames(resid(model), dat1$species)
phylosig(pruned.tree, res, method="lambda", test=TRUE)
library(nlme)
summary(gls(rsq ~ median.trans, 
            data = dat1))
library(piecewiseSEM)
rsquared(gls(rsq ~ median.trans, 
             data = dat1))














# old stuff

library(viridis)
# transform results
dat <- read.csv("../results/parsed.csv")
d <- dat
dat <- dat[!duplicated(dat$species), ]
dat <- dat[!is.na(dat$chromnum.1n), ]
dat$median.trans <- 1 - (dat$median/70)
dat$totalrep.prop <- dat$totalrep.pct * 0.01
contigs <- c()
for (h in dat$species) {
  sub <- d[d$species == h, ]
  contigs <- c(contigs, nrow(sub))
}
dat$contigs <- contigs

# subset data
dat <- na.omit(dat[, c("species", "rsq", "class", "order", "family", "clade", "median.trans", "totalrep.prop", "chromnum.1n", "est.gnsz.Mbp", "asmblysize.Mbp", "contigs")])
dat <- dat[dat$clade == "Mammalia", ]

# marsupials?
mars <- c()
for (g in 1:nrow(dat)) {
  if (dat[g, ]$order %in% c("Didelphimorphia", 
                            "Paucituberculata", 
                            "Microbiotheria",
                            "Dasyuromorphia", 
                            "Notoryctemorphia",
                            "Peramelemorphia", 
                            "Diprotodontia")){
    mars <- c(mars, TRUE)
    
  } else {
    mars <- c(mars, FALSE)
  }
}
dat$mars <- mars

# remove assembly with bloated assembly size
dat <- dat[dat$species != "Callithrix jacchus", ]

# assign weights; weights approach 1 as assembly sizes approach genome sizes. weights tend toward 0 as assembly sizes deviate from genome sizes
dat$w <- 1 - (abs(dat$asmblysize.Mbp - dat$est.gnsz.Mbp) / dat$est.gnsz.Mbp)

# marsupials
mars <- dat[, c("species", "mars", "rsq", "order")]
mars <- mars[mars$order != "Monotrema", ]
obs.diff <- mean(mars[mars$mars == TRUE, ]$rsq) - mean(mars[mars$mars == FALSE, ]$rsq)
null.diff <- c()
for (f in 1:10000) {
  null <- sample(mars)
  null.diff <- c(null.diff, mean(null[null$mars == TRUE, ]$rsq) - mean(null[null$mars == FALSE, ]$rsq))
}
p <- mean(abs(null.diff) >= abs(obs.diff))
print(paste0("p = ", p))
# nope




# median
cols <- viridis(length(unique(dat$w)), alpha = 0.45)[as.factor(dat$w)]
plot(dat$median.trans, 
     dat$rsq,
     col = cols,
     pch = 16)
abline(glm(dat$rsq ~ dat$median.trans, weights = dat$w), col = "blue")
abline(glm(dat$rsq ~ dat$median.trans))
legend("bottomright", 
       legend = round(seq(min(dat$median.trans), max(dat$median.trans), length.out = 5), 2), 
       fill = viridis(5), 
       title = "Legend")

# totsl repeat
plot(dat$totalrep.prop, 
     dat$rsq,
     col = cols,
     pch = 16)
abline(glm(dat$rsq ~ dat$totalrep.prop, weights = dat$w), col = "blue")
abline(glm(dat$rsq ~ dat$totalrep.prop))
legend("bottomright", 
       legend = round(seq(min(dat$totalrep.prop), max(dat$totalrep.prop), length.out = 5), 2), 
       fill = viridis(5), 
       title = "Legend")

# product of total repeat and median
plot(dat$median.trans * dat$totalrep.prop,  
     dat$rsq,
     col = cols,
     pch = 16)
term <- dat$median.trans * dat$totalrep.prop
abline(glm(dat$rsq ~ term, weights = dat$w), col = "blue")
abline(glm(dat$rsq ~ term))
legend("bottomright", 
       legend = round(seq(min(dat$median.trans * dat$totalrep.prop), max(dat$median.trans * dat$totalrep.prop), length.out = 5), 2), 
       fill = viridis(5), 
       title = "Legend")









# average number of complete genomes
sum(dat$w)

# add an exponent to emphasize high quality genomes
dat$w <- dat$w^5

sum(dat$w)
hist(dat$w)

# model
model <- glm(rsq ~ chromnum.1n + median.trans + median.trans:chromnum.1n, weights = dat$w, data = dat)

# generate model predictions and format results
x <- seq(min(dat$median.trans), max(dat$median.trans), length.out = 100)
y <- seq(min(dat$chromnum.1n), max(dat$chromnum.1n), length.out = 100)
grid <- expand.grid(median.trans = x, chromnum.1n = y)
grid$rsq <- predict(model, newdata = grid, type = "response")
z <- matrix(grid$rsq, nrow = length(x), ncol = length(y))

# plot
par(mar = c(4, 4, 3, 8))
image(x = x, 
      y = y, 
      z = z, 
      col = viridis(100), 
      xlab = "",
      ylab = "",
      main = "Title")
mtext("Expansion Recency", side=1, line=2.5)
mtext("Chromosome Number", side=2, line=2.5)
contour(x, 
        y, 
        z, 
        add = TRUE,
        col = "black", 
        lwd = 1,
        drawlabels = FALSE)
par(new = TRUE)
par(mar = c(4, 25, 3, 6))
z_range <- seq(min(z), max(z), length.out = 100)
image(1, 
      z_range, 
      t(matrix(z_range)), 
      col = viridis(100), 
      xaxt = "n",
      yaxt = "n", 
      xlab = "", 
      ylab = "")
axis(4, at = pretty(z_range), labels = round(pretty(z_range), 2))
mtext("Predicted Consistency", side=4, line=2.5)
par(mar = c(5.1, 4.1, 4.1, 2.1))





# group species by very low chromosome numbers or very complete genomes
wgroup <- chromgroup <- c()
for (i in 1:nrow(dat)) {
  if (dat[i, ]$w >= 0.95) {
    wgroup <- c(wgroup, T)
  } else {
    wgroup <- c(wgroup, F)
  }
  if (dat[i, ]$chromnum.1n <= 10) {
    chromgroup <- c(chromgroup, T)
  } else {
    chromgroup <- c(chromgroup, F)
  }
}

# tests
dat$wgroup <- wgroup
dat$chromgroup <- chromgroup
wilcox.test(w ~ chromgroup, data = dat)



library(phytools)
tree <- read.tree("../data/formatted_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
int <- intersect(tree$tip.label, dat$species)
pruned.tree <- keep.tip(tree, int)
dat1 <- dat[dat$species %in% int, ]
dat1 <- dat1[match(pruned.tree$tip.label, dat1$species), ]
model <- glm(rsq ~ median.trans, weights = w, data = dat1)
res <- setNames(resid(model), dat1$species)
phylosig(pruned.tree, res, method="lambda", test=TRUE)
library(nlme)
summary(gls(rsq ~ median.trans, 
             weights = varFixed(~w), 
             data = dat1))
library(piecewiseSEM)
rsquared(gls(rsq ~ median.trans, 
             weights = varFixed(~w), 
             data = dat1))
