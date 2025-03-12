

library(MuMIn) 
library(viridis)
models <- readRDS("../results/all.models.rds")
models <- models[1:length(which(cumsum(models$weight) <= 0.95))] # confidence set
avg <- model.avg(models) # average

terms <- c("age.dna", "prop.dna", "age.dna:prop.dna", 
           "age.line", "prop.line", "age.line:prop.line", 
           "age.ltr", "prop.ltr", "age.ltr:prop.ltr", 
           "age.sine", "prop.sine", "age.sine:prop.sine", 
           "age.unknown", "prop.unknown", "age.unknown:prop.unknown", 
           "age.others", "prop.others", "age.others:prop.others"
           )
est <- coef(avg)[match(terms, names(coef(avg)))] # parameter estimates
int <- confint(avg, full = FALSE)[terms, ] # confidence intervals
imp <- sw(models)[1:length(terms)] # importance
imp <- imp[match(terms, names(imp))] # reorder importance
x <- 1:length(terms) # x positions

# color mapping
res <- 1000 # resolution
palette <- magma(res) # palette
cols <- palette[round(((imp - min(imp)) / diff(range(imp))) * (res-1)) + 1] # colors

# x labels
labels <- c("DNA age", "DNA proportion", "DNA interaction", 
           "LINE age", "LINE proportion", "LINE interaction", 
           "LTR age", "LTR proportion", "LTR interaction", 
           "SINE age", "SINE proportion", "SINE interaction", 
           "Unidtf age", "Unidtf proportion", "Unidtf interaction", 
           "Other age", "Other proportion", "Other interaction" 
)

par(oma = c(0, 0, 3, 0))
layout(matrix(1:2, ncol = 2), widths = c(4, 1)) # make 2 plots

# main plot
par(mar = c(8, 4, 1, 0))
plot(y = est, x = x, type = "n", ylim = 1.05 * range(int), 
     ylab = "Parameter estimate", xlab = NA, axes = FALSE,
     xlim = c(0.75, length(terms) + 0.2)) # plot
abline(h = 0, lty = 2, col = "gray") # line at y = 0
segments(x, int[, 1], x, int[, 2], lwd = 2) # confidence bars
points(x, est, pch = 16, col = cols) # colored points
axis(2) # y axis
axis(1, at = x, labels = labels, las = 2) # x axis
box()

# color bar
par(mar = c(8, 1, 1, 4))
height <- seq(min(imp), max(imp), length.out = res + 1) # y values
z <- matrix(seq(min(imp), max(imp), length.out = res), nrow = 1, ncol = res) # color gradient
image(x = c(0, 1), y = height, z = z, col = palette, 
      axes = FALSE, xlab = "", ylab = "") # make color bar
ticks <- seq(min(imp), max(imp), length.out = 5) # ticks
axis(4, at = ticks, labels = round(ticks, 2), las = 1) # y axis

# title
mtext("Parameter estimates for averaged model of all species", 
      outer = TRUE, cex = 1.1, line = 0, font = 2, family = "sans", 
      adj = 0.35)
par(mar = c(5, 4, 4, 2) + 0.1) 





























library(MuMIn) 
library(viridis)
models <- readRDS("../results/Mammalia.models.rds")
models <- models[1:length(which(cumsum(models$weight) <= 0.95))] # confidence set
avg <- model.avg(models) # average

terms <- c("age.dna", "prop.dna", "age.dna:prop.dna", 
           "age.line", "prop.line", "age.line:prop.line", 
           "age.ltr", "prop.ltr", "age.ltr:prop.ltr", 
           "age.sine", "prop.sine", "age.sine:prop.sine", 
           "age.unknown", "prop.unknown", "age.unknown:prop.unknown", 
           "age.others", "prop.others", "age.others:prop.others"
)
est <- coef(avg)[match(terms, names(coef(avg)))] # parameter estimates
int <- confint(avg, full = FALSE)[terms, ] # confidence intervals
imp <- sw(models)[1:length(terms)] # importance
imp <- imp[match(terms, names(imp))] # reorder importance
x <- 1:length(terms) # x positions

# color mapping
res <- 1000 # resolution
palette <- magma(res) # palette
cols <- palette[round(((imp - min(imp)) / diff(range(imp))) * (res-1)) + 1] # colors

# x labels
labels <- c("DNA age", "DNA proportion", "DNA interaction", 
            "LINE age", "LINE proportion", "LINE interaction", 
            "LTR age", "LTR proportion", "LTR interaction", 
            "SINE age", "SINE proportion", "SINE interaction", 
            "Unidtf age", "Unidtf proportion", "Unidtf interaction", 
            "Other age", "Other proportion", "Other interaction" 
)

par(oma = c(0, 0, 3, 0))
layout(matrix(1:2, ncol = 2), widths = c(4, 1)) # make 2 plots

# main plot
par(mar = c(8, 4, 1, 0))
plot(y = est, x = x, type = "n", ylim = 1.05 * range(int), 
     ylab = "Parameter estimate", xlab = NA, axes = FALSE,
     xlim = c(0.75, length(terms) + 0.2)) # plot
abline(h = 0, lty = 2, col = "gray") # line at y = 0
segments(x, int[, 1], x, int[, 2], lwd = 2) # confidence bars
points(x, est, pch = 16, col = cols) # colored points
axis(2) # y axis
axis(1, at = x, labels = labels, las = 2) # x axis
box()

# color bar
par(mar = c(8, 1, 1, 4))
height <- seq(min(imp), max(imp), length.out = res + 1) # y values
z <- matrix(seq(min(imp), max(imp), length.out = res), nrow = 1, ncol = res) # color gradient
image(x = c(0, 1), y = height, z = z, col = palette, 
      axes = FALSE, xlab = "", ylab = "") # make color bar
ticks <- seq(min(imp), max(imp), length.out = 5) # ticks
axis(4, at = ticks, labels = round(ticks, 2), las = 1) # y axis

# title
mtext("Parameter estimates for averaged model of mammals", 
      outer = TRUE, cex = 1.1, line = 0, font = 2, family = "sans", 
      adj = 0.35)
par(mar = c(5, 4, 4, 2) + 0.1) 



















































library(MuMIn) 
library(viridis)
models <- readRDS("../results/Actinopterygii.models.rds")
models <- models[1:length(which(cumsum(models$weight) <= 0.95))] # confidence set
avg <- model.avg(models) # average

terms <- c("age.dna", "prop.dna", "age.dna:prop.dna", 
           "age.line", "prop.line", "age.line:prop.line", 
           "age.ltr", "prop.ltr", "age.ltr:prop.ltr", 
           "age.sine", "prop.sine", "age.sine:prop.sine", 
           "age.unknown", "prop.unknown", "age.unknown:prop.unknown", 
           "age.others", "prop.others", "age.others:prop.others"
)
est <- coef(avg)[match(terms, names(coef(avg)))] # parameter estimates
int <- confint(avg, full = FALSE)[terms, ] # confidence intervals
imp <- sw(models)[1:length(terms)] # importance
imp <- imp[match(terms, names(imp))] # reorder importance
x <- 1:length(terms) # x positions

# color mapping
res <- 1000 # resolution
palette <- magma(res) # palette
cols <- palette[round(((imp - min(imp)) / diff(range(imp))) * (res-1)) + 1] # colors

# x labels
labels <- c("DNA age", "DNA proportion", "DNA interaction", 
            "LINE age", "LINE proportion", "LINE interaction", 
            "LTR age", "LTR proportion", "LTR interaction", 
            "SINE age", "SINE proportion", "SINE interaction", 
            "Unidtf age", "Unidtf proportion", "Unidtf interaction", 
            "Other age", "Other proportion", "Other interaction" 
)

par(oma = c(0, 0, 3, 0))
layout(matrix(1:2, ncol = 2), widths = c(4, 1)) # make 2 plots

# main plot
par(mar = c(8, 4, 1, 0))
plot(y = est, x = x, type = "n", ylim = 1.05 * range(int), 
     ylab = "Parameter estimate", xlab = NA, axes = FALSE,
     xlim = c(0.75, length(terms) + 0.2)) # plot
abline(h = 0, lty = 2, col = "gray") # line at y = 0
segments(x, int[, 1], x, int[, 2], lwd = 2) # confidence bars
points(x, est, pch = 16, col = cols) # colored points
axis(2) # y axis
axis(1, at = x, labels = labels, las = 2) # x axis
box()

# color bar
par(mar = c(8, 1, 1, 4))
height <- seq(min(imp), max(imp), length.out = res + 1) # y values
z <- matrix(seq(min(imp), max(imp), length.out = res), nrow = 1, ncol = res) # color gradient
image(x = c(0, 1), y = height, z = z, col = palette, 
      axes = FALSE, xlab = "", ylab = "") # make color bar
ticks <- seq(min(imp), max(imp), length.out = 5) # ticks
axis(4, at = ticks, labels = round(ticks, 2), las = 1) # y axis

# title
mtext("Parameter estimates for averaged model of ray-finned fish", 
      outer = TRUE, cex = 1.1, line = 0, font = 2, family = "sans", 
      adj = 0.35)
par(mar = c(5, 4, 4, 2) + 0.1) 


















terms <- c("age.dna", "prop.dna", "age.dna:prop.dna", 
           "age.line", "prop.line", "age.line:prop.line", 
           "age.ltr", "prop.ltr", "age.ltr:prop.ltr", 
           "age.sine", "prop.sine", "age.sine:prop.sine", 
           "age.unknown", "prop.unknown", "age.unknown:prop.unknown", 
           "age.others", "prop.others", "age.others:prop.others"
)


a <- readRDS("../results/dna.line.sauria.models.rds")
b <- readRDS("../results/dna.ltr.sauria.models.rds")
c <- readRDS("../results/dna.sine.sauria.models.rds")
d <- readRDS("../results/dna.unknown.sauria.models.rds")
e <- readRDS("../results/dna.others.sauria.models.rds")
f <- readRDS("../results/line.ltr.sauria.models.rds")
g <- readRDS("../results/line.sine.sauria.models.rds")
h <- readRDS("../results/line.unknown.sauria.models.rds")
i <- readRDS("../results/line.others.sauria.models.rds")
j <- readRDS("../results/ltr.sine.sauria.models.rds")
k <- readRDS("../results/ltr.unknown.sauria.models.rds")
l <- readRDS("../results/ltr.others.sauria.models.rds")
m <- readRDS("../results/sine.unknown.sauria.models.rds")
n <- readRDS("../results/sine.others.sauria.models.rds")
o <- readRDS("../results/unknown.others.sauria.models.rds")

a <- a[1:length(which(cumsum(a$weight) <= 0.95))]
b <- b[1:length(which(cumsum(b$weight) <= 0.95))]
c <- c[1:length(which(cumsum(c$weight) <= 0.95))]
d <- d[1:length(which(cumsum(d$weight) <= 0.95))]
e <- e[1:length(which(cumsum(e$weight) <= 0.95))]
f <- f[1:length(which(cumsum(f$weight) <= 0.95))]
g <- g[1:length(which(cumsum(g$weight) <= 0.95))]
h <- h[1:length(which(cumsum(h$weight) <= 0.95))]
i <- i[1:length(which(cumsum(i$weight) <= 0.95))]
j <- j[1:length(which(cumsum(j$weight) <= 0.95))]
k <- k[1:length(which(cumsum(k$weight) <= 0.95))]
l <- l[1:length(which(cumsum(l$weight) <= 0.95))]
m <- m[1:length(which(cumsum(m$weight) <= 0.95))]
n <- n[1:length(which(cumsum(n$weight) <= 0.95))]
o <- o[1:length(which(cumsum(o$weight) <= 0.95))]



df <- data.frame(
  dna.line = sw(a)[terms],
  dna.ltr = sw(b)[terms],
  dna.sine = sw(c)[terms], 
  dna.unknown = sw(d)[terms], 
  dna.others = sw(e)[terms], 
  line.ltr = sw(f)[terms], 
  line.sine = sw(g)[terms], 
  line.unknown = sw(h)[terms], 
  line.others = sw(i)[terms], 
  ltr.sine = sw(j)[terms], 
  ltr.unknown = sw(k)[terms], 
  ltr.others = sw(l)[terms], 
  sine.unknown = sw(m)[terms], 
  sine.others = sw(n)[terms], 
  unknown.others = sw(o)[terms]
)
row.names(df) <- terms

vec <- c()
for (i in 1:nrow(df)) {
  vec <- c(vec, mean(df[i, ][!is.na(df[i, ])]))
}
names(vec) <- terms
sort(vec)













































library(MuMIn) 
library(viridis)
models <- readRDS("../results/dna.unknown.Sauria.models.rds")
models <- models[1:length(which(cumsum(models$weight) <= 0.95))] # confidence set
avg <- model.avg(models) # average

terms <- c("age.dna", "prop.dna", 
           "age.line", "prop.line", "age.line:prop.line", 
           "age.ltr", "prop.ltr", "age.ltr:prop.ltr", 
           "age.sine", "prop.sine", "age.sine:prop.sine", 
           "age.unknown", "prop.unknown", 
           "age.others", "prop.others", "age.others:prop.others"
)
est <- coef(avg)[match(terms, names(coef(avg)))] # parameter estimates
int <- confint(avg, full = FALSE)[terms, ] # confidence intervals
imp <- sw(models)[1:length(terms)] # importance
imp <- imp[match(terms, names(imp))] # reorder importance
x <- 1:length(terms) # x positions

# color mapping
res <- 1000 # resolution
palette <- magma(res) # palette
cols <- palette[round(((imp - min(imp)) / diff(range(imp))) * (res-1)) + 1] # colors

# x labels
labels <- c("DNA age", "DNA proportion", 
            "LINE age", "LINE proportion", "LINE interaction", 
            "LTR age", "LTR proportion", "LTR interaction", 
            "SINE age", "SINE proportion", "SINE interaction", 
            "Unidtf age", "Unidtf proportion", 
            "Other age", "Other proportion", "Other interaction" 
)

par(oma = c(0, 0, 3, 0))
layout(matrix(1:2, ncol = 2), widths = c(4, 1)) # make 2 plots

# main plot
par(mar = c(8, 4, 1, 0))
plot(y = est, x = x, type = "n", ylim = 1.05 * range(int), 
     ylab = "Parameter estimate", xlab = NA, axes = FALSE,
     xlim = c(0.75, length(terms) + 0.2)) # plot
abline(h = 0, lty = 2, col = "gray") # line at y = 0
segments(x, int[, 1], x, int[, 2], lwd = 2) # confidence bars
points(x, est, pch = 16, col = cols) # colored points
axis(2) # y axis
axis(1, at = x, labels = labels, las = 2) # x axis
box()

# color bar
par(mar = c(8, 1, 1, 4))
height <- seq(min(imp), max(imp), length.out = res + 1) # y values
z <- matrix(seq(min(imp), max(imp), length.out = res), nrow = 1, ncol = res) # color gradient
image(x = c(0, 1), y = height, z = z, col = palette, 
      axes = FALSE, xlab = "", ylab = "") # make color bar
ticks <- seq(min(imp), max(imp), length.out = 5) # ticks
axis(4, at = ticks, labels = round(ticks, 2), las = 1) # y axis

# title
mtext("Parameter estimates for averaged model of reptiles", 
      outer = TRUE, cex = 1.1, line = 0, font = 2, family = "sans", 
      adj = 0.4)
par(mar = c(5, 4, 4, 2) + 0.1) 
