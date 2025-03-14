



library(MuMIn) 
library(viridis)
models <- readRDS("../results/Actinopterygii.rds")
models <- models[1:length(which(cumsum(models$weight) <= 0.95))] # confidence set
avg <- model.avg(models) # average

imp <- sort(sw(models), decreasing = T) # importance
terms <- names(imp) # terms
est <- coef(avg)[match(terms, names(coef(avg)))] # parameter estimates
int <- as.data.frame(confint(avg, full = FALSE)) [terms, ] # confidence intervals
x <- 1:length(terms) # x positions

# color mapping
res <- 1000 # resolution
palette <- viridis(res, begin = 0, end = 0.8, option = "A") # palette
cols <- palette[round(((imp - min(imp)) / diff(range(imp))) * (res-1)) + 1] # colors

# x labels
labels <- c()
for (i in terms) {
  rep <- toupper(regmatches(i, regexpr("(?<=\\.)[a-zA-Z]+", i, perl = TRUE)))
  if (rep == "OTHERS") {
    rep <- "Others"
  } else if (rep == "UNKNOWN") {
    rep <- "Unidtf"
  }
  
  if (grepl(":", i)) {
    type <- "int."
  } else if (sub("\\..*", "", i) == "prop") {
    type <- "prop."
  } else {
    type <- "age"
  }
  labels <- c(labels, paste(rep, type))
}

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


















library(MuMIn) 
library(viridis)
models <- readRDS("../results/Mammalia.rds")
models <- models[1:length(which(cumsum(models$weight) <= 0.95))] # confidence set
avg <- model.avg(models) # average

imp <- sort(sw(models), decreasing = T) # importance
terms <- names(imp) # terms
est <- coef(avg)[match(terms, names(coef(avg)))] # parameter estimates
int <- as.data.frame(confint(avg, full = FALSE)) [terms, ] # confidence intervals
x <- 1:length(terms) # x positions

# color mapping
res <- 1000 # resolution
palette <- viridis(res, begin = 0, end = 0.8, option = "A") # palette
cols <- palette[round(((imp - min(imp)) / diff(range(imp))) * (res-1)) + 1] # colors

# x labels
labels <- c()
for (i in terms) {
  rep <- toupper(regmatches(i, regexpr("(?<=\\.)[a-zA-Z]+", i, perl = TRUE)))
  if (rep == "OTHERS") {
    rep <- "Others"
  } else if (rep == "UNKNOWN") {
    rep <- "Unidtf"
  }
  
  if (grepl(":", i)) {
    type <- "int."
  } else if (sub("\\..*", "", i) == "prop") {
    type <- "prop."
  } else {
    type <- "age"
  }
  labels <- c(labels, paste(rep, type))
}

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
models <- readRDS("../results/Mammalia.models.rds")
models <- models[1:length(which(cumsum(models$weight) <= 0.95))] # confidence set
avg <- model.avg(models) # average

imp <- sort(sw(models), decreasing = T) # importance
terms <- names(imp) # terms
est <- coef(avg)[match(terms, names(coef(avg)))] # parameter estimates
int <- as.data.frame(confint(avg, full = FALSE)) [terms, ] # confidence intervals
x <- 1:length(terms) # x positions

# color mapping
res <- 1000 # resolution
palette <- viridis(res, begin = 0, end = 0.8, option = "A") # palette

cols <- palette[round(((imp - min(imp)) / diff(range(imp))) * (res-1)) + 1] # colors

# x labels
labels <- terms

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

imp <- sort(sw(models), decreasing = T) # importance
terms <- names(imp) # terms
est <- coef(avg)[match(terms, names(coef(avg)))] # parameter estimates
int <- as.data.frame(confint(avg, full = FALSE)) [terms, ] # confidence intervals
x <- 1:length(terms) # x positions

# color mapping
res <- 1000 # resolution
palette <- viridis(res, begin = 0, end = 0.8, option = "A") # palette

cols <- palette[round(((imp - min(imp)) / diff(range(imp))) * (res-1)) + 1] # colors

# x labels
labels <- terms

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
















































