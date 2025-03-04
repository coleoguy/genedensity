


library(MuMIn)
dat <- readRDS("../results/all.8.rds")
dat <- dat[!duplicated(dat)]
assignInNamespace(
  ".modelNames",
  function(allTerms, uqTerms) {
    paste0("model", seq_along(allTerms))
  },
  ns = "MuMIn"
)
avg <- model.avg(dat, subset = delta <= 10)

int <- confint(avg, full = FALSE)
estimates <- coef(avg)
terms <- c("Intercept", "Age", "Proportion", "Interaction")
y <- seq_along(estimates)
plot(y = estimates, x = y, ylim = range(int), pch = 19, 
     ylab = "Parameter estimate", xlab = NA, axes = FALSE,
     main = "Parameter estimates for model of all species", 
     xlim = c(0.5, 4.5))
segments(y, int[, 1], y, int[, 2], lwd = 2)
axis(2)
axis(1, at = y, labels = terms, las = 2)
abline(h = 0, lty = 2, col = "gray")
box()














dat <- readRDS("../results/mammalia.8.rds")
dat <- dat[!duplicated(dat)]
assignInNamespace(
  ".modelNames",
  function(allTerms, uqTerms) {
    paste0("model", seq_along(allTerms))
  },
  ns = "MuMIn"
)
avg <- model.avg(dat, subset = delta <= 10)

int <- confint(avg, full = FALSE)
estimates <- coef(avg)
terms <- c("Intercept", "Age", "Proportion", "Interaction")
y <- seq_along(estimates)
plot(y = estimates, x = y, ylim = range(int), pch = 19, 
     ylab = "Parameter estimate", xlab = NA, axes = FALSE,
     main = "Parameter estimates for model of mammals", 
     xlim = c(0.5, 4.5))
segments(y, int[, 1], y, int[, 2], lwd = 2)
axis(2)
axis(1, at = y, labels = terms, las = 2)
abline(h = 0, lty = 2, col = "gray")
box()












dat <- readRDS("../results/actinopterygii.8.rds")
dat <- dat[!duplicated(dat)]
assignInNamespace(
  ".modelNames",
  function(allTerms, uqTerms) {
    paste0("model", seq_along(allTerms))
  },
  ns = "MuMIn"
)
avg <- model.avg(dat, subset = delta <= 10)

int <- confint(avg, full = FALSE)
estimates <- coef(avg)
terms <- c("Intercept", "Age", "Proportion", "Interaction")
y <- seq_along(estimates)
plot(y = estimates, x = y, ylim = range(int), pch = 19, 
     ylab = "Parameter estimate", xlab = NA, axes = FALSE,
     main = "Parameter estimates for model of fish", 
     xlim = c(0.5, 4.5))
segments(y, int[, 1], y, int[, 2], lwd = 2)
axis(2)
axis(1, at = y, labels = terms, las = 2)
abline(h = 0, lty = 2, col = "gray")
box()




















dat <- readRDS("../results/sauria.8.rds")
dat <- dat[!duplicated(dat)]
assignInNamespace(
  ".modelNames",
  function(allTerms, uqTerms) {
    paste0("model", seq_along(allTerms))
  },
  ns = "MuMIn"
)
avg <- model.avg(dat, subset = delta <= 10)

int <- confint(avg, full = FALSE)
estimates <- coef(avg)
terms <- c("Intercept", "Age", "Proportion", "Interaction")
y <- seq_along(estimates)
plot(y = estimates, x = y, ylim = range(int), pch = 19, 
     ylab = "Parameter estimate", xlab = NA, axes = FALSE,
     main = "Parameter estimates for model of reptiles", 
     xlim = c(0.5, 4.5))
segments(y, int[, 1], y, int[, 2], lwd = 2)
axis(2)
axis(1, at = y, labels = terms, las = 2)
abline(h = 0, lty = 2, col = "gray")
box()

