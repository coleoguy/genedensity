

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
           "UnID age", "UnID proportion", "UnID interaction", 
           "Other age", "Other proportion", "Other interaction" 
)

# plot
par(mar = c(8, 4, 4, 2))
plot(y = est, x = x, type = "n", ylim = 1.05 * range(int), 
     ylab = "Parameter estimate", xlab = NA, axes = FALSE,
     main = "Parameter estimates for averaged model of all species", 
     xlim = c(0.75, length(terms)+0.2))
abline(h = 0, lty = 2, col = "gray") # dashed line at y = 0
segments(x, int[, 1], x, int[, 2], lwd = 2) # confidence bars
points(x, est, pch = 16, col = cols) # colored points
axis(2) # y axis
axis(1, at = x, labels = labels, las = 2) # x axis
box() # box
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
            "UnID age", "UnID proportion", "UnID interaction", 
            "Other age", "Other proportion", "Other interaction" 
)

# plot
par(mar = c(8, 4, 4, 2))
plot(y = est, x = x, type = "n", ylim = 1.05 * range(int), 
     ylab = "Parameter estimate", xlab = NA, axes = FALSE,
     main = "Parameter estimates for averaged model of mammals", 
     xlim = c(0.75, length(terms)+0.2))
abline(h = 0, lty = 2, col = "gray") # dashed line at y = 0
segments(x, int[, 1], x, int[, 2], lwd = 2) # confidence bars
points(x, est, pch = 16, col = cols) # colored points
axis(2) # y axis
axis(1, at = x, labels = labels, las = 2) # x axis
box() # box
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
            "UnID age", "UnID proportion", "UnID interaction", 
            "Other age", "Other proportion", "Other interaction" 
)

# plot
par(mar = c(8, 4, 4, 2))
plot(y = est, x = x, type = "n", ylim = 1.05 * range(int), 
     ylab = "Parameter estimate", xlab = NA, axes = FALSE,
     main = "Parameter estimates for averaged model of ray-finned fish", 
     xlim = c(0.75, length(terms)+0.2))
abline(h = 0, lty = 2, col = "gray") # dashed line at y = 0
segments(x, int[, 1], x, int[, 2], lwd = 2) # confidence bars
points(x, est, pch = 16, col = cols) # colored points
axis(2) # y axis
axis(1, at = x, labels = labels, las = 2) # x axis
box() # box
par(mar = c(5, 4, 4, 2) + 0.1) 
