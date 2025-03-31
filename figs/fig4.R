


library(MuMIn) 
library(viridis)
all <- readRDS("../results/all.rds")
mam <- readRDS("../results/mammalia.rds")
fish <- readRDS("../results/actinopterygii.rds")
rep <- readRDS("../results/sauria.rds")

all <- all[order(all$AICc), ]
mam <- mam[order(mam$AICc), ]
fish <- fish[order(fish$AICc), ]
rep <- rep[order(rep$AICc), ]

# remove seemingly overparametrized model and recalculate AICc
rep <- rep[-1, ]
rep <- rbind(rep, rep[1,])
rep <- rep[-1, ]

all <- all[cumsum(all$weight) <= 0.95, ]
mam <- mam[cumsum(mam$weight) <= 0.95, ]
fish <- fish[cumsum(fish$weight) <= 0.95, ]
rep <- rep[cumsum(rep$weight) <= 0.95, ]

combined.df <- data.frame()
for (i in c("all", "mam", "fish", "rep")) {
  assign(i, model.avg(get(i))) # average
  imp <- sort(sw(get(i)), decreasing = TRUE) # get importance
  ci <- confint(get(i)) # get ci
  ci <- ci[match(names(imp), row.names(ci)), ] #match ci
  ci <- as.data.frame(ci)
  
  idx <- which(sign(ci[, 1]) == sign(ci[, 2])) # idx where 0 is not in ci
  
  ci <- ci[idx, ]# subset ci
  imp <- imp[idx]# subset importance
  
  df <- data.frame(i, 
                   names(imp), 
                   sapply(1:nrow(ci), function(x) mean(unlist(ci[x, ]))), 
                   imp, 
                   ci[, 1], 
                   ci[, 2])
  colnames(df) <- c("clade", "model", "estimate", "importance", "lower", "upper")
  
  if (is.null(combined.df)) {
    df <- combined.df
  } else {
    combined.df <- rbind(combined.df, df)
  }
}

combined.df <- combined.df[combined.df$importance >= 0.5, ]

# x positions
x <- c()
for (i in 1:length(combined.df$clade)) {
  if (is.null(x)) {
    x <- c(1)
  } else {
    x <- c(x, ifelse(combined.df$clade[i] == prev, tail(x, 1) + 0.55, tail(x, 1) + 1))
  }
  prev <- combined.df$clade[i]
}
imp <- combined.df$importance

# color mapping
res <- 10000 # resolution
palette <- viridis(res, begin = 0, end = 0.8, option = "A") # palette
cols <- palette[round(((imp - min(imp)) / diff(range(imp))) * (res-1)) + 1] # colors

# x labels
labels <- c()
for (i in combined.df$model) {
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
int.range <- range(as.matrix(combined.df[c("lower", "upper")]))
plot(y = combined.df$estimate, x = x, type = "n", ylim = 1.05 * int.range, 
     ylab = "Parameter estimate", xlab = NA, axes = FALSE,
     xlim = c(min(x)-0.25, max(x)+0.25)) # plot
abline(h = 0, lty = 1, col = "black") # line at y = 0
for (l in -100:100) {
  abline(h = l, lty = 2, col = "grey") # line at y = 0
}
segments(x, combined.df$lower, x, combined.df$upper, lwd = 2) # confidence bars
segments(x-0.1, combined.df$upper, x+0.1, combined.df$upper, lwd = 2)
segments(x-0.1, combined.df$lower, x+0.1, combined.df$lower, lwd = 2)
points(x, combined.df$estimate, pch = 16, cex = 1.5, col = cols) # colored points
axis(2) # y axis
axis(2, at = seq(-10, 10, by = 0.5), labels = FALSE, tcl = -0.2)
axis(2, at = seq(-10, 10, by = 1), labels = FALSE, tcl = -0.5)
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
mtext("Parameter estimates for averaged models", 
      outer = TRUE, cex = 1.1, line = 0, font = 2, family = "sans", 
      adj = 0.35)
par(mar = c(5, 4, 4, 2) + 0.1) 

















































