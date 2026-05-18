library(viridis)
library(devEMF) # svg doesn't work in powerpoint

combined.df <- read.csv("../../results/model-averaging.csv")

# keep rows where importance > 50 and CI doesn't include 0
combined.df <- combined.df[!is.na(combined.df$importance) & 
                            combined.df$importance > 0.5 & 
                            !is.na(combined.df$lower) & 
                            !is.na(combined.df$upper) & 
                            sign(combined.df$lower) == sign(combined.df$upper), ]

# flip the signs of gini and CV so larger number -> more homogeneous
to.flip <- combined.df$response %in% c("gini", "cv")
combined.df$estimate.flip <- ifelse(to.flip, -combined.df$estimate, combined.df$estimate)
combined.df$lower.flip <- ifelse(to.flip, -combined.df$upper, combined.df$lower)
combined.df$upper.flip <- ifelse(to.flip, -combined.df$lower, combined.df$upper)


clades <- c("All", "Mammalia", "Actinopterygii", "Sauropsida")
combined.df$clade <- factor(combined.df$clade, levels = clades)
combined.df$response <- factor(combined.df$response, levels = c("rsq", "gini", "cv"))

# response colors
response.colors <- c(rsq = "#1f78b4", gini = "#e31a1c", cv = "#33a02c")
response.offset <- c(rsq = 0.25, gini = 0, cv = -0.25)


# response shapes (circle, square, triangle — matches CI bar colors)
response.pch <- c(rsq = 21, gini = 22, cv = 24)
combined.df$pch <- response.pch[as.character(combined.df$response)]

# colors for importance — fixed scale 0.5 to 1 for consistency across figures
ncol <- 10000
pal <- viridis(ncol, begin = 0, end = 0.8, option = "A")
imp.min <- 0.5
imp.max <- 1.0
imp <- combined.df$importance
idx <- round(((imp - imp.min) / (imp.max - imp.min)) * (ncol - 1)) + 1
combined.df$point.col <- pal[pmax(1, pmin(ncol, idx))]

# per-clade plot
par(mfrow = c(1, 1))
for (cl in clades) {
  emf(paste0(cl, ".emf"), width = 7, height = 4)
  dat <- combined.df[combined.df$clade == cl, ]
  
  par(mar = c(5, 8, 3, 10) + 0.1)
  
  dat <- dat[order(dat$model, dat$response), ]
  unique.models <- unique(dat$model)
  y.base <- setNames(rev(seq_len(length(unique.models))), unique.models)
  dat$y.base <- y.base[as.character(dat$model)]
  dat$y <- dat$y.base + response.offset[as.character(dat$response)]
  
  ylabs <- c()
  for (i in unique.models) {
    rep <- toupper(regmatches(i, regexpr("(?<=\\.)[a-zA-Z]+", i, perl = TRUE)))
    if (length(rep) == 0) rep <- ""
    if (rep == "OTHERS") rep <- "Others"
    if (rep == "UNKNOWN") rep <- "Unknown"
    if (rep == "TOTAL") rep <- "All repeats"
    if (grepl(":", i)) {
      ylabs <- c(ylabs, paste(rep, "age x prop."))
      next
    }
    type <- if (sub("\\..*", "", i) == "prop") "proportion" else "age"
    ylabs <- c(ylabs, paste(rep, type))
  }
  
  n.models <- dat$num.models[!is.na(dat$num.models)][1]
  
  x.range <- range(c(dat$lower.flip, dat$upper.flip), na.rm = TRUE) + c(-0.5, 0.5)
  y.range <- range(dat$y) + c(-0.7, 0.7)
  
  plot(dat$estimate.flip, dat$y, type = "n", 
       xlim = x.range, ylim = y.range, 
       xlab = "Sign-aligned β", 
       ylab = NA, axes = FALSE, mgp = c(2.5, 0, 0))
  
  title.str <- if (!is.na(n.models)) paste0(cl, "  (", n.models, " models)") else cl
  title(main = title.str, line = 0.8, adj = 0, cex.main = 1.1)
  
  abline(v = 0, lty = 5, lwd = 1, col = "black")
  
  # subtle horizontal separators between predictors
  for (p in seq_len(length(unique.models) - 1)) {
    abline(h = p + 0.5, lty = 3, lwd = 0.5, col = "grey85")
  }
  
  # CI bars
  segments(dat$lower.flip, dat$y, dat$upper.flip, 
           dat$y, lwd = 1.4, col = "black")
  segments(dat$lower.flip, dat$y - 0.1, dat$lower.flip, 
           dat$y + 0.1, lwd = 1.4, col = "black")
  segments(dat$upper.flip, dat$y - 0.1, dat$upper.flip, 
           dat$y + 0.1, lwd = 1.4, col = "black")
  
  points(dat$estimate.flip, dat$y, pch = dat$pch, cex = 2.2, 
         bg = dat$point.col, col = "black", lwd = 0.8)
  
  axis(1, at = pretty(x.range), mgp = c(1, 0.8, 0))
  axis(1, at = seq(-10, 10, by = 0.5), labels = FALSE, tcl = -0.2)
  axis(1, at = seq(-10, 10, by = 1), labels = FALSE, tcl = -0.5)
  axis(2, at = y.base, labels = ylabs, las = 2, cex.axis = 0.9)
  box()
  
  usr <- par("usr")
  par(xpd = NA)
  lx <- usr[2] + 0.04 * diff(usr[1:2])
  
  legend(lx, usr[4],
         legend = c("R²", "Gini", "CV"), 
         col = "black", 
         pch = response.pch[c("rsq", "gini", "cv")], 
         lwd = 1.4, pt.cex = 1.6, 
         title = "Response metric", title.adj = 0, 
         bty = "n", cex = 0.82, xjust = 0, yjust = 1)
  
  # color bar
  margin.r.user <- diff(usr[1:2]) / par("pin")[1] * par("mai")[4]
  bx1 <- usr[2] + 0.20 * margin.r.user
  bx2 <- usr[2] + 0.92 * margin.r.user
  by1 <- usr[3] + 0.04 * diff(usr[3:4])
  by2 <- usr[3] + 0.08 * diff(usr[3:4])
  bar_img <- array(t(col2rgb(pal) / 255), c(1, ncol, 3))
  rasterImage(bar_img, bx1, by1, bx2, by2)
  rect(bx1, by1, bx2, by2, border = "black", lwd = 0.8)
  tick_at <- c(bx1, mean(c(bx1, bx2)), bx2)
  tick_lab <- c("0.5", "0.75", "1.0")
  axis(1, at = tick_at, labels = tick_lab, 
       pos = by1 - 0.015 * diff(usr[3:4]), 
       tck = -0.012, cex.axis = 0.65, mgp = c(1, -0.01, 1))
  text(mean(c(bx1, bx2)), by2 + 0.03 * diff(usr[3:4]), 
       "Variable importance", adj = c(0.5, 0), cex = 0.75)
  
  par(xpd = FALSE)
  dev.off()
}
