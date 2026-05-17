library(viridis)

sub.df <- read.csv("../results/subsample.csv")

sub.df$predictor <- sapply(sub.df$predictor, function(t) {
  if (grepl(":", t)) paste(sort(strsplit(t, ":")[[1]]), collapse = ":") else t
})

# flip the signs of gini and CV so larger number -> more homogeneous
to.flip <- sub.df$response %in% c("gini", "cv")
sub.df$lower.flip <- ifelse(to.flip, -sub.df$upper, sub.df$lower)
sub.df$upper.flip <- ifelse(to.flip, -sub.df$lower, sub.df$upper)

sub.df$sig.iter <- !is.na(sub.df$importance) & sub.df$importance > 0.5 & 
                   !is.na(sub.df$lower.flip) & !is.na(sub.df$upper.flip) & 
                   sign(sub.df$lower.flip) == sign(sub.df$upper.flip)

# summarize per (clade, response, predictor) -> % significant and median importance
pct.df <- do.call(rbind, lapply(
  split(sub.df, list(sub.df$clade, sub.df$response, sub.df$predictor), drop = TRUE), 
  function(g) {
    if (nrow(g) == 0) return(NULL)
    data.frame(
      clade = g$clade[1], 
      response = g$response[1], 
      predictor = g$predictor[1], 
      pct.sig = mean(g$sig.iter, na.rm = TRUE) * 100, 
      med.imp = median(g$importance, na.rm = TRUE)
    )
  }))
rownames(pct.df) <- NULL

clades <- c("All", "Mammalia", "Actinopterygii", "Sauropsida")
response.pch <- c(rsq = 21, gini = 22, cv = 24)

pct.df$clade <- factor(pct.df$clade, levels = clades)
pct.df$response <- factor(pct.df$response, levels = c("rsq", "gini", "cv"))
pct.df$pch <- response.pch[as.character(pct.df$response)]

ncol <- 10000
pal <- viridis(ncol, begin = 0, end = 0.8, option = "A")
imp.min <- 0.5
imp.max <- 1.0
imp <- pct.df$med.imp
idx <- round(((imp - imp.min) / (imp.max - imp.min)) * (ncol - 1)) + 1
pct.df$point.col <- pal[pmax(1, pmin(ncol, idx))]

par(mfrow = c(1, 1))
for (cl in clades) {
  dat <- pct.df[pct.df$clade == cl, ]

  par(mar = c(7, 5, 3, 10) + 0.1)

  unique.preds <- unique(dat$predictor)
  responses <- levels(droplevels(dat$response))
  n.pred <- length(unique.preds)
  n.resp <- length(responses)

  group.w <- 0.7
  offsets <- setNames(seq(-group.w/2, group.w/2, length.out = n.resp), responses)

  plot(1, 1, type = "n",
       xlim = c(0.5, n.pred + 0.5), ylim = c(-5, 105), 
       xlab = NA, ylab = "% of subsamples where predictor is significant", 
       axes = FALSE, mgp = c(3.5, 0, 0))

  title(main = paste0("Subsample robustness (n=1000): ", cl), 
        line = 0.8, adj = 0, cex.main = 1.1)

  abline(h = c(0, 50, 100), lty = c(1, 3, 1), lwd = c(0.8, 1, 0.8), 
         col = c("grey80", "grey50", "grey80"))

  for (p in seq_len(n.pred - 1)) {
    abline(v = p + 0.5, lty = 3, lwd = 0.5, col = "grey85")
  }

  # labels
  xlabs <- c()
  for (i in unique.preds) {
    rep <- toupper(regmatches(i, regexpr("(?<=\\.)[a-zA-Z]+", i, perl = TRUE)))
    if (length(rep) == 0) rep <- ""
    if (rep == "OTHERS") rep <- "Others"
    if (rep == "UNKNOWN") rep <- "Unknown"
    if (rep == "TOTAL") rep <- "All repeats"
    if (grepl(":", i)) {
      xlabs <- c(xlabs, paste(rep, "age x prop."))
      next
    }
    type <- if (sub("\\..*", "", i) == "prop") "proportion" else "age"
    xlabs <- c(xlabs, paste(rep, type))
  }

  for (p.idx in seq_along(unique.preds)) {
    pred <- unique.preds[p.idx]
    for (r in responses) {
      row <- dat[dat$response == r & dat$predictor == pred, ]
      if (nrow(row) == 0) next
      x <- p.idx + offsets[r]
      segments(x, 0, x, row$pct.sig, lwd = 1, col = "grey60")
      points(x, row$pct.sig, pch = row$pch, cex = 1.8, 
             bg = row$point.col, col = "black", lwd = 0.7)
    }
  }

  axis(1, at = seq_along(unique.preds), labels = xlabs, las = 2, cex.axis = 0.88)
  axis(2, at = seq(0, 100, by = 25), mgp = c(1, 0.8, 0))
  box()

  usr <- par("usr")
  par(xpd = NA)
  lx <- usr[2] + 0.04 * diff(usr[1:2])

  legend(lx, usr[4],
         legend = c("R²", "Gini", "CV"), 
         col = "black", 
         pch = response.pch[c("rsq", "gini", "cv")], 
         pt.cex = 1.6, 
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
       "Median importance", adj = c(0.5, 0), cex = 0.75)

  par(xpd = FALSE)
}
