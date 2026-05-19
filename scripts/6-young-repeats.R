# Zhaobo Hu
# zhaobohu2002@gmail.com

# For each species, compute the fraction of young repeats (Kimura < cutoff)
# across a sweep of cutoffs. For each (cutoff, response, group) combination,
# fit a null model and a slope model and record their AICc comparison plus
# the slope coefficient and its CI.

exclude.classes <- c("Simple_repeat", "Satellite", "Low_complexity",
                     "ARTEFACT", "rRNA", "scRNA", "snRNA", "tRNA")

# AICc helper
aicc <- function(fit) {
  n <- length(resid(fit))
  k <- length(coef(fit)) + 1
  ll <- logLik(fit)[1]
  -2 * ll + 2 * k + (2 * k * (k + 1)) / (n - k - 1)
}
rsq.df <- read.csv("../results/rsq.csv")
rsq.df[[1]] <- tolower(gsub(" ", "_", rsq.df[[1]]))
species.col <- colnames(rsq.df)[1]

files <- list.files("../results/divsums", pattern = "\\.divsum$", 
                    full.names = T)

clades <- c("Mammalia", "Actinopterygii", "Sauropsida")
responses <- c("rsq", "gini.flip", "cv.flip")
resp.names <- c("rsq", "gini", "cv")
cutoffs<- seq(2, 20, by = 2)


sens.df <- data.frame()
# loop
for (i in cutoffs) {
  
  v <- c()
  for (f in files) {
    lines <- readLines(f)
    hist.start <- grep("^Coverage for each repeat class and divergence", lines)
    if (length(hist.start) == 0) return(NA)
    header.line <- lines[hist.start + 1]
    data.lines <- lines[(hist.start + 2):length(lines)]
    data.lines <- data.lines[nzchar(trimws(data.lines))]
    cols <- strsplit(trimws(header.line), "\\s+")[[1]]
    mat <- do.call(rbind, lapply(data.lines, function(l) {
      as.numeric(strsplit(trimws(l), "\\s+")[[1]])
    }))
    colnames(mat) <- cols
    div.bin <- mat[, "Div"]
    repeat.cols <- setdiff(colnames(mat), "Div")
    if (length(exclude.classes) > 0) {
      cls.prefix  <- sub("/.*$", "", repeat.cols)
      repeat.cols <- repeat.cols[!cls.prefix %in% exclude.classes]
    }
    if (length(repeat.cols) == 0) return(NA)
    repeat.bp <- rowSums(mat[, repeat.cols, drop = F], na.rm = T)
    total <- sum(repeat.bp)
    if (total == 0) return(NA)
    v <- c(v, sum(repeat.bp[div.bin < i]) / total)
  }

  yf <- data.frame(
    species = tolower(sub("\\.divsum$", "", basename(files))), 
    young.frac = as.numeric(v), 
    stringsAsFactors = F
  )
  yf <- yf[!is.na(yf$young.frac), ]
  yf$species <- tolower(gsub(" ", "_", yf$species))

  d <- merge(rsq.df, yf, by.x = species.col, by.y = "species")
  d$gini.flip <- -d$gini
  d$cv.flip <- -d$cv

  for (ri in seq_along(responses)) {
    resp <- responses[ri]
    for (grp in c("All", clades)) {
      sub <- if (grp == "All") d else d[d$clade == grp, ]
      if (nrow(sub) < 5) next
      y <- sub[[resp]]
      x <- sub$young.frac
      keep <- complete.cases(x, y)
      y <- y[keep]; x <- x[keep]
      # null vs slope
      m0 <- lm(y ~ 1)
      m1 <- lm(y ~ x)
      daic <- aicc(m0) - aicc(m1)
      w1   <- exp(0.5 * daic) / (1 + exp(0.5 * daic))
      beta <- coef(m1)[2]
      # 95% CI on the slope coefficient
      ci <- confint(m1)
      lo <- ci[2, 1]
      hi <- ci[2, 2]
      ci.signif <- !is.na(lo) & !is.na(hi) & sign(lo) == sign(hi)
      sens.df <- rbind(sens.df, data.frame(
        cutoff = i, 
        response = resp.names[ri], 
        group = grp, 
        n = length(y), 
        beta = beta, 
        beta.lower = lo, 
        beta.upper = hi, 
        delta.aicc = daic, 
        w1 = w1, 
        ci.signif = ci.signif, 
        supported = (w1 > 0.5) & ci.signif
      ))
    }
  }
}

write.csv(sens.df, "../results/young-repeats-sensitivity.csv", 
          row.names = FALSE)
