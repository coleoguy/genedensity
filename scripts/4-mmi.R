library(MuMIn)
library(phytools)
library(caper)
source("functions.R")

rsq <- read.csv("../results/rsq.csv")
repeats <- read.csv("../results/repeat-results.csv")
tree <- read.tree("../data/formatted-tree.nwk")

combined.df <- data.frame()
clades <- c("All", "Mammalia", "Actinopterygii", "Sauropsida")
# loop for each clade
for (i in 1:4) {
  # subset results
  clade <- clades[i]
  dat <- merge(rsq, repeats, by.x = "species", by.y = "species", all.x = T, all.y = T)
  if (clade %in% c("Mammalia", "Actinopterygii", "Sauropsida")) {
    dat <- dat[dat$clade %in% clade, ]
  }
  variables <- colnames(dat)[grep("^(prop|age)\\.", colnames(dat))]
  dat <- na.omit(dat[, c("species", "clade", "rsq", variables)])
  
  # rescale
  for (j in variables) {
    dat[[j]] <- (dat[[j]]-min(dat[[j]])) / diff(range(dat[[j]]))
  }
  
  # identify all predictors
  repeats <- unique(sub("^[^.]*\\.", "", variables))
  interactions <- paste0("age.", repeats, ":prop.", repeats)
  all.terms <- c(
    variables, 
    interactions
  )
  
  # pgls
  # cd <- comparative.data(tree, dat, names.col = "species", vcv = TRUE)
  
  # n <- nrow(cd$data)
  n <- nrow(dat)
  max.vars <- n-2 # leave 2 degrees of freedom
  
  # fit models
  model.list <- list()
  model.idx <- 1
  for (k in 0:length(all.terms)) {
    if (k == 0) {
      fml <- as.formula("rsq ~ 1")
      # fit <- pgls(fml, data = cd)
      fit <- glm(fml, data = dat)
      model.name <- paste0("M", model.idx)
      model.list[[model.name]] <- fit
      model.idx <- model.idx + 1
      next
    }
    combos <- combn(all.terms, k, simplify = FALSE)
    for (m in combos) {
      fml <- as.formula(paste("rsq ~", paste(m, collapse = " + ")))
      predictor.num <- length(attr(terms(fml), "term.labels"))
      if (predictor.num >= max.vars) next
      if (!model.marginal(m)) next
      # fit <- pgls(fml, data = cd)
      fit <- glm(fml, data = dat)
      model.name <- paste0("M", model.idx)
      model.list[[model.name]] <- fit
      model.idx <- model.idx + 1
    }
  }
  
  # create model selection table
  models <- model.sel(model.list)
  models <- models[order(models$AICc), ]
  models <- models[cumsum(models$weight) <= 0.95, ]
  rm(model.list)
  gc()
  
  # average and calculate CIs
  num <- nrow(models)
  imp <- sort(sw(models), decreasing = TRUE)
  avg <- model.avg(models)
  ci <- confint(avg)
  ci <- ci[match(names(imp), row.names(ci)), ] # match ci
  ci <- as.data.frame(ci)
  
  # subset significant models
  idx <- which(sign(ci[, 1]) == sign(ci[, 2])) # idx where 0 is not in ci
  ci <- ci[idx, ] # subset ci
  imp <- imp[idx] # subset importance
  if (length(imp) == 0) {
    next
  }
  
  # append to results
  df <- data.frame(clade, 
                   names(imp), 
                   num, 
                   sapply(1:nrow(ci), function(x) mean(unlist(ci[x, ]))), 
                   imp, 
                   ci[, 1], 
                   ci[, 2])
  colnames(df) <- c("clade", "model", "num.models", "estimate", "importance", "lower", "upper")
  if (nrow(combined.df) == 0) {
    combined.df <- df
  } else {
    combined.df <- rbind(combined.df, df)
  }
  
  # add lambda p values and save table
  new <- c()
  for (n in 1:nrow(models)) {
    l <- get.models(models, subset = n)[[1]]
    new[n] <- lambda.test(l)
  }
  models$lambda.p <- new
  write.csv(models, paste0("../results/", tolower(clade), "-models.csv"), row.names = F)
}
combined.df <- combined.df[combined.df$importance > 0.5, ]
write.csv(combined.df, "../results/model-averaging.csv", row.names = F)

