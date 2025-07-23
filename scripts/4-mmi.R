# Zhaobo Hu
# zhaobohu2002@gmail.com

# Finds significant models and predictors of R2 using IT-based model averaging. 
# Safe to run on low observation count. Can choose between glm() and pgls() models
# and optionally generate a csv of all significant models


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
  repeat.types <- unique(sub("^[^.]*\\.", "", variables))
  interactions <- paste0("age.", repeat.types, ":prop.", repeat.types)
  all.terms <- c(
    variables, 
    interactions
  )
  
  # pgls
  # cd <- comparative.data(tree, dat, names.col = "species", vcv = TRUE)
  
  # n <- nrow(cd$data)
  n <- nrow(dat)
  max.vars <- n-2 # leave 1 residual degree of freedom
  
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
  
  # calculate CIs
  ci.table <- data.frame(matrix(NA, nrow = length(all.terms), ncol = 3))
  colnames(ci.table) <- c("importance", "lower", "upper")
  rownames(ci.table) <- all.terms
  ci.table[names(sw(models)), 1] <- sw(models)
  ci <- confint(model.avg(models))
  ci.found <- intersect(row.names(ci), all.terms)
  ci.table[ci.found, 2:3] <- ci[ci.found, ]
  
  # subset significant models
  idx <- which(sign(ci.table[, 2]) == sign(ci.table[, 3])) # idx where 0 is not in ci
  ci.table <- ci.table[idx, ] # subset ci
  
  # append to results
  df <- data.frame(clade, 
                   rownames(ci.table), 
                   nrow(models), 
                   sapply(1:nrow(ci.table), 
                          function(x) mean(unlist(ci.table[x, 2:3]))), 
                   ci.table)
  colnames(df) <- c("clade", "model", "num.models", "estimate", "importance", "lower", "upper")
  if (nrow(combined.df) == 0) {
    combined.df <- df
  } else {
    combined.df <- rbind(combined.df, df)
  }
  
  # add lambda p values and save table
  # lambda.p <- c()
  # for (n in 1:nrow(models)) {
  #   l <- get.models(models, subset = n)[[1]]
  #   lambda.p[n] <- lambda.test(l)
  # }
  # models$lambda.p <- lambda.p
  # write.csv(data.frame(models), paste0("../results/", tolower(clade), "-models.csv"), row.names = F)
}
combined.df <- combined.df[combined.df$importance > 0.5, ]
write.csv(combined.df, "../results/model-averaging.csv", row.names = F)

