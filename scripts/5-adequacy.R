# Zhaobo Hu
# zhaobohu2002@gmail.com

# A permutation-based approach to assess model adequacy


library(MuMIn) # multimodel inference
library(phytools) # load and prune tree
library(caper) # PGLS
options(na.action = "na.fail")  

# functions
source("functions.R")

rsq <- read.csv("../results/rsq.csv")
repeats <- read.csv("../results/repeat-results.csv")
tree <- read.tree("../data/formatted-tree.nwk")

main <- merge(rsq, repeats, by.x = "species", by.y = "species", all.x = T, all.y = T)
variables <- colnames(main)[grep("^(prop|age)\\.", colnames(main))]
main <- na.omit(main[, c("species", "clade", "rsq", variables)])
clades <- c("All", "Mammalia", "Actinopterygii", "Sauropsida")
constant.cols <- c("species", "clade", "rsq")

for (h in 1:3) { # for each run
  
  print(h)
  run.results <- c()
  
  # permute
  perm.colname <- setdiff(names(main), constant.cols)
  block <- main[, perm.colname]
  
  # block
  block <- block[sample(nrow(block)), ]
  main[, perm.colname] <- block
  
  for (i in 1:4) { # for each clade
    dat <- main
    clade <- clades[i]
    
    # subset results
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
    for (k in 1:length(all.terms)) {
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
    ci <- ci[match(all.terms, row.names(ci)), ] # match ci
    lower <- ci[, 1]
    upper <- ci[, 2]
    run.results <- c(run.results, imp, upper, lower)
    
    # make initial df 
    if (h == 1 && i == 1) {
      df <- data.frame(
        dataset = rep(clades, each = length(imp) * 3), 
        stat = rep(rep(c("importance", "upper", "lower"), each = length(imp)), 4), 
        variable = rep(rep(names(imp), 3), 4)
      )
    }
  }
  df[[paste0("run", h)]] <- run.results
  # write.csv(df, "../results/permute.csv", row.names = F)
  gc()
}
write.csv(df, "../results/permute.csv", row.names = F)

