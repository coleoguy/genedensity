
library(MuMIn) # multimodel inference
library(phytools) # load and prune tree
library(caper) # PGLS
options(na.action = "na.fail")  

# functions
source("functions.R")

rsq <- read.csv("../results/rsq.csv")
repeats <- read.csv("../results/repeat-results.csv")
tree <- read.tree("../data/formatted-tree.nwk")

dat <- merge(rsq, repeats, by.x = "species", by.y = "species", all.x = T, all.y = T)
variables <- colnames(dat)[grep("^(prop|age)\\.", colnames(dat))]
dat <- na.omit(dat[, c("species", "clade", "rsq", variables)])
clades <- c("All", "Mammalia", "Actinopterygii", "Sauropsida")

for (h in 1:10) { # for each run
  
  print(h)
  run.results <- c()
  
  # permute
  constant.cols <- c("species", "clade", "rsq")
  perm.colname <- setdiff(names(dat), constant.cols)
  block <- dat[, perm.colname]
  for (l in 1:length(block)) {
    block[, l] <- sample(block[, l])
  }
  dat[, perm.colname] <- block
  
  for (i in 1:4) { # for each clade
    clade <- clades[i]
    sub <- dat
    
    # subset results
    if (clade %in% c("Mammalia", "Actinopterygii", "Sauropsida")) {
      sub <- sub[sub$clade %in% clade, ]
    }
    
    # rescale
    for (j in variables) {
      sub[[j]] <- (sub[[j]]-min(sub[[j]])) / diff(range(sub[[j]]))
    }
    
    # set up interactions
    rep <- unique(sub("^[^.]*\\.", "", variables))
    interactions <- paste0("age.", rep, ":prop.", rep)
    all.terms <- c(
      variables, 
      interactions
    )
    
    # gls
    global.model <- glm(reformulate(all.terms, response = "rsq"), data = sub)
    
    # pgls
    # cd <- comparative.data(tree, sub, names.col = "species", vcv = TRUE)
    # global.model <- pgls(reformulate(all.terms, response = "rsq"), data = cd)
    
    # set constraints
    model.terms <- unlist(strsplit(as.character(global.model$formula)[3], " \\+ "))
    model.interactions <- grep(":", model.terms, value = TRUE)
    constraints <- character(length(model.interactions))
    for (k in seq_along(model.interactions)) {
      parts <- strsplit(model.interactions[k], ":")[[1]]
      constraints[k] <- sprintf("((!`%s`) | (%s & %s))", model.interactions[k], parts[1], parts[2])
    }
    subset.expr <- parse(text = paste(constraints, collapse = " & "))[[1]]
    
    # dredge
    num <- ifelse(nrow(sub) > 21, 16, 15) # get rid of error in reptiles
    models <- dredge(global.model, 
                     subset = subset.expr, 
                     m.lim = c(0, num)
                     # extra = list(shapirowilk.p = sw.test, lambda.p = lambda.test)
    )
    models <- models[order(models$AICc), ]
    models <- models[cumsum(models$weight) <= 0.95, ]
    
    # get CIs and int in a vector
    imp <- sw(models)
    imp <- imp[match(all.terms, names(imp))]
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
  write.csv(df, "../results/permute.csv", row.names = F)
  gc()
}

