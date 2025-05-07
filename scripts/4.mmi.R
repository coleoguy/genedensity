
library(MuMIn) # multimodel inference
library(phytools) # load and prune tree
library(caper) # PGLS
options(na.action = "na.fail")  

# functions
source("functions.R")

rsq <- read.csv("../results/rsq.csv")
repeats <- read.csv("../results/repeat.results.csv")

combined.df <- data.frame()
# loop for each clade
for (i in c("All", "Mammalia", "Actinopterygii", "Sauropsida")) {
  
  # subset results
  dat <- merge(rsq, repeats, by.x = "species", by.y = "species", all.x = T, all.y = T)
  if (i %in% c("Mammalia", "Actinopterygii", "Sauropsida")) {
    dat <- dat[dat$clade %in% i, ]
  }
  variables <- colnames(dat)[grep("^(prop|age)\\.", colnames(dat))]
  dat <- na.omit(dat[, c("species", "clade", "rsq", variables)])
  
  # prune tree
  # tree <- read.tree("../data/formatted.tree.nwk")
  # int <- intersect(dat$species, tree$tip.label)
  
  # rescale
  for (j in variables) {
    dat[[j]] <- (dat[[j]]-min(dat[[j]])) / diff(range(dat[[j]]))
  }
  # cd <- comparative.data(pruned.tree, dat, names.col = "species", vcv = TRUE)
  
  # set up interactions
  rep <- unique(sub("^[^.]*\\.", "", variables))
  interactions <- paste0("age.", rep, ":prop.", rep)
  all.terms <- c(
    variables, 
    interactions
  )
  
  # global model
  global.model <- glm(reformulate(all.terms, response = "rsq"), data = dat)
  
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
  models <- dredge(global.model, 
                   subset = subset.expr
                   # extra = list(shapirowilk.p = sw.test, lambda.p = lambda.test)
  )
  models <- models[order(models$AICc), ]
  models <- models[cumsum(models$weight) <= 0.95, ]
  imp <- sort(sw(models), decreasing = TRUE)
  avg <- model.avg(models)
  ci <- confint(avg)
  ci <- ci[match(names(imp), row.names(ci)), ] #match ci
  ci <- as.data.frame(ci)
  
  idx <- which(sign(ci[, 1]) == sign(ci[, 2])) # idx where 0 is not in ci
  ci <- ci[idx, ]# subset ci
  imp <- imp[idx]# subset importance
  if (length(imp) == 0) {
    next
  }
  df <- data.frame(i, 
                   names(imp), 
                   sapply(1:nrow(ci), function(x) mean(unlist(ci[x, ]))), 
                   imp, 
                   ci[, 1], 
                   ci[, 2])
  colnames(df) <- c("clade", "model", "estimate", "importance", "lower", "upper")
  
  if (nrow(combined.df) == 0) {
    combined.df <- df
  } else {
    combined.df <- rbind(combined.df, df)
  }
}
combined.df <- combined.df[combined.df$importance > 0.5, ]
write.csv(combined.df, "../results/model.averaging.csv", row.names = F)


