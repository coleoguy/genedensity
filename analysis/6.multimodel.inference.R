

library(MuMIn) # multimodel inference
library(phytools) # load and prune tree
library(caper) # PGLS
library(parallel)
options(na.action = "na.fail")  

# functions
sw.test <- function(model) {
  res <- residuals(model)
  sw.p <- shapiro.test(res)$p.value
  return(sw.p)
}
lambda.test <- function(model) {
  res <- setNames(residuals(model), dat$species)
  lambda.p <- phylosig(tree, res, method = "lambda", test = TRUE, niter = 100)$P
  return(lambda.p)
}

# loop for each clade
for (i in c("All", "Mammalia", "Actinopterygii", "Sauria")) {
  
  # subset results
  tree <- read.tree("../data/formatted.tree.nwk")
  tree$tip.label <- gsub("_", " ", tree$tip.label)
  dat <- read.csv("../results/parsed.csv")
  dat <- dat[!is.na(dat$chromnum.1n) & !duplicated(dat$species), ]
  if (i %in% c("Mammalia", "Actinopterygii", "Sauria")) {
    dat <- dat[dat$clade %in% i, ]
  }
  int <- intersect(dat$species, tree$tip.label)
  variables <- colnames(dat)[grep("^(prop|age)\\.", colnames(dat))]
  dat <- na.omit(dat[, c("species", 
                         "clade", 
                         "rsq", 
                         variables
  )])
  
  # normalize
  for (j in variables) {
    dat[[j]] <- (max(dat[[j]])-dat[[j]]) / diff(range(dat[[j]]))
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
                   subset = subset.expr, 
                   extra = list(shapirowilk.p = sw.test, lambda.p = lambda.test)
  )
  
  # write
  saveRDS(models, paste0("../results/", tolower(i), ".rds"))
}








