# Zhaobo Hu
# zhaobohu2002@gmail.com

# Finds significant models and predictors of gene density homogeneity with
# model averaging
# For each response variable, do:
# (1) main analysis 
# (2) run with combined "all-repeats" predictors (prop.total, age.total)
# (3) drop top 10% species by Cook's D then re-dredge


library(MuMIn)
library(future.apply)
source("functions.R")

options(na.action = "na.fail")
options(future.globals.maxSize = 4 * 1024^3)
rsq <- read.csv("../results/rsq.csv")
repeats <- read.csv("../results/repeat-results.csv")
clades <- c("All", "Mammalia", "Actinopterygii", "Sauropsida")
responses <- c("rsq", "gini", "cv")

# do not fit interaction terms in Sauropsida; too few species relative to 
# number of models
no.interactions <- c("Sauropsida")

# the analysis in function form for future.apply
run.analysis <- function(rsq, repeats, clade, response = "rsq", 
                              variables = NULL, custom.terms = NULL, 
                              interactions = TRUE) {
  
  # merge, subset, rescale
  dat <- merge(rsq, repeats, by = "species", all = TRUE)
  if (clade != "All") {
    dat <- dat[dat$clade %in% clade, ]
  }
  
  if (is.null(variables)) {
    variables <- colnames(dat)[grep("^(prop|age)\\.", colnames(dat))]
    # get rid of the combined variables
    variables <- setdiff(variables, c("prop.total", "age.total"))
  }
  dat <- na.omit(dat[, c("species", "clade", response, variables)])
  
  for (j in variables) {
    if (diff(range(dat[[j]])) > 0) {
      dat[[j]] <- (dat[[j]] - min(dat[[j]])) / diff(range(dat[[j]]))
    }
  }

  if (is.null(custom.terms)) {
    if (interactions) {
      repeat.types <- unique(sub("^[^.]*\\.", "", variables))
      interaction.terms <- paste0("age.", repeat.types, ":prop.", repeat.types)
      all.terms <- c(variables, interaction.terms)
    } else {
      all.terms <- variables
    }
  } else {
    all.terms <- custom.terms
  }

  # get CI table from model.sel object
  models <- fitAllModels(dat, all.terms, response = response)
  ci.table <- data.frame(matrix(NA, nrow = length(all.terms), ncol = 3))
  colnames(ci.table) <- c("importance", "lower", "upper")
  rownames(ci.table) <- all.terms
  # bypass the sw() error when the confidence set collapses to one model
  if (nrow(models) < 2) { 
    top.fit <- get.models(models, 1)[[1]]
    top.terms <- attr(terms(top.fit), "term.labels")
    ci.table[intersect(top.terms, all.terms), "importance"] <- 1
    ci <- confint(top.fit)
    ci.found <- intersect(rownames(ci), all.terms)
    ci.table[ci.found, c("lower", "upper")] <- ci[ci.found, ]
  } else {
    ci.table[names(sw(models)), 1] <- sw(models)
    ci <- confint(model.avg(models))
    ci.found <- intersect(row.names(ci), all.terms)
    ci.table[ci.found, 2:3] <- ci[ci.found, ]
  }

  # build output data frame
  out <- data.frame(
    response = response, 
    clade = clade, 
    model = rownames(ci.table), 
    num.models = nrow(models), 
    estimate = sapply(seq_len(nrow(ci.table)), 
                      function(x) mean(unlist(ci.table[x, 2:3]))
                      ), 
    importance = ci.table$importance, 
    lower = ci.table$lower, 
    upper = ci.table$upper
  )

  strip_env <- function(f) { environment(f) <- baseenv(); f }
  top.formula <- strip_env(formula(get.models(models, subset = 1)[[1]]))
  
  set.formulas <- lapply(get.models(models, subset = TRUE), 
                         function(m) strip_env(formula(m))
                         )

  list(
    out = out, 
    models = models,  
    dat = dat, 
    all.terms = all.terms, 
    response = response, 
    top.formula = top.formula, 
    set.formulas = set.formulas
  )
}


plan(multisession, workers = 4)
combined.df <- data.frame()
allrepeats.df <- data.frame()
highinf.df <- data.frame()

# run the loop
for (response in responses) {
  print(paste0("Response: ", response))

  print("the main analysis")
  main.results <- future_lapply(clades, function(clade) {
    run.analysis(rsq, repeats, clade, response = response, 
                 interactions = !(clade %in% no.interactions))
  }, future.seed = TRUE)
  names(main.results) <- clades

  cur.main <- do.call(rbind, lapply(main.results, function(x) x$out))
  # cur.main <- cur.main[cur.main$importance > 0.5, ]
  combined.df <- rbind(combined.df, cur.main)


  print("standalone analysis with all repeats")
  allrepeats.results <- future_lapply(clades, function(clade) {
    run.analysis(rsq, repeats, clade, response = response, 
                      variables    = c("prop.total", "age.total"), 
                      custom.terms = c("prop.total", "age.total", 
                                       "prop.total:age.total")
                      )
  }, future.seed = TRUE)
  names(allrepeats.results) <- clades

  cur.ar <- do.call(rbind, lapply(allrepeats.results, function(x) x$out))
  allrepeats.df <- rbind(allrepeats.df, cur.ar)

  print("drop 10% top species and refit")
  highinf.results <- future_lapply(clades, function(clade) {
    dat <- main.results[[clade]]$dat
    top.formula <- main.results[[clade]]$top.formula

    top.model <- lm(top.formula, data = dat)
    cd <- cooks.distance(top.model)

    n.drop <- ceiling(0.10 * nrow(dat))
    drop.idx <- order(cd, decreasing = TRUE)[1:n.drop]
    dropped.species <- dat$species[drop.idx]

    rsq.filtered <- rsq[!rsq$species %in% dropped.species, ]
    result <- run.analysis(rsq.filtered, repeats, clade, response = response, 
                           interactions = !(clade %in% no.interactions))

    if (nrow(result$out) > 0) {
      result$out$num.dropped <- n.drop
      result$out$dropped.species <- paste(dropped.species, collapse = ";")
    }
    return(result$out)
  }, future.seed = TRUE)

  cur.hi <- do.call(rbind, highinf.results)
  highinf.df <- rbind(highinf.df, cur.hi)
}

write.csv(combined.df, "../results/model-averaging.csv", row.names = FALSE)
write.csv(allrepeats.df, "../results/model-averaging-allrepeats.csv", row.names = FALSE)
write.csv(highinf.df, "../results/model-averaging-highinfluence.csv", row.names = FALSE)