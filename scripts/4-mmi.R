# Zhaobo Hu
# zhaobohu2002@gmail.com

# Finds significant models and predictors of gene density homogeneity with
# model averaging
# For each response variable, do:
# (1) main analysis 
# (2) run with combined "all-repeats" predictors (prop.total, age.total)
# (3) refit the top model after dropping each species
# (4) refit the original 95% confidence set with a random 80% subsample (1000x)
# (5) drop top 10% species by Cook's D then re-dredge


library(MuMIn)
library(phytools)
library(caper)
library(future.apply)
source("functions.R")

options(na.action = "na.fail")
options(future.globals.maxSize = 4 * 1024^3)

rsq <- read.csv("../results/rsq.csv")
repeats <- read.csv("../results/repeat-results.csv")
tree <- read.tree("../data/formatted-tree.nwk")

clades <- c("All", "Mammalia", "Actinopterygii", "Sauropsida")
responses <- c("rsq", "gini", "cv")


# get CI table from model.sel object
get.table <- function(models, all.terms) {
  ci.table <- data.frame(matrix(NA, nrow = length(all.terms), ncol = 3))
  colnames(ci.table) <- c("importance", "lower", "upper")
  rownames(ci.table) <- all.terms
  ci.table[names(sw(models)), 1] <- sw(models)
  ci <- confint(model.avg(models))
  ci.found <- intersect(row.names(ci), all.terms)
  ci.table[ci.found, 2:3] <- ci[ci.found, ]
  return(ci.table)
}

# run the analysis for one clade and one response
run.analysis <- function(rsq, repeats, clade, response = "rsq",
                              variables = NULL, custom.terms = NULL) {
  
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
    repeat.types <- unique(sub("^[^.]*\\.", "", variables))
    interactions <- paste0("age.", repeat.types, ":prop.", repeat.types)
    all.terms <- c(variables, interactions)
  } else {
    all.terms <- custom.terms
  }

  models   <- fit_all_models(dat, all.terms, response = response)
  ci.table <- get.table(models, all.terms)

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


plan(multisession, workers = 12)
combined.df <- data.frame()
allrepeats.df <- data.frame()
#loso.df <- data.frame()
#subsample.df <- data.frame()
highinf.df <- data.frame()

# run the loop
for (response in responses) {
  print(paste0("Response: ", response))

  print("1. main analysis")
  main.results <- future_lapply(clades, function(clade) {
    run.analysis(rsq, repeats, clade, response = response)
  }, future.seed = TRUE)
  names(main.results) <- clades

  cur.main <- do.call(rbind, lapply(main.results, function(x) x$out))
  # cur.main <- cur.main[cur.main$importance > 0.5, ]
  combined.df <- rbind(combined.df, cur.main)


  print("2. all repeats")
  allrepeats.results <- future_lapply(clades, function(clade) {
    run.analysis(rsq, repeats, clade, response = response, 
                      variables    = c("prop.total", "age.total"), 
                      custom.terms = c("prop.total", "age.total", 
                                       "prop.total:age.total")
                      )
  }, future.seed = TRUE)
  names(allrepeats.results) <- clades

  cur.ar <- do.call(rbind, lapply(allrepeats.results, function(x) x$out))
  # cur.ar <- cur.ar[cur.ar$importance > 0.5, ]
  allrepeats.df <- rbind(allrepeats.df, cur.ar)

  #print("3. leaving out each species")
  #for (clade in clades) {
  #  dat <- main.results[[clade]]$dat
  #  top.formula <- main.results[[clade]]$top.formula

  #  p <- progressor(steps = nrow(dat))
  #  loso.list <- future_lapply(seq_len(nrow(dat)), function(i) {
  #    fit <- glm(top.formula, data = dat[-i, ])
  #    coefs <- coef(fit)
  #    p(sprintf("[%s] dropped %s", clade, dat$species[i]))
  #    data.frame(
  #      response = response, 
  #      clade = clade, 
  #      dropped.species = dat$species[i], 
  #      predictor = names(coefs), 
  #      coef = coefs, 
  #      row.names = NULL
  #    )
  #  }, future.seed = TRUE)

  #  loso.df <- rbind(loso.df, do.call(rbind, loso.list))
  #}


  #print("4. drop 20% of species and refit 95% of models")
  #n.iter <- 1000
  #for (clade in clades) {
  #  dat <- main.results[[clade]]$dat
  #  set.formulas <- main.results[[clade]]$set.formulas
  #  all.terms <- main.results[[clade]]$all.terms

  #  sub.list <- future_lapply(seq_len(n.iter), function(b) {
  #    idx <- sample(seq_len(nrow(dat)), size = round(0.8 * nrow(dat)))
  #    sub.dat <- dat[idx, ]

  #    res <- tryCatch({
  #      fits <- lapply(set.formulas, function(f) glm(f, data = sub.dat))
  #      names(fits) <- paste0("M", seq_along(fits))
  #      models.sub <- model.sel(fits)
  #      ci.table <- get.table(models.sub, all.terms)
  #      ci.table$response <- response
  #      ci.table$clade <- clade
  #      ci.table$iter <- b
  #      ci.table$predictor <- rownames(ci.table)
  #      rownames(ci.table) <- NULL
  #      ci.table
  #    }, error = function(e) NULL)
  #    return(res)
  #  }, future.seed = TRUE)

  #  subsample.df <- rbind(subsample.df, do.call(rbind, sub.list))
  #}


  print("5. drop 10% top species and refit")
  highinf.results <- future_lapply(clades, function(clade) {
    dat <- main.results[[clade]]$dat
    top.formula <- main.results[[clade]]$top.formula

    top.model <- glm(top.formula, data = dat)
    cd <- cooks.distance(top.model)

    n.drop <- ceiling(0.10 * nrow(dat))
    drop.idx <- order(cd, decreasing = TRUE)[1:n.drop]
    dropped.species <- dat$species[drop.idx]

    rsq.filtered <- rsq[!rsq$species %in% dropped.species, ]
    result <- run.analysis(rsq.filtered, repeats, clade, response = response)

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
#write.csv(loso.df, "../results/loso.csv", row.names = FALSE)
#write.csv(subsample.df, "../results/subsample.csv", row.names = FALSE)
write.csv(highinf.df, "../results/model-averaging-highinfluence.csv", row.names = FALSE)

print("all done")
