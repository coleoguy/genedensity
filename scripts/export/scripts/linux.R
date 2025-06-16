library(MuMIn)       # multimodel inference
library(phytools)    # load and prune tree
library(caper)       # PGLS
library(foreach)
library(doParallel)
options(na.action = "na.fail")

# functions
source("functions.R")

# read data once
rsq     <- read.csv("../results/rsq.csv")
repeats <- read.csv("../results/repeat-results.csv")
tree    <- read.tree("../data/formatted-tree.nwk")

clades <- c("All", "Mammalia", "Actinopterygii", "Sauropsida")

# set up parallel backend
ncores <- parallel::detectCores() - 1
cl     <- makeCluster(ncores)
registerDoParallel(cl)

# run h = 1:1000 in parallel
res_mat <- foreach(
  h = 1:1000,
  .combine = "cbind",
  .packages = c("MuMIn", "phytools", "caper")
) %dopar% {
  run.results <- c()
  
  for (i in seq_along(clades)) {
    clade <- clades[i]
    
    # subset results
    dat <- merge(rsq, repeats, by = "species", all.x = TRUE, all.y = TRUE)
    if (clade != "All") {
      dat <- dat[dat$clade == clade, ]
    }
    
    variables <- grep("^(prop|age)\\.", names(dat), value = TRUE)
    dat <- na.omit(dat[, c("species", "clade", "rsq", variables)])
    
    # rescale
    for (j in variables) {
      dat[[j]] <- (dat[[j]] - min(dat[[j]])) / diff(range(dat[[j]]))
    }
    
    # permute block of variables
    constant.cols <- c("species", "clade", "rsq")
    perm.cols     <- setdiff(names(dat), constant.cols)
    block         <- dat[, perm.cols]
    for (l in seq_along(block)) {
      block[, l] <- sample(block[, l])
    }
    dat[, perm.cols] <- block
    
    # build terms + interactions
    reps         <- unique(sub("^[^.]*\\.", "", variables))
    interactions <- paste0("age.", reps, ":prop.", reps)
    all.terms    <- c(variables, interactions)
    
    # fit global model
    global.model <- glm(reformulate(all.terms, response = "rsq"), data = dat)
    
    # build constraint for dredge
    model.terms       <- attr(terms(global.model), "term.labels")
    model.interactions <- grep(":", model.terms, value = TRUE)
    constraints <- vapply(model.interactions, function(mi) {
      parts <- strsplit(mi, ":", fixed = TRUE)[[1]]
      sprintf("((!`%1$s`) | (`%2$s` & `%3$s`))", mi, parts[1], parts[2])
    }, character(1))
    subset.expr <- parse(text = paste(constraints, collapse = " & "))[[1]]
    
    # dredge + select top 95%
    max.terms <- if (nrow(dat) > 21) 16 else 15
    models   <- dredge(global.model,
                       subset = subset.expr,
                       m.lim   = c(0, max.terms))
    models   <- models[order(models$AICc), ]
    models   <- models[cumsum(models$weight) <= 0.95, ]
    
    # extract importance & CI
    imp    <- sw(models)
    imp    <- imp[match(all.terms, names(imp))]
    avg    <- model.avg(models)
    ci_mat <- confint(avg)[match(all.terms, rownames(confint(avg))), , drop = FALSE]
    lower  <- ci_mat[, 1]
    upper  <- ci_mat[, 2]
    
    run.results <- c(run.results, imp, upper, lower)
  }
  
  # return a named numeric vector of length (length(clades) * length(all.terms) * 3)
  run.results
}

stopCluster(cl)

# Build your final df
# First create the metadata columns once:
imp_len <- length(grep("^(prop|age)\\.", names(read.csv("../results/rsq.csv"))))*2  # approx
# (better to capture the true length from one of the parallel results:)
L        <- nrow(res_mat) / 1000          # rows per run
terms    <- names(res_mat[,1])            # if named; otherwise reconstruct
df_meta  <- data.frame(
  dataset  = rep(clades, each = length(imp)*3),
  stat     = rep(rep(c("importance","upper","lower"), each = length(imp)), times = length(clades)),
  variable = rep(rep(names(imp), 3), times = length(clades)),
  stringsAsFactors = FALSE
)

# Combine metadata with runs:
df_final <- cbind(df_meta, as.data.frame(res_mat))
colnames(df_final)[-(1:3)] <- paste0("run", 1:1000)
