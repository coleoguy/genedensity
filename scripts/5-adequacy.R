# Zhaobo Hu
# zhaobohu2002@gmail.com

# A permutation-based approach to assess model adequacy.

library(MuMIn)
library(phytools)
options(na.action = "na.fail")

source("functions.R")

rsq <- read.csv("../results/rsq.csv")
repeats <- read.csv("../results/repeat-results.csv")
responses <- c("rsq", "gini", "cv")
clades <- c("All", "Mammalia", "Actinopterygii", "Sauropsida")
constant.cols <- c("species", "clade", responses)

main <- merge(rsq, repeats, by = "species", all = TRUE)
variables <- colnames(main)[grep("^(prop|age)\\.", colnames(main))]
variables <- setdiff(variables, c("prop.total", "age.total"))
main <- na.omit(main[, c(constant.cols, variables)])

rep.types <- unique(sub("^[^.]*\\.", "", variables))
all.terms <- c(variables, paste0("age.", rep.types, ":prop.", rep.types))

df <- data.frame(
  response = rep(responses, each = length(clades) * 3 * length(all.terms)), 
  dataset = rep(rep(clades, each = 3 * length(all.terms)), length(responses)), 
  stat = rep(
    rep(c("importance", "upper", "lower"), each = length(all.terms)), 
    length(responses) * length(clades)
  ), 
  variable = rep(all.terms, length(responses) * length(clades) * 3)
)

n.runs <- 2000

for (i in seq_len(n.runs)) {
  
  perm.main <- main
  perm.cols <- setdiff(names(perm.main), constant.cols)
  perm.main[, perm.cols] <- perm.main[sample(nrow(perm.main)), perm.cols]
  
  run.results <- c()
  for (cur.response in responses) {
    for (clade in clades) {
      # subset to clade
      dat <- if (clade == "All") {
        perm.main 
        } else perm.main[perm.main$clade == clade, ]
      
      # build term list for this clade
      cur.vars <- colnames(dat)[grep("^(prop|age)\\.", colnames(dat))]
      cur.rep.types <- unique(sub("^[^.]*\\.", "", cur.vars))
      cur.interact <- paste0("age.", cur.rep.types, ":prop.", cur.rep.types)
      cur.terms <- c(cur.vars, cur.interact)
      max.vars <- nrow(dat) - 2
      # rescae predictors to [0, 1]
      for (j in cur.vars) {
        rng <- diff(range(dat[[j]]))
        if (rng > 0) dat[[j]] <- (dat[[j]] - min(dat[[j]])) / rng
      }
      
      # fit all valid models
      model.list <- list()
      for (k in seq_along(cur.terms)) {
        for (m in combn(cur.terms, k, simplify = FALSE)) {
          fml <- as.formula(paste(cur.response, "~", paste(m, collapse = " + ")))
          if (length(attr(terms(fml), "term.labels")) >= max.vars) next
          if (!model.marginal(m)) next
          model.list[[length(model.list) + 1L]] <- lm(fml, data = dat)
        }
      }
      names(model.list) <- paste0("M", seq_along(model.list))
      
      # keep models within the 95% cumulative weight confidence set
      models <- model.sel(model.list)
      models <- models[cumsum(models$weight[order(models$AICc)]) <= 0.95, ]
      
      # extract importance and model-averaged CIs
      ct <- data.frame(matrix(NA, nrow = length(all.terms), ncol = 3))
      colnames(ct) <- c("importance", "lower", "upper")
      rownames(ct) <- all.terms
      
      if (nrow(models) < 2) {
        # single-model confidence set: sw() and model.avg() would error
        top.fit <- get.models(models, 1)[[1]]
        top.terms <- attr(terms(top.fit), "term.labels")
        ct[intersect(top.terms, all.terms), "importance"] <- 1
        ci <- confint(top.fit)
        ci.found <- intersect(rownames(ci), all.terms)
        ct[ci.found, c("lower", "upper")] <- ci[ci.found, ]
      } else {
        ct[names(sw(models)), "importance"] <- sw(models)
        ci <- confint(model.avg(models))
        ci.found <- intersect(rownames(ci), all.terms)
        ct[ci.found, c("lower", "upper")] <- ci[ci.found, ]
      }
      
      run.results <- c(run.results, unlist(ct))
    }
  }
  df[[paste0("run", i)]] <- run.results
  cat("run", i, "of", n.runs, "complete\n")
}

write.csv(df, "permute.csv", row.names = FALSE)
