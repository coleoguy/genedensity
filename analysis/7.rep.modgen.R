

library(MuMIn) # multimodel inference
library(phytools) # load and prune tree
library(caper) # PGLS
options(na.action = "na.fail")  

# functions
sw.test <- function(model) {
  res <- residuals(model)
  sw.p <- shapiro.test(res)$p.value
  return(sw.p)
}
lambda.test <- function(model) {
  cur.terms <- unlist(strsplit(as.character(model$formula)[3], " \\+ "))
  nophylo.formula <- reformulate(cur.terms, response = "rsq")
  nophylo.model <- glm(nophylo.formula, data = cd$data)
  res <- residuals(nophylo.model)
  lambda.p <- phylosig(pruned.tree, res, method = "lambda", test = TRUE, niter = 100)$P
  return(lambda.p)
}


blocks <- list(
  dna     = c("age.dna",     "prop.dna",     "age.dna:prop.dna"),
  line    = c("age.line",    "prop.line",    "age.line:prop.line"),
  ltr     = c("age.ltr",     "prop.ltr",     "age.ltr:prop.ltr"),
  sine    = c("age.sine",    "prop.sine",    "age.sine:prop.sine"),
  unknown = c("age.unknown", "prop.unknown", "age.unknown:prop.unknown"),
  others  = c("age.others",  "prop.others",  "age.others:prop.others")
)

term.list <- list()

# remove two interactions
to.drop <- combn(names(blocks), 2, simplify = FALSE)
for (i in to.drop) {
  model <- c()
  for (rep in names(blocks)) {
    if (rep %in% i) {
      for (j in 1:2) {
        model <- c(model, blocks[[rep]][j])
      }
    } else {
      for (j in 1:3) {
        model <- c(model, blocks[[rep]][j])
      }
    }
  }
  term.list[[length(term.list) + 1]] <- model
}

# remove one interaction and one main effect 
for (rep in names(blocks)) {
  for (choice in 1:2) {
    model <- c()
    for (rep2 in names(blocks)) {
      if (rep2 == rep) {
        if (choice == 1) {
          model <- c(model, blocks[[rep2]][2])
        } else {
          model <- c(model, blocks[[rep2]][1])
        }
      } else {
        for (j in 1:3) {
          model <- c(model, blocks[[rep2]][j])
        }
      }
    }
    term.list[[length(term.list) + 1]] <- model
  }
}



for (i in 1:length(term.list)) {
  
  # subset results
  tree <- read.tree("../data/formatted.tree.nwk")
  tree$tip.label <- gsub("_", " ", tree$tip.label)
  dat <- read.csv("../results/parsed.csv")
  dat <- dat[!is.na(dat$chromnum.1n) & !duplicated(dat$species), ]
  dat <- dat[dat$clade %in% "Sauria", ]
  int <- intersect(dat$species, tree$tip.label)
  dat <- dat[dat$species %in% int, ]
  pruned.tree <- keep.tip(tree, int)
  terms <- term.list[[i]]
  main <- terms[!grepl(":", terms)]
  interactions <- terms[grepl(":", terms)]
  dat <- na.omit(dat[, c("species", 
                         "clade", 
                         "rsq", 
                         main
  )])
  
  # normalize
  for (j in main) {
    dat[[j]] <- (max(dat[[j]])-dat[[j]]) / diff(range(dat[[j]]))
  }
  cd <- comparative.data(pruned.tree, dat, names.col = "species", vcv = TRUE)
  
  # global model
  global.model <<- pgls(reformulate(terms, response = "rsq"), data = cd)
  
  # set constraints
  constraints <- character(length(interactions))
  for (k in seq_along(interactions)) {
    parts <- strsplit(interactions[k], ":")[[1]]
    constraints[k] <- sprintf("((!`%s`) | (%s & %s))", interactions[k], parts[1], parts[2])
  }
  subset.expr <- parse(text = paste(constraints, collapse = " & "))[[1]]
  
  # dredge
  models <- dredge(global.model, 
                   subset = subset.expr, 
                   extra = list(shapirowilk.p = sw.test, lambda.p = lambda.test)
  )
  
  # write
  saveRDS(models, paste0("../results/Sauria.", i, ".rds"))
}



