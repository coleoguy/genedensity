


library(MuMIn)
library(phytools)
library(caper)
library(parallel)
options(na.action = "na.fail")  

# function to apply shapiro-wilk test to model residuals and return the p value
sw.test <- function(model) {
  res <- residuals(model)
  shapiro.test(res)$p.value
}

# function to apply pagel's lambda test to model residuals and return the p value
lambda.test <- function(model) {
  if(model$df.residual < 2) {
    return(NA)
  }
  res <- setNames(residuals(model), dat$species)
  phylosig(tree, res, method = "lambda", test = TRUE, niter = 100)$P
}

# function to map var count to var combinations within repeat types
map.varnum <- function(block, state) {
  if (state == "no.vars") {
    list(terms = character(0), count = 0)
  } else if (state == "age") {
    list(terms = block[1], count = 1)
  } else if (state == "prop") {
    list(terms = block[2], count = 1)
  } else if (state == "both") {
    list(terms = block[1:2], count = 2)
  } else if (state == "interac") {
    list(terms = block, count = 3)
  }
}

# function to dredge each global model
run <- function(candidate_info, candidate_index) {
  terms <- candidate_info$predictors
  main <- terms[!grepl(":", terms)]
  interactions <- terms[grepl(":", terms)]
  
  subse <- na.omit(dat[, c("species", "clade", "rsq", main)]) # subset data
  mod <- glm(reformulate(terms, response = "rsq"), data = subse) # fit model
  
  # prohibit interactions without [resence of both main effects
  if(length(interactions) > 0) {
    constraints <- sapply(interactions, function(int) {
      parts <- strsplit(int, ":")[[1]]
      sprintf("((!`%s`) | (`%s` & `%s`))", int, parts[1], parts[2])
    })
    subset.expr <- parse(text = paste(constraints, collapse = " & "))[[1]]
  } else {
    subset.expr <- TRUE
  }
  
  # run dredge
  dredged <- dredge(mod, 
                    subset = subset.expr, 
                    extra = list(shapirowilk.p = sw.test,
                                 lambda.p = lambda.test))
  # save
  saveRDS(dredged, paste0("../../results/sauria.09/Sauria.", candidate_index, ".rds"))
  return(candidate_index)
}

# group by repeat type
blocks <- list(
  dna     = c("age.dna", "prop.dna", "age.dna:prop.dna"),
  line    = c("age.line", "prop.line", "age.line:prop.line"),
  ltr     = c("age.ltr", "prop.ltr", "age.ltr:prop.ltr"),
  sine    = c("age.sine", "prop.sine", "age.sine:prop.sine"),
  unknown = c("age.unknown", "prop.unknown", "age.unknown:prop.unknown"),
  others  = c("age.others", "prop.others", "age.others:prop.others")
)

# build a grid with every possible model
states <- c("no.vars", "age", "prop", "both", "interac")
grid <- expand.grid(rep(list(states), length(blocks)), stringsAsFactors = FALSE)
colnames(grid) <- names(blocks)

# from all models, find ones with 12 terms
model.list <- list()
for (i in 1:nrow(grid)) {
  model.terms <- character(0)
  total.varcount <- 0
  for (b in names(blocks)) {
    mapping <- map.varnum(blocks[[b]], grid[i, b])
    model.terms <- c(model.terms, mapping$terms)
    total.varcount <- total.varcount + mapping$count
  }
  if (total.varcount == 12) {
    model.list[[length(model.list) + 1]] <- list(
      predictors = model.terms,
      num.vars = total.varcount
    )
  }
}


# read data and prune tree
tree <- read.tree("../../data/formatted.tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
dat <- read.csv("../../results/parsed.09.csv")
dat <- dat[!is.na(dat$chromnum.1n) & !duplicated(dat$species), ]
dat <- dat[dat$clade %in% "Sauria", ]
dat <- dat[dat$species != "Strigops habroptila", ] # remove problematic species
int <- intersect(dat$species, tree$tip.label)
pruned.tree <- keep.tip(tree, int)

# normalize data
all.main <- c("age.dna", "prop.dna", "age.line", "prop.line", "age.ltr", 
              "prop.ltr", "age.sine", "prop.sine", "age.unknown", 
              "prop.unknown", "age.others", "prop.others")
for (var in all.main) {
  if (var %in% names(dat)) {
    dat[[var]] <- (max(dat[[var]]) - dat[[var]]) / diff(range(dat[[var]]))
  }
}

ncores <- round(0.7 * detectCores())
cl <- makeCluster(ncores)

# export libraries and variables
clusterEvalQ(cl, {
  library(MuMIn)
  library(phytools)
  library(caper)
})
clusterExport(cl, 
              varlist = c("model.list", "dat", "pruned.tree", "sw.test", 
                          "lambda.test", "run", "blocks"),
              envir = environment())

# run
results <- parLapply(cl, seq_along(model.list), function(i) {
  run(model.list[[i]], i)
})
stopCluster(cl)


