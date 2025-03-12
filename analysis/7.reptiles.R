

number <- 1

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
  lambda.p <- phylosig(pruned.tree, res, method = "lambda", test = TRUE, niter = 10)$P
  return(lambda.p)
}

# subset results
tree <- read.tree("../data/formatted.tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
dat <- read.csv("../results/parsed.csv")
dat <- dat[!is.na(dat$chromnum.1n) & !duplicated(dat$species), ]
dat <- dat[dat$clade %in% "Sauria", ]
int <- intersect(dat$species, tree$tip.label)
dat <- dat[dat$species %in% int, ]
pruned.tree <- keep.tip(tree, int)
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
cd <- comparative.data(pruned.tree, dat, names.col = "species", vcv = TRUE)

# set up interactions
rep <- unique(sub("^[^.]*\\.", "", variables))
interactions <- paste0("age.", rep, ":prop.", rep)
all.terms <- c(
  variables, 
  interactions
)
all.terms <- all.terms[!grepl(all.terms[number], all.terms)]


# global model
global.model <- pgls(reformulate(all.terms, response = "rsq"), data = cd)

# set constraints
model.terms <- unlist(strsplit(as.character(global.model$formula)[3], " \\+ "))
model.interactions <- grep(":", all.terms, value = TRUE)
constraints <- character(length(model.interactions))
for (k in seq_along(model.interactions)) {
  parts <- strsplit(model.interactions[k], ":")[[1]]
  constraints[k] <- sprintf("((!`%s`) | (%s & %s))", model.interactions[k], parts[1], parts[2])
}
subset.expr <- parse(text = paste(constraints, collapse = " & "))[[1]]

# dredge
models <- dredge(global.model, 
                 subset = subset.expr, 
                 extra = list(shapirowilk.p = sw.test, lambda.p = lambda.test),
                 m.lim = c(1,16)
)

# write
saveRDS(models, paste0("../results/rep.", i, ".models.rds"))

