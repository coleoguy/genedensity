
library(MuMIn)
library(phytools)

# functions
source("functions.R")

rsq <- read.csv("../results/rsq.csv")
repeats <- read.csv("../results/repeat-results.csv")
tree <- read.tree("../data/formatted-tree.nwk")

clade <- "Sauria"

# subset results
dat <- merge(rsq, repeats, by.x = "species", by.y = "species", all.x = T, all.y = T)
if (clade %in% c("Mammalia", "Actinopterygii", "Sauropsida")) {
  dat <- dat[dat$clade %in% clade, ]
}
variables <- colnames(dat)[grep("^(prop|age)\\.", colnames(dat))]
dat <- na.omit(dat[, c("species", "clade", "rsq", variables)])

# rescale
for (j in variables) {
  dat[[j]] <- (dat[[j]]-min(dat[[j]])) / diff(range(dat[[j]]))
}

# set up interactions
rep <- unique(sub("^[^.]*\\.", "", variables))
interactions <- paste0("age.", rep, ":prop.", rep)
all.terms <- c(
  variables, 
  interactions
)

# pgls
# cd <- comparative.data(tree, dat, names.col = "species", vcv = TRUE)


response <- "rsq"
predictors <- all.terms
n <- nrow(dat$data)
max_k <- min(n - 2, length(predictors))

model_list <- list()
model_index <- 1


for (k in 0:max_k) {
  if (k == 0) {
    fml <- as.formula(paste(response, "~ 1"))
    # fit <- pgls(fml, data = cd)
    fit <- glm(fml, data = dat)
    model_name <- paste0("M", model_index)
    model_list[[model_name]] <- fit
    model_index <- model_index + 1
    next
  }
  
  combos <- combn(predictors, k, simplify = FALSE)
  for (m in combos) {
    if (!model.marginal(m)) next
    fml <- as.formula(paste(response, "~", paste(m, collapse = " + ")))
    # fit <- try(pgls(fml, data = cd))
    fit <- glm(fml, data = dat)
    model_name <- paste0("M", model_index)
    model_list[[model_name]] <- fit
    model_index <- model_index + 1
  }
}

models <- model.sel(model_list)
models <- models[order(models$AICc), ]
models <- models[cumsum(models$weight) <= 0.95, ]

new <- c()
for (n in 1:nrow(models)) {
  l <- get.models(models, subset = n)[[1]]
  new[n] <- lambda.test(l)
}
models$lambda.p <- new
