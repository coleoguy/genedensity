

library(MuMIn) # multimodel inference
library(phytools) # load and prune tree
library(caper) # PGLS
library(lmtest) # Breusch-Pagan test

# subset results
tree <- read.tree("../data/formatted.tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
dat <- read.csv("../results/parsed.csv")
dat <- dat[!is.na(dat$chromnum.1n) & !duplicated(dat$species), ]
terms <- gsub("rep.prop.", "", names(dat)[grep("^rep.prop", names(dat))])
combs <- gsub("rep.age.", "", names(dat)[grep("^rep.age", names(dat))])
dat <- na.omit(dat[, c("species", 
                       "clade", 
                       "rsq", 
                       paste0("rep.prop.", terms), 
                       paste0("rep.age.", combs))])

# loop for each clade
for (i in c("All", "Mammalia", "Actinopterygii", "Sauria")) {
  
  # subset clades
  if (i %in% c("Mammalia", "Actinopterygii", "Sauria")) {
    dat <- dat[dat$clade %in% c("Mammalia"), ]
  }
  
  # prune tree and dataset
  int <- intersect(dat$species, tree$tip.label)
  dat <- dat[dat$species %in% int, ]
  pruned.tree <- keep.tip(tree, int)
  
  # loop through combinations of repeats; fit PGLS model for each
  pgls.models <- other.stats <- data.frame()
  for (j in combs) {
    
    # get repeat age and proportion
    rep <- unlist(strsplit(j, "\\."))
    prop <- as.numeric(rowSums(as.data.frame(dat[, c(paste0("rep.prop.", rep))])))
    age <- dat[, c(paste0("rep.age.", j))]
    
    # create new dataframe to make comparative data object for PGLS; normalize data
    sub <- dat[, c("species", "clade", "rsq")]
    sub$age.norm <- (age - range(age)[1]) / diff(range(age))
    sub$prop.norm <- (prop - range(prop)[1]) / diff(range(prop))
    
    # create list with presence/absence for use in "extra" argument in dredge()
    is.present <- as.numeric(terms %in% rep)
    extra.list <- eval(parse(text = paste0(
      "alist(",
      paste(terms, 
            "= function(x) is.present[", 
            seq_along(terms), 
            "]", 
            collapse = ", "), 
            ")")
      ))
    
    # add empty columns to extra.list
    na.cols <- alist(
      clade = function(x) NA,
      sw.p = function(x) NA,
      bp.p = function(x) NA,
      lambda.p = function(x) NA
    )
    extra.list <- c(extra.list, na.cols)
    
    # fit model
    cd <- comparative.data(pruned.tree, sub, names.col = "species", vcv = TRUE)
    models <- dredge(pgls(rsq ~ age.norm * prop.norm, data = cd), 
                     subset = dc(x1, x2, x1:x2), 
                     extra = extra.list)
    
    # record estimates
    pgls.models <- rbind(pgls.models, models)
  }
  
  # for each model do: shapiro-welk, breusch-pagen, pagel's lambda
  for (k in 1:nrow(pgls.models)) {
    
    # define terms and formula
    all.terms <- c("age.norm", "prop.norm", "age.norm:prop.norm")
    model.terms <- as.data.frame(pgls.models)[k, all.terms]
    if (all(is.na(model.terms))) {
      model.terms <- "1"
    } else {
      model.terms <- all.terms[which(!is.na(model.terms))]
    }
    
    # lambda (phylogenetic signal)
    formula <- reformulate(model.terms, response = "rsq")
    gls.model <- glm(formula, data = cd$data)
    gls.res <- resid(gls.model)
    pgls.models$lambda.p[k] <- phylosig(tree = pruned.tree, 
                                        x = gls.res, 
                                        method = "lambda", 
                                        test = TRUE, 
                                        niter = 300)$P
    
    # shapiro-welk (residual normality)
    pgls.model <- get.models(pgls.models, k)[[1]]
    pgls.res <- resid(pgls.model)
    pgls.models$sw.p[k] <- shapiro.test(pgls.res)$p.value
    
    # breusch-pagan (homoskedasticity)
    if (all(is.na(model.terms))) {
      pgls.models$bp.p[k] <- NA
    } else {
      pgls.models$bp.p[k] <- bptest(pgls.model)$p.value[[1]]
    }
  }
  
  # convert presence/absence to true/false
  for (term in terms) {
    pgls.models[[term]] <- as.logical(unlist(pgls.models[[term]]))
  }
  
  # set clade
  pgls.models$clade <- i
  
  # write
  saveRDS(pgls.models, paste0("../results/", tolower(i), ".models.rds"))
}











































# snippets

confint(model.avg(pgls.models))

# remove models involving misc/unidentified repeats; re-calculate deltas
pgls.models <- pgls.models[
  apply(pgls.models[, terms], 1, function(x) all(x == TRUE)) |  # All columns are TRUE
    (pgls.models[, terms[length(terms)-1]] == FALSE & pgls.models[, terms[length(terms)]] == FALSE),  # Last two columns are FALSE
  , ]
dummy <- pgls.models[1,]
pgls.models <- rbind(pgls.models, dummy)
pgls.models <- pgls.models[-1,]