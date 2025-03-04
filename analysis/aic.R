

library(MuMIn) # multimodel inference
library(phytools) # load and prune tree
library(caper) # PGLS
library(lmtest) # Breusch-Pagan test

# loop for each clade
for (i in c("All", "Mammalia", "Actinopterygii", "Sauria")) {
  
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
  
  
  # subset clades
  if (i %in% c("Mammalia", "Actinopterygii", "Sauria")) {
    dat <- dat[dat$clade %in% i, ]
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
    
    # determine repeats present in model
    rep <- terms[which(pgls.models[k, terms] == 1)]
    prop <- rowSums(dat[, paste0("rep.prop.", rep)])
    age <- dat[, paste0("rep.age.", paste(terms[which(pgls.models[k, terms] == 1)], collapse = "."))]
    
    # normalize
    sub <- dat[, c("species", "clade", "rsq")]
    sub$age.norm <- (age - range(age)[1]) / diff(range(age))
    sub$prop.norm <- (prop - range(prop)[1]) / diff(range(prop))
    
    # make models
    gls.model <- glm(formula, data = sub)
    cd <- comparative.data(pruned.tree, data = sub, names.col = "species", vcv = TRUE)
    pgls.model <- pgls(formula, data = cd)
    
    # define terms and formula
    all.terms <- c("age.norm", "prop.norm", "age.norm:prop.norm")
    model.terms <- as.data.frame(pgls.models)[k, all.terms]
    if (all(is.na(model.terms))) {
      model.terms <- "1"
    } else {
      model.terms <- all.terms[which(!is.na(model.terms))]
    }
    formula <- reformulate(model.terms, response = "rsq")
    
    # lambda (phylogenetic signal)
    gls.res <- setNames(resid(gls.model), sub$species)
    pgls.models$lambda.p[k] <- phylosig(tree = pruned.tree, 
                                        x = gls.res, 
                                        method = "lambda", 
                                        test = TRUE, 
                                        niter = 1000)$P
    phylosig(pruned.tree, gls.res)
    
    # shapiro-welk (residual normality)
    pgls.res <- resid(pgls.model)
    pgls.models$sw.p[k] <- shapiro.test(pgls.res)$p.value
    
    # breusch-pagan (homoskedasticity)
    if (length(model.terms) == 1 & all(model.terms == "1")) {
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


library(MuMIn)
dat <- readRDS("../results/all.8.rds")
dat <- dat[!duplicated(dat)]
assignInNamespace(
  ".modelNames",
  function(allTerms, uqTerms) {
    paste0("model", seq_along(allTerms))
  },
  ns = "MuMIn"
)
avg <- model.avg(dat, subset = delta <= 10)




int <- confint(avg, full = FALSE)
estimates <- coef(avg)
terms <- rownames(int)
y <- seq_along(estimates)
plot(y = estimates, x = y, ylim = range(int), pch = 19, 
     ylab = "Estimate", xlab = NA, axes = FALSE,
     main = "Coefficient estimates for model using all species")
segments(y, int[, 1], y, int[, 2], lwd = 2)
axis(2)
axis(1, at = y, labels = terms, las = 2)
abline(h = 0, lty = 2, col = "gray")
box()
sw(avg)

















############### pseudo residual and qq


library(MuMIn) # multimodel inference
library(phytools) # load and prune tree
library(caper) # PGLS
library(lmtest) # Breusch-Pagan test

# loop for each clade
for (i in c("All")) {
  
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
  
  
  # subset clades
  if (i %in% c("Mammalia", "Actinopterygii", "Sauria")) {
    dat <- dat[dat$clade %in% i, ]
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
    
    # determine repeats present in model
    rep <- terms[which(pgls.models[k, terms] == 1)]
    prop <- rowSums(dat[, paste0("rep.prop.", rep)])
    age <- dat[, paste0("rep.age.", paste(terms[which(pgls.models[k, terms] == 1)], collapse = "."))]
    
    # normalize
    sub <- dat[, c("species", "clade", "rsq")]
    sub$age.norm <- (age - range(age)[1]) / diff(range(age))
    sub$prop.norm <- (prop - range(prop)[1]) / diff(range(prop))
    
    # make models
    gls.model <- glm(formula, data = sub)
    cd <- comparative.data(pruned.tree, data = sub, names.col = "species", vcv = TRUE)
    pgls.model <- pgls(formula, data = cd)
    
    # define terms and formula
    all.terms <- c("age.norm", "prop.norm", "age.norm:prop.norm")
    model.terms <- as.data.frame(pgls.models)[k, all.terms]
    if (all(is.na(model.terms))) {
      model.terms <- "1"
    } else {
      model.terms <- all.terms[which(!is.na(model.terms))]
    }
    formula <- reformulate(model.terms, response = "rsq")
    
    # lambda (phylogenetic signal)
    gls.res <- setNames(resid(gls.model), sub$species)
    pgls.models$lambda.p[k] <- phylosig(tree = pruned.tree, 
                                        x = gls.res, 
                                        method = "lambda", 
                                        test = TRUE, 
                                        niter = 1000)$P
    phylosig(pruned.tree, gls.res)
    
    # shapiro-welk (residual normality)
    pgls.res <- resid(pgls.model)
    pgls.models$sw.p[k] <- shapiro.test(pgls.res)$p.value
    
    # breusch-pagan (homoskedasticity)
    if (length(model.terms) == 1 & all(model.terms == "1")) {
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






################# snippets





library(MuMIn)


# load results and filter for terms of interest
dat <- readRDS("../results/actinopterygii.8.rds")
# terms <- read.csv("../results/parsed.csv")
# terms <- gsub("rep.prop.", "", names(terms)[grep("^rep.prop", names(terms))])
# dat <- dat[c(which(rowSums(dat[, terms]) == length(terms)), 
             #which(rowSums(dat[, c("others", "unknown")]) == 0)), ]

# filter for models that meet assumptions
# dat <- dat[which(dat$lambda.p < 0.05)]
# dat <- dat[which(dat$bp.p < 0.05)]
# dat <- dat[which(dat$sw.p < 0.05)]

# subset and average
# dat <- subset(dat, delta < 9.0)
# dat <- dat[!duplicated(dat), ]

# recalculate
# dat <- rbind(dat, dat[1,])
# dat <- dat[-1,]

rownames(dat) <- make.names(seq_len(nrow(dat)), unique = TRUE)
dat <- dat[, c(1:4, 15:19)]
avg <- model.avg(dat, subset = T)
summary(avg)
confint(avg, full = F) 

avg$coefficients
avg$coefArray



# Ensure row names are unique
newNames <- paste0("model", seq_len(nrow(dat)))

# Assign new row names
rownames(dat) <- newNames

# Modify the internal "model.names" attribute if it exists
if (!is.null(attr(dat, "model.names"))) {
  attr(dat, "model.names") <- newNames
}

assignInNamespace(
  ".modelNames",
  function(allTerms, uqTerms) {
    paste0("model", seq_along(allTerms))
  },
  ns = "MuMIn"
)

# Now try model averaging
avg <- model.avg(dat, subset = TRUE)













##############3333 pseudo-residuals

# Assuming 'avgmod' is your model.avg object and 'original_data' is your dataset

# Generate predictions
predictions <- predict(avg, newdata = cd$data[3:4])

# Calculate pseudo-residuals
pseudo_residuals <- original_data$response_variable - predictions

# Create QQ plot
qqnorm(pseudo_residuals)
qqline(pseudo_residuals)




l <- str(dat)


