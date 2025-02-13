


# for each repeat in each clade in each threshold, fit a glm model where the 
# response is r2 and the predictors are (divergence containing 1/2 repeat sum), 
# repeat proportion, and the interaction of these two terms. then, select
# for the best model using step(). if the residuals of this model show 
# phylogenetic signal, fit a PGLS model instead and select the best model using
# AIC (step() doesn't work for PGLS objects)

# helper functions
# universal ordering of a table of 5 models
sortModels <- function(x) {
  first <- which(!is.na(x[, terms[-1]][1]) & 
                   !is.na(x[, terms[-1]][2]) & 
                   !is.na(x[, terms[-1]][3])) 
  second <- which(!is.na(x[, terms[-1]][1]) & 
                    !is.na(x[, terms[-1]][2]) & 
                    is.na(x[, terms[-1]][3])) 
  third <- which(is.na(x[, terms[-1]][1]) & 
                   !is.na(x[, terms[-1]][2]) & 
                   is.na(x[, terms[-1]][3])) 
  fourth <- which(!is.na(x[, terms[-1]][1]) & 
                    is.na(x[, terms[-1]][2]) & 
                    is.na(x[, terms[-1]][3])) 
  fifth <- which(is.na(x[, terms[-1]][1]) & 
                   is.na(x[, terms[-1]][2]) & 
                   is.na(x[, terms[-1]][3])) 
  order <- c(first, second, third, fourth, fifth)
  return(x[order, ])
}



library(caper) # apply PGLS
library(data.table) # quickly read data
library(MuMIn) # calculate AIC
library(piecewiseSEM) # calculate R2 for PGLS objects
library(phytools) # phy stuff
terms <- c(
  "(Intercept)", 
  "median.trans", 
  "rep.prop", 
  "median.trans:rep.prop"
)
dat <- fread("../results/parsed.csv")
dat <- as.data.frame(dat)
dat <- dat[!is.na(dat$chromnum.1n), ]
dat <- dat[dat$thrs == 0.8, ]

# store results a list to save memory
lis <- list()

# for each clade
for (cl in c("Total", "Mammalia", "Actinopterygii", "Sauria")) {
  # for each repeat type
  for (rep in c("total", "line", "sine", "ltr", "dna", "rc")) {
    head <- as.data.frame(t(setNames(c(cl, rep), c("cl", "rep"))))
    # subset data for normalization
    sub <- dat[!duplicated(dat$species), ]
    sub <- na.omit(sub[, c("species", "rsq", "clade", paste0(rep, ".rep.median"), paste0(rep, ".rep.pct"))])
    if (cl %in% c("Mammalia", "Actinopterygii", "Sauria")) {
      sub <- sub[sub$clade == cl, ]
    }
    sub <- sub[sub$species != "Callithrix jacchus", ]
    if (nrow(sub) <= 2) {
      next
    }
    # transform and normalize
    median <- 1 - (sub[[paste0(rep, ".rep.median")]]/70)
    median.range <- range(na.omit(median))
    sub$median.trans <- (median - median.range[1]) / diff(median.range)
    vol <- sub[[paste0(rep, ".rep.pct")]] / 100
    vol.range <- range(na.omit(vol))
    sub$rep.prop <- (vol - vol.range[1]) / diff(vol.range)
    sub <- na.omit(sub[, c("species", "rsq", "clade", "median.trans", "rep.prop")])
    
    # prune tree and dataset based on species intersection 
    tree <- read.tree("../data/formatted.tree.nwk")
    tree$tip.label <- gsub("_", " ", tree$tip.label)
    int <- intersect(tree$tip.label, sub$species)
    if (length(int) <= 2) {
      next
    }
    pruned.tree <- keep.tip(tree, int)
    sub <- sub[sub$species %in% int, ]
    
    # find all possible GLMs
    formula <- reformulate(terms[-1], "rsq")
    glm.models <- dredge(glm(formula, data = sub, na.action = na.fail), subset = dc(x1, x2, x1:x2))
    glm.models <- sortModels(glm.models)
    
    # for each GLM
    for (i in 1:nrow(glm.models)) {
      model <- get.models(glm.models, i)[[1]]
      p <- as.data.frame(t(setNames(data.frame(summary(model)$coefficients)[terms, ]$Pr...t.., 
                    c("intercept.p", "age.p", "rep.p", "interact.p"))))
      res <- setNames(resid(model), sub$species)
      res.physig.p <- as.data.frame(t(setNames(phylosig(pruned.tree, res, method="lambda", test=TRUE)[[4]], "res.physig.p")))
      out <- cbind(head, as.data.frame(glm.models)[i, ], res.physig.p, p)
      lis <- c(lis, list(out))
    }
    
    # find all PGLS models
    cd <- comparative.data(pruned.tree, sub, names.col = "species", vcv = T)
    pgls.models <- dredge(pgls(formula, data = cd), subset = dc(x1, x2, x1:x2))
    pgls.models <- sortModels(pgls.models)
    
    # for each PGLS 
    for (i in 1:nrow(pgls.models)) {
      model <- get.models(pgls.models, i)[[1]]
      p <- as.data.frame(t(setNames(data.frame(summary(model)$coefficients)[terms, ]$Pr...t.., 
                                    c("intercept.p", "age.p", "rep.p", "interact.p"))))
      res.physig.p <- as.data.frame(t(setNames(NA, "res.physig.p")))
      out <- cbind(head, as.data.frame(pgls.models)[i, ], res.physig.p, p)
      lis <- c(lis, list(out))
    }
  }
}

    


# convert list into dataframe
df <- as.data.frame(do.call(rbind, lis))

# write
write.csv(df, file = "../results/AICs.csv", row.names = F)



df <- read.csv("../results/AICs.csv")



