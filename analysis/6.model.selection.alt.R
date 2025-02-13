


# for each repeat in each clade in each threshold, fit a glm model where the 
# response is r2 and the predictors are (divergence containing 1/2 repeat sum), 
# repeat proportion, and the interaction of these two terms. then, select
# for the best model using step(). if the residuals of this model show 
# phylogenetic signal, fit a PGLS model instead and select the best model using
# AIC (step() doesn't work for PGLS objects)


library(caper) # apply PGLS
library(MuMIn) # calculate AIC
library(piecewiseSEM) # calculate R2 for PGLS objects
library(phytools) # phy stuff
terms <- c(
  "(Intercept)", 
  "median.trans", 
  "rep.prop", 
  "median.trans:rep.prop"
)
dat <- read.csv("../results/parsed.csv")
dat <- dat[!is.na(dat$chromnum.1n), ]

# use a list to save memory
lis <- list()

# for each threshold
for (thrs in (0:100)*0.01) {
  # for each clade
  for (cl in c("Total", "Mammalia", "Actinopterygii", "Sauria")) {
    classes <- c("total", "line", "sine", "ltr", "dna", "rc")
    # for each repeat
    for (rep in classes) {
      
      # subset data
      sub <- dat[dat$thrs == thrs, ]
      sub <- sub[!duplicated(sub$species), ]
      sub <- na.omit(sub[, c("species", "rsq", "clade", paste0(rep, ".rep.median"), paste0(rep, ".rep.pct"))])
      if (cl %in% c("Mammalia", "Actinopterygii", "Sauria")) {
        sub <- sub[sub$clade == cl, ]
      }
      sub <- sub[sub$species != "Callithrix jacchus", ]
      if (nrow(sub) <= 2) {
        next
      }
      
      # transform data
      median <- 1 - (sub[[paste0(rep, ".rep.median")]]/70)
      median.range <- range(na.omit(median))
      sub$median.trans <- (median - median.range[1]) / diff(median.range)
      vol <- sub[[paste0(rep, ".rep.pct")]] / 100
      vol.range <- range(na.omit(vol))
      sub$rep.prop <- (vol - vol.range[1]) / diff(vol.range)
      sub <- na.omit(sub[, c("species", "rsq", "clade", "median.trans", "rep.prop")])
      
      # find intersection between tree tips and results 
      tree <- read.tree("../data/formatted.tree.nwk")
      tree$tip.label <- gsub("_", " ", tree$tip.label)
      int <- intersect(tree$tip.label, sub$species)
      if (length(int) <= 2) {
        next
      }
      pruned.tree <- keep.tip(tree, int)
      
      # initial model selection
      formula <- reformulate(terms[-1], "rsq")
      model <- step(glm(formula, data = sub))
      effects <- data.frame(summary(model)$coefficients)[terms, ]
      # model <- get.models(dredge(glm(formula, data = dat), subset = dc(x1, x2, x1:x2)), 1)[[1]]
      
      # if phylogenetic signals are present in model residuals, apply PGLS
      res <- setNames(resid(model), sub$species)
      if (all(res == 0)) {
        next
      }
      signal <- phylosig(pruned.tree, res, method="lambda", test=TRUE)[[4]]
      if (signal < 0.05) {
        sig <- TRUE
        cd <- comparative.data(pruned.tree, sub, names.col = "species", vcv = T)
        model <- get.models(dredge(pgls(formula, data = cd), 
                                   subset = dc(x1, x2, x1:x2)), 1)[[1]]
        effects <- data.frame(summary(model)$coefficients)[terms, ]
        lis <- c(lis, list(c(cl, thrs, rep, sig, nrow(cd$data), "p", effects$Pr...t.., rsquared(model)[[5]]/(1-rsquared(model)[[5]]))))
        lis <- c(lis, list(c(cl, thrs, rep, sig, nrow(cd$data), "beta", effects$Estimate, rsquared(model)[[5]]/(1-rsquared(model)[[5]]))))
        
      # if model residuals don't show signal, use initial model
      } else {
        sig <- FALSE
        lis <- c(lis, list(c(cl, thrs, rep, sig, nrow(sub), "p", effects$Pr...t.., rsquared(model)[[5]]/(1-rsquared(model)[[5]]))))
        lis <- c(lis, list(c(cl, thrs, rep, sig, nrow(sub), "beta", effects$Estimate, rsquared(model)[[5]]/(1-rsquared(model)[[5]]))))
      }
    }
  }
}

# convert list into dataframe
df <- as.data.frame(do.call(rbind, lis))
names(df) <- c("clade", "thrs", "repeat", "phylosig", "n", "stat", "intercept", terms[-1], "f2")
num <- c("thrs", "intercept", terms[-1], "f2")
df[, num] <- lapply(df[, num], as.numeric)
df[] <- lapply(df, function(x) if(is.list(x)) sapply(x, paste, collapse=",") else x)

# write
write.csv(df, file = "../results/models1.csv", row.names = F)







