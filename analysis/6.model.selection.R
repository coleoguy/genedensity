


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
      sub$median.trans <- 1 - (sub[[paste0(rep, ".rep.median")]]/70)
      sub$rep.prop <- sub[[paste0(rep, ".rep.pct")]] / 100
      sub <- na.omit(sub[, c("species", "rsq", "clade", "median.trans", "rep.prop")])
      if (cl %in% c("Mammalia", "Actinopterygii", "Sauria")) {
        sub <- sub[sub$clade == cl, ]
      }
      sub <- sub[sub$species != "Callithrix jacchus", ]
      if (nrow(sub) < 2) {
        next
      }
      
      # find intersection between tree tips and results 
      tree <- read.tree("../data/formatted.tree.nwk")
      tree$tip.label <- gsub("_", " ", tree$tip.label)
      int <- intersect(tree$tip.label, sub$species)
      if (length(int) < 2) {
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
write.csv(df, file = "../results/models.csv", row.names = F)


terms <- c(
  "intercept", 
  "median.trans", 
  "rep.prop", 
  "median.trans:rep.prop"
)
dat <- read.csv("../results/models.csv")
names(dat) <- c("clade", "thrs", "repeat", "phylosig", "n", "stat", terms, "f2")


# models where every term is either significant or NA
df <- dat[apply(dat[, terms[-1]], 1, function(x) all(is.na(x) | x < 0.05)), ]

# remove models with all NA terms
df <- df[rowSums(is.na(df[, terms[-1]])) < length(terms[-1]), ]

# filter for significant clade-threshold-repeat combinations
df <- df[df$stat == "p", ]
hit <- c("clade", "thrs", "repeat")
df <- df[, hit]

# get beta coefficients for significant combinations
df <- merge(dat, df, by = intersect(names(dat), names(df)))
df <- df[df$stat == "beta", ]
df <- df[!df$`repeat` %in% "rc", ]
df <- df[order(df$thrs), ]
df <- df[order(df$`repeat`), ]
df <- df[order(df$clade), ]
rownames(df) <- c(1:nrow(df))
df <- df[df$thrs >= 0.8 & df$thrs <= 0.95, ]


#combinations of clades and repeats
comb <- unique(df[, c("clade", "repeat"), ])
beta <- c()
for (i in 1:nrow(comb)) {
  cl <- comb[i, ][[1]]
  rep <- comb[i, ][[2]]
  sub <- df[df$clade == cl, ]
  sub <- sub[sub$`repeat` == rep, ]
  sub <- sub[, c(8, 9)]
  beta <- c(beta, mean(as.matrix(sub)[which(!is.na(sub))]))
}
comb$beta <- beta
comb$abs <- abs(beta)

# define tick marks
n <- 5  
log.max <- log10(max(comb$abs) + 1)  
ticks <- seq(0, log.max, length.out = n)
labels <- 10^ticks - 1  

# color mapping
library(viridis)
dat.trans <- log10(comb$abs + 1)
norm <- dat.trans / log.max  
res <- 1000  
palette <- viridis(res)  
comb$cols <- palette[round(norm * (res - 1)) + 1]  
