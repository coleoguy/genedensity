

# repeat all calculations for 93 contig size cutoff thresholds and 
# find the best model for each repeat type in each clade at each
# cutoff threshold using an exhaustive approach. results stored in
# csv



library(performance)
library(caper)
library(MuMIn)
library(piecewiseSEM)
library(phytools)
options(na.action = "na.fail")
terms <- c(
  "(Intercept)", 
  "median.trans", 
  "rep.prop", 
  "median.trans:rep.prop"
)
lis <- list()
dat1 <- read.csv("../results/parsed.csv")
dat1 <- dat1[!duplicated(dat1$species), ]
dat1 <- dat1[!is.na(dat1$chromnum.1n), ]
thrs <- 0.9

for (cl in c("Total", "Mammalia", "Actinopterygii", "Sauria")) {
  classes <- c("total", "line", "sine", "ltr", "dna", "rc")
  for (qwe in classes) {
    # subset relevant results for analysis
    dat <- dat1
    dat$median.trans <- 1 - (dat[[paste0(qwe, ".rep.median")]]/70)
    dat$rep.prop <- dat[[paste0(qwe, ".rep.pct")]] / 100
    dat <- na.omit(dat[, c("species", "rsq", "clade", "median.trans", "rep.prop")])
    if (cl %in% c("Mammalia", "Actinopterygii", "Sauria")) {
      dat <- dat[dat$clade == cl, ]
    }
    dat <- dat[dat$species != "Callithrix jacchus", ]
    
    # find intersection between tree tips and results and create compdata obj
    tree <- read.tree("../data/formatted_tree.nwk")
    tree$tip.label <- gsub("_", " ", tree$tip.label)
    int <- intersect(tree$tip.label, dat$species)
    pruned.tree <- keep.tip(tree, int)
    
    # initial model selection
    formula <- reformulate(terms[-1], "rsq")
    model <- step(glm(formula, data = dat))
    effects <- data.frame(summary(model)$coefficients)[terms, ]
    # model <- get.models(dredge(glm(formula, data = dat), subset = dc(x1, x2, x1:x2)), 1)[[1]]
    
    # if phylogenetic signals are present in the residuals, use PGLS
    res <- setNames(resid(model), dat$species)
    signal <- phylosig(pruned.tree, res, method="lambda", test=TRUE)[[4]]
    if (signal < 0.05) {
      sig <- TRUE
      cd <- comparative.data(pruned.tree, dat, names.col = "species", vcv = T)
      model <- get.models(dredge(pgls(formula, data = cd), subset = dc(x1, x2, x1:x2)), 1)[[1]]
      effects <- data.frame(summary(model)$coefficients)[terms, ]
      lis <- c(lis, list(c(cl, thrs, qwe, sig, nrow(cd$data), "p", effects$Pr...t.., rsquared(model)[[5]]/(1-rsquared(model)[[5]]))))
      lis <- c(lis, list(c(cl, thrs, qwe, sig, nrow(cd$data), "beta", effects$Estimate, rsquared(model)[[5]]/(1-rsquared(model)[[5]]))))
    } else {
      sig <- FALSE
      lis <- c(lis, list(c(cl, thrs, qwe, sig, nrow(dat), "p", effects$Pr...t.., rsquared(model)[[5]]/(1-rsquared(model)[[5]]))))
      lis <- c(lis, list(c(cl, thrs, qwe, sig, nrow(dat), "beta", effects$Estimate, rsquared(model)[[5]]/(1-rsquared(model)[[5]]))))
    }
  }
}







df <- as.data.frame(do.call(rbind, lis))
names(df) <- c("clade", "thrs", "repeat", "phylosig", "n", "stat", "intercept", terms[-1], "f2")
num <- c("thrs", "intercept", terms[-1], "f2")
df[, num] <- lapply(df[, num], as.numeric)
# some stuff are lists
df[] <- lapply(df, function(x) if(is.list(x)) sapply(x, paste, collapse=",") else x)
write.csv(df, file = "../results/models.chromlevel.csv", row.names = F)


terms <- c(
  "intercept", 
  "median.trans", 
  "rep.prop", 
  "median.trans:rep.prop"
)
df <- read.csv("../results/models.chromlevel.csv")
names(df) <- c("clade", "thrs", "repeat", "phylosig", "n", "stat", terms, "f2")


# models where every term is either significant or NA
df1 <- df[apply(df[, terms[-1]], 1, function(x) all(is.na(x) | x < 0.05)), ]

# remove models with all NA terms
df1 <- df1[rowSums(is.na(df1[, terms[-1]])) < length(terms[-1]), ]

# filter for significant clade-threshold-repeat combinations
df1 <- df1[df1$stat == "p", ]
hit <- c("clade", "thrs", "repeat")
df1 <- df1[, hit]

# get beta coefficients for significant combinations
df2 <- merge(df, df1, by = intersect(names(df), names(df1)))
df2 <- df2[df2$stat == "beta", ]
df2 <- df2[!df2$`repeat` %in% "rc", ]
df2 <- df2[order(df2$thrs), ]
df2 <- df2[order(df2$`repeat`), ]
df2 <- df2[order(df2$clade), ]
rownames(df2) <- c(1:nrow(df2))

