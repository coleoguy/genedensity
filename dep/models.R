

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
  #"chromnum.1n", 
  "median.trans", 
  "rep.prop", 
  #"chromnum.1n:median.trans", 
  #"chromnum.1n:rep.prop", 
  "median.trans:rep.prop"
  #"chromnum.1n:median.trans:rep.prop"
)
lis <- list()
dat <- read.csv("../results/parsed.csv")
tree <- read.tree("../data/formatted.tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
# for each threshold
for (thrs in (80)*0.01) {
  
  sub <- dat[dat$thrs == thrs, ]
  
  for (cl in c("Total", "Mammalia", "Actinopterygii", "Sauria")) {
    classes <- c("total", "line", "sine", "ltr", "dna", "rc")
    for (qwe in classes) {
      
      # subset relevant results for analysis
      sub <- sub[!duplicated(sub$species), ]
      sub <- sub[!is.na(dat$chromnum.1n), ]
      sub$median.trans <- 1 - (sub[[paste0(qwe, ".rep.median")]]/70)
      sub$rep.prop <- sub[[paste0(qwe, ".rep.pct")]] / 100
      sub <- na.omit(sub[, c("species", "rsq", "clade", "median.trans", "rep.prop", "chromnum.1n")])
      if (cl %in% c("Mammalia", "Actinopterygii", "Sauria")) {
        dat <- dat[dat$clade == cl, ]
      }
      sub <- sub[sub$species != "Callithrix jacchus", ]
      
      # find intersection between tree tips and results and create compdata obj
      int <- intersect(tree$tip.label, dat$species)
      pruned.tree <- keep.tip(tree, int)
      
      # initial model selection
      formula <- reformulate(terms[-1], "rsq")
      model <- get.models(dredge(glm(formula, data = dat)), 1)[[1]]
      
      # if phylogenetic signals are present in the residuals, use PGLS
      res <- setNames(resid(model), dat$species)
      signal <- phylosig(pruned.tree, res, method="lambda", test=TRUE)[[4]]
      if (signal < 0.05) {
        sig <- TRUE
        cd <- comparative.data(pruned.tree, dat, names.col = "species", vcv = T)
        model <- get.models(dredge(pgls(formula, data = cd)), 1)[[1]]
        effects <- data.frame(summary(model)$coefficients)[terms, ]
        lis <- c(lis, list(c(cl, thrs, qwe, sig, nrow(cd$data), "p", effects$Pr...t.., rsquared(model)[[5]])))
        lis <- c(lis, list(c(cl, thrs, qwe, sig, nrow(cd$data), "beta", effects$Estimate, rsquared(model)[[5]])))
      } else {
        sig <- FALSE
        lis <- c(lis, list(c(cl, thrs, qwe, sig, nrow(dat), "p", effects$Pr...t.., rsquared(model)[[5]])))
        lis <- c(lis, list(c(cl, thrs, qwe, sig, nrow(dat), "beta", effects$Estimate, rsquared(model)[[5]])))
      }
    }
  }
}









df <- as.data.frame(do.call(rbind, lis))
names(df) <- c("clade", "thrs", "repeat", "phylosig", "n", "stat", "intercept", terms[-1], "R2")
num <- c("thrs", "intercept", terms[-1], "R2")
df[, num] <- lapply(df[, num], as.numeric)
# some stuff are lists
df[] <- lapply(df, function(x) if(is.list(x)) sapply(x, paste, collapse=",") else x)
write.csv(df, file = "../results/models.csv", row.names = F)


terms <- c(
  "intercept", 
  "chromnum.1n", 
  "median.trans", 
  "rep.prop", 
  "chromnum.1n:median.trans", 
  "chromnum.1n:rep.prop", 
  "median.trans:rep.prop", 
  "chromnum.1n:median.trans:rep.prop"
)
df <- read.csv("../results/models.csv")
names(df) <- c("clade", "thrs", "repeat", "phylosig", "n", "stat", terms, "R2")


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



df2 <- df2[!df2$clade %in% "Sauria", ]
df2 <- df2[!df2$`repeat` %in% "rc", ]

#combinations of clades and repeats
df3 <- unique(df2[, c("clade", "repeat"), ])
l <- list()
for (i in 1:nrow(df3)) {
  cl <- df3[i, ][[1]]
  rep <- df3[i, ][[2]]
  sub <- df2[df2$clade == cl, ]
  sub <- sub[sub$`repeat` == rep, ]
  sub <- sub[sub$thrs == round(mean(sub$thrs), 2), ][, 8:14]
  beta <- sub[which(!is.na(sub))]
  l <- c(l, list(c(cl, rep, beta)))
}
df4 <- as.data.frame(do.call(rbind, l))
vec <- abs(unlist(df4$median.trans))
library(viridis)
cols <- viridis(300)
image(1, seq(0, 3, length.out = 300), t(seq(0, 3, length.out = 300)), col = cols, axes = FALSE, xlab = "", ylab = "")
axis(4)



num <- vec / 3
cols <- viridis(300)
ramp <- colorRamp(cols, interpolate = "linear")
rgb(ramp(num), maxColorValue = 255)


"#1FA286" "#E1E318" "#23898D" "#1F998A" "#23888E" "#23888E"






library(phytools)
library(caper)
dat <- read.csv("../results/parsed.csv")
dat <- dat[dat$thrs == 0.85, ]
dat <- dat[!duplicated(dat$species), ]
dat <- dat[dat$clade == "Mammalia", ]
dat <- na.omit(dat[, c("species", "clade", "rsq", "total.rep.median")])
tree <- read.tree("../data/formatted_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
int <- intersect(dat$species, tree$tip.label)
pruned.tree <- keep.tip(tree, int)
model <- glm(dat$rsq ~ dat$total.rep.median)
res <- resid(model)
sig <- phylosig(pruned.tree, 
                setNames(res, dat$species), 
                method = "lambda", 
                nsim = 10000, 
                test = TRUE)[[4]]
if (sig < 0.05) {
  cd <- comparative.data(pruned.tree, dat, names.col = "species", vcv = TRUE)
  model <- pgls(rsq ~ total.rep.median, data = cd)
}
plot(x = dat$total.rep.median, 
     y = dat$rsq,
     xlab = "Repeat Age",
     ylab = "Gene Density Variation",
     pch = 16)
abline(model)

