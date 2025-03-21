

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

# for each threshold
for (thrs in (0:100)*0.01) {
  
  # read data
  dat <- read.csv("../data/data.csv")
  library(data.table)
  dir <- "../results/indiv.contigs"
  files <- paste0(dir, "/",  list.files(dir))
  contigs <- lapply(files, fread)
  contigs <- as.data.frame(rbindlist((contigs), fill = TRUE))
  
  # parse by contig size
  contigs <- contigs[contigs$size.Mbp >= 10, ]
  
  # remove species with less than 2 contigs
  rm <- names(table(contigs$species)[table(contigs$species) < 2])
  contigs <- contigs[!(contigs$species %in% rm), ]
  
  # test new method
  parsed <- data.frame()
  for (z in unique(contigs$species)) {
    sub <- contigs[contigs$species == z, ]
    cont <- sum(sub$size.Mbp)
    total <- contigs[contigs$species == z, ]$asmblysize[1]
    if (cont/total >= thrs) {
      parsed <- rbind(parsed, sub)
    }
  }
  
  # calculate stats based on parsed results
  sp <- unique(parsed$species)
  final <- data.frame()
  for (species in sp) {
    i <- species
    sub <- parsed[which(parsed$species == i), ]
    if (nrow(sub) > 0){
      fit <- summary(glm(sub$genecount ~ sub$size.Mbp))
      beta <- fit$coefficients[2, 1]
      pval.beta <- fit$coefficients[2, 4]
      rsq <- summary(lm(sub$genecount ~ sub$size.Mbp))$r.squared
      weightmean <- sum(sub$genedens * sub$size.Mbp) / sum(sub$size.Mbp)
      weightsd <- sqrt(sum(sub$size.Mbp * (sub$genedens - weightmean)^2) / sum(sub$size.Mbp))
      weightcv <- weightsd / weightmean
      contig.stats <- data.frame(species, beta, pval.beta, rsq, weightmean, weightsd, weightcv)
    } else {
      beta <- pval.beta <- rsq <- weightmean <- weightsd <- weightcv <- NA
      contig.stats <- data.frame(species, beta, pval.beta, rsq, weightmean, weightsd, weightcv)
    }
    final <- rbind(final, merge(merge(dat[dat$species == species, ], contig.stats, by = "species"), sub, by = "species", all = TRUE))
  }
  
  #assign clades
  final$clade <- final$class
  final[final$clade %in% "Aves", ]$clade <- "Sauria"
  final[final$clade %in% "Reptilia", ]$clade <- "Sauria"
  final[!(final$clade %in% c("Actinopterygii", "Mammalia", "Sauria")), ]$clade <- "Others"
  
  final$cont.asmb.rat.cutoff <- thrs
  
  files <- list.files("../results/divsums")
  sp <- gsub("_", " ", gsub(".divsum$", "", files))
  asmbsz <- final[!duplicated(final$species), ]
  asmbsz <- asmbsz[asmbsz$species %in% sp, ]
  asmbsz <- setNames(asmbsz$asmblysize.Mbp*1000000, asmbsz$species)
  repstats <- data.frame()
  for (i in 1:length(sp)) {
    species <- sp[i]
    # read text file into lines
    lines <- readLines(paste0("../results/divsums/", files[i]))
    # look for the start of relevant information
    phrase <- "Coverage for each repeat class and divergence (Kimura)"
    start.index <- match(phrase, lines) + 1
    # condense relevant lines into a table
    lines <- lines[start.index:length(lines)]
    table <- read.table(textConnection(lines), 
                        sep = " ", 
                        header = TRUE)
    # drop NA columns
    table <- table[-c(which(sapply(table, function(col) all(is.na(col)))))]
    
    # condense table
    classes <- c("LINE", "SINE", "LTR", "DNA", "RC", "Div", "Unknown")
    for (j in classes) {
      pat <- paste0("^", j, "(\\.|$)")
      headers <- grep(pat, names(table), value = TRUE)
      sub <- table[, headers]
      sums <- rowSums(as.matrix(sub))
      table <- table[, !names(table) %in% headers]
      assign(j, sums)
    }
    Others <- rowSums(as.matrix(table))
    table <- data.frame(Div, LINE, SINE, LTR, DNA, RC, Others, Unknown)
    
    # all repeat total and median
    rep.bp <- rowSums(table[, !names(table) == "Div"])
    total.rep.pct <- sum((rep.bp / asmbsz[sp[i]]) * 100)
    
    total.rep.median <- which(cumsum(rep.bp) > sum(rep.bp)/2)[1]
    
    for (k in classes) {
      assign(paste0(tolower(k), ".rep.pct"), sum(table[k] /  asmbsz[sp[i]] * 100))
      assign(paste0(tolower(k), ".rep.median"), which(cumsum(table[k]) > sum(table[k])/2)[1])
    }
    
    # build dataframe
    df <- data.frame(species, 
                     total.rep.pct, 
                     total.rep.median, 
                     line.rep.pct, 
                     line.rep.median, 
                     sine.rep.pct, 
                     sine.rep.median, 
                     ltr.rep.pct, 
                     ltr.rep.median, 
                     dna.rep.pct, 
                     dna.rep.median, 
                     rc.rep.pct, 
                     rc.rep.median
    )
    repstats <- rbind(repstats, df)
  }
  final <- merge(final, repstats, by = "species", all.x = TRUE)

  for (cl in c("Total", "Mammalia", "Actinopterygii", "Sauria")) {
    classes <- c("total", "line", "sine", "ltr", "dna", "rc")
    for (qwe in classes) {
      
      # subset relevant results for analysis
      dat <- final[!duplicated(final$species), ]
      dat <- dat[!is.na(dat$chromnum.1n), ]
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
}









df <- as.data.frame(do.call(rbind, lis))
names(df) <- c("clade", "thrs", "repeat", "phylosig", "n", "stat", "intercept", terms[-1], "f2")
num <- c("thrs", "intercept", terms[-1], "f2")
df[, num] <- lapply(df[, num], as.numeric)
# some stuff are lists
df[] <- lapply(df, function(x) if(is.list(x)) sapply(x, paste, collapse=",") else x)
write.csv(df, file = "../results/models.no.chromnum.csv", row.names = F)


terms <- c(
  "intercept", 
  "median.trans", 
  "rep.prop", 
  "median.trans:rep.prop"
)
df <- read.csv("../results/models.no.chromnum.csv")
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
df2 <- df2[df2$thrs >= 0.7 & df2$thrs <= 0.9, ]


#combinations of clades and repeats
df3 <- unique(df2[, c("clade", "repeat"), ])
l <- list()
for (i in 1:nrow(df3)) {
  cl <- df3[i, ][[1]]
  rep <- df3[i, ][[2]]
  sub <- df2[df2$clade == cl, ]
  sub <- sub[sub$`repeat` == rep, ]
  sub <- sub[, c(8, 9)]
  beta <- mean(as.matrix(sub)[which(!is.na(sub))])
  l <- c(l, list(c(cl, rep, beta)))
}
df4 <- as.data.frame(do.call(rbind, l))
vec <- as.numeric(df4$V3)
library(viridis)
cols <- viridis(300)
image(1, seq(0, 3, length.out = 300), t(seq(0, 3, length.out = 300)), col = cols, axes = FALSE, xlab = "", ylab = "")
axis(4)



num <- (abs(vec)-min(abs(vec))) / diff(range(abs(vec)))
cols <- viridis(300)
ramp <- colorRamp(cols, interpolate = "linear")
rgb(ramp(num), maxColorValue = 255)


[1] "#440155" "#FDE725" "#2A768E" "#471264" "#450558" "#450357" "#481A6C"
[8] "#440154" "#440256"



# Install and load necessary packages
if (!requireNamespace("viridis", quietly = TRUE)) install.packages("viridis")
if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")
library(viridis)
library(scales)

# Your data vector
vec <- c(1.492454, 21.466283, 9.323701, 2.367490, 1.667414, 1.576235, 2.788643, 1.431416, 1.542669)

# Generate colors using a log-transformed scale
log_colors <- viridis(length(vec))[rank(log(vec))]

# Example of plotting points with log-transformed color scale
plot(vec, pch = 16, col = log_colors, cex = 2, main = "Log-Transformed Color Scale", xlab = "Index", ylab = "Value")

# Optional: Add a legend
legend("topright", legend = signif(range(vec), 2), fill = viridis(2), title = "Log Scale")











# Your data vector
vec <- c(1.492454, 21.466283, 9.323701, 2.367490, 1.667414, 1.576235, 2.788643, 1.431416, 1.542669)

# Define the range of values for the color scale
vec_range <- range(vec)

# Generate a sequence of values for the scale (log-transformed)
scale_values <- seq(log(vec_range[1]), log(vec_range[2]), length.out = 100)

# Generate corresponding colors using viridis
scale_colors <- viridis(length(scale_values))

# Create tick positions on the log scale
tick_positions <- seq(log(vec_range[1]), log(vec_range[2]), length.out = 5)

# Convert log tick positions back to original scale for labeling
tick_labels <- signif(exp(tick_positions), 2)

# Plot the color bar
par(mar = c(5, 1, 1, 5)) # Adjust margins for color bar
image(
  1, scale_values, t(as.matrix(scale_values)), 
  col = scale_colors, 
  axes = FALSE, 
  xlab = "", 
  ylab = "Log Scale", 
  main = "Color Scale"
)

# Add axis labels at the specified tick positions
axis(4, at = tick_positions, labels = tick_labels)
box()


"#472D7BFF" "#FDE725FF" "#AADC32FF" "#27AD81FF" "#21908CFF" "#2C728EFF"
[7] "#5DC863FF" "#440154FF" "#3B528BFF"






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

