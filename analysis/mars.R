
terms <- c(
  "(Intercept)", 
  "chromnum.1n", 
  "median.trans", 
  "rep.prop", 
  "chromnum.1n:median.trans", 
  "chromnum.1n:rep.prop", 
  "median.trans:rep.prop", 
  "chromnum.1n:median.trans:rep.prop"
)
lis <- list()

# for each threshold
for (thrs in (0:100)*0.01) {
  
  # read data
  dat <- read.csv("../data/data.csv")
  library(data.table)
  dir <- "../results/individual_species_results"
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
  
  # reorder columns
  final <- final[, c(1, 23, 2:8, 12:17, 22, 24, 9:11, 18:21)]
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
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
  
  # reorganize and save results
  final <- final[, c(1:20, 25:36, 21:24)]
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # LINEs for all clades
  library(caper)
  library(phytools)
  library(MuMIn)

  dat <- final
  dat <- dat[!duplicated(dat$species), ]
  dat <- dat[!is.na(dat$chromnum.1n), ]
  dat <- dat[dat$clade %in% "Mammalia", ]
  mars <- c()
  for (g in 1:nrow(dat)) {
    if (dat[g, ]$order %in% c("Didelphimorphia", 
                              "Paucituberculata", 
                              "Microbiotheria",
                              "Dasyuromorphia", 
                              "Notoryctemorphia",
                              "Peramelemorphia", 
                              "Diprotodontia")){
      mars <- c(mars, TRUE)
      
    } else {
      mars <- c(mars, FALSE)
    }
  }
  dat$mars <- mars
  dat <- na.omit(dat[, c("species", "rsq", "clade", "mars")])
  dat <- dat[dat$species != "Callithrix jacchus", ]
  
  tree <- read.tree("../data/formatted_tree.nwk")
  tree$tip.label <- gsub("_", " ", tree$tip.label)
  int <- intersect(tree$tip.label, dat$species)
  pruned.tree <- keep.tip(tree, int)
  model <- glm(rsq ~ as.factor(mars),data = dat)
  
  rsq <- setNames(dat$rsq, dat$species)
  signal <- phylosig(pruned.tree, rsq, method="lambda", test=TRUE)[[4]]
  
  if (signal < 0.05) {
    sig <- TRUE
    x <- setNames(as.factor(dat$mars), dat$species)
    y <- setNames(dat$rsq, dat$species)
    p <- phylANOVA(pruned.tree, x, y, nsim = 10000, posthoc = TRUE)$Pf
    lis <- c(lis, list(c(thrs, sig, p)))
  } else {
    sig <- FALSE
    p <- summary(aov(rsq ~ mars, data = dat))[[1]][1, 5]
    lis <- c(lis, list(c(thrs, sig, p)))
  }
}



df <- as.data.frame(do.call(rbind, lis))
names(df) <- c("thrs", "p")
write.csv(df, "../results/mars.csv", row.names = F)
l <- read.csv("../results/mars.csv")