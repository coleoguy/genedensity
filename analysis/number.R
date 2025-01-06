
# slope sig rsq
aaa <- bbb <- ccc <- c()

for (qw in (5:95)*0.01) {
  
  # Zhaobo Hu
  # zhaobohu2002@gmail.com
  
  # Description: Parses results and calculates additional statistics
  # to summarize contigs for each species
  
  dat <- read.csv("../data/data.csv")
  
  # combine raw contig results
  library(data.table)
  dir <- "../results/individual_species_results"
  files <- paste0(dir, "/",  list.files(dir))
  raw <- lapply(files, fread)
  raw <- as.data.frame(rbindlist((raw), fill = TRUE))
  
  # parse by contig size
  raw <- raw[raw$size.Mbp >= 10, ]
  
  # test new method
  parsed <- data.frame()
  for (z in unique(raw$species)) {
    sub <- raw[raw$species == z, ]
    cont <- sum(sub$size.Mbp)
    total <- raw[raw$species == z, ]$asmblysize[1]
    if (cont/total >= qw) {
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
      sdgd <- sd(sub$genedens)
      meangd <- mean(sub$genedens)
      cor <- cor(sub$size.Mbp, sub$genecount)
      weightmean <- sum(sub$genedens * sub$size.Mbp) / sum(sub$size.Mbp)
      weightsd <- sqrt(sum(sub$size.Mbp * (sub$genedens - weightmean)^2) / sum(sub$size.Mbp))
      weightcv <- weightsd / weightmean
      chromsd <- sd(sub$size.Mbp)/mean(sub$size.Mbp)
      contig.stats <- data.frame(species, beta, meangd, sdgd, pval.beta, rsq, cor, weightmean, weightsd, weightcv, chromsd)
    } else {
      beta <- meangd <- sdgd <- pval.beta <- rsq <- cor <- weightmean <- weightsd <- weightcv <- chromsd <- NA
      contig.stats <- data.frame(species, beta, meangd, sdgd, pval.beta, rsq, cor, weightmean, weightsd, weightcv, chromsd)
    }
    final <- rbind(final, merge(merge(dat[dat$species == species, ], contig.stats, by = "species"), sub, by = "species", all = TRUE))
  }
  
  #assign clades
  final$clade <- final$class
  final[final$clade %in% "Aves", ]$clade <- "Sauria"
  final[final$clade %in% "Reptilia", ]$clade <- "Sauria"
  final[!(final$clade %in% c("Actinopterygii", "Mammalia", "Sauria")), ]$clade <- "Others"
  
  # reorder columns
  final <- final[, c(1, 27, 2:11, 26, 12:25)]
  
  # write csv
  write.csv(final, "../results/parsed.csv", row.names = FALSE)
  
  # Zhaobo Hu
  # zhaobohu2002@gmail.com
  
  # Description: Calculates statistics to summarize repeat landscape
  # characteristics for each species. saves one file with parsed
  # contigs and another file for parsed contigs
  
  # calculate stats
  # library(e1071)
  files <- list.files("../results/divsums")
  sp <- gsub("_", " ", gsub(".divsum$", "", files))
  
  dat <- read.csv("../results/parsed.csv")
  asmbsz <- dat[!duplicated(dat$species), ]
  asmbsz <- asmbsz[asmbsz$species %in% sp, ]
  asmbsz <- asmbsz[order(asmbsz$species == sp), ]
  asmbsz <- asmbsz$asmblysize.Mbp*1000000
  
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
    # vector of divergence scores
    div <- table$Div
    # vector of the repeat content of each divergence score
    rep.bp <- rowSums(table[, !names(table) == "Div"])
    # repeat content in Mbp
    totalrep.Mbp <- sum(rep.bp) / 1000000
    # repeat content in percent coverage
    rep.pct <- (rep.bp / asmbsz[i]) * 100
    totalrep.pct <- sum(rep.pct)
    # divergence bin with median repeat
    median <- which(cumsum(rep.pct) > sum(rep.pct)/2)[1]
    
    # unused
    # mean <- sum(divergence*rep.pct)/sum(rep.pct)
    # k <- kurtosis(rep.pct)
    # s <- skewness(rep.pct)
    # max <- max(rep.pct)
    # which <- which.max(rep.pct)
    # s <- skewness(rep.pct)
    # k <- kurtosis(rep.pct)
    
    # build dataframe
    df <- data.frame(species, 
                     totalrep.Mbp,
                     totalrep.pct, 
                     median
    )
    repstats <- rbind(repstats, df)
  }
  df <- merge(dat, repstats, by = "species", all.x = TRUE)
  
  # reorganize and save results
  df <- df[, c(1:23, 28:30, 24:27)]
  write.csv(df,
            "../results/parsed.csv", 
            row.names = FALSE)
  
  
  # before pgls: new parsing method kept 21 mammals and rsq~transformed median has slope -1.6988 and p value 0.016207; no need to exponentiate weights
  # after pgls: 20 species, beta = 1.747489, p = 0.0183, r2 = 0.2731145, predicted rsq diff between highest and lowest medians: 0.37446190
  # what about other clades?
  
  
  # new model
  library(viridis)
  dat <- read.csv("../results/parsed.csv")
  d <- dat
  dat <- dat[!duplicated(dat$species), ]
  dat <- dat[!is.na(dat$chromnum.1n), ]
  dat$median.trans <- 1 - (dat$median/70)
  dat$totalrep.prop <- dat$totalrep.pct * 0.01
  dat <- na.omit(dat[, c("species", "rsq", "class", "order", "family", "clade", "median.trans", "totalrep.prop", "chromnum.1n", "est.gnsz.Mbp", "asmblysize.Mbp")])
  dat <- dat[dat$clade == "Mammalia", ]
  dat <- dat[dat$species != "Callithrix jacchus", ]
  library(phytools)
  tree <- read.tree("../data/formatted_tree.nwk")
  tree$tip.label <- gsub("_", " ", tree$tip.label)
  int <- intersect(tree$tip.label, dat$species)
  pruned.tree <- keep.tip(tree, int)
  dat1 <- dat[dat$species %in% int, ]
  dat1 <- dat1[match(pruned.tree$tip.label, dat1$species), ]
  model <- glm(rsq ~ median.trans, data = dat1)
  res <- setNames(resid(model), dat1$species)
  phylosig(pruned.tree, res, method="lambda", test=TRUE)
  library(nlme)
  aaa <- c(aaa, summary(gls(rsq ~ median.trans, 
              data = dat1))$tTable[2, 1])
  bbb <- c(bbb, summary(gls(rsq ~ median.trans, 
              data = dat1))$tTable[2, 4])
  library(piecewiseSEM)
  ccc <- c(ccc, as.numeric(rsquared(gls(rsq ~ median.trans, data = dat1))[5]))
  

}
df <- data.frame(c((5:95)*0.01), aaa, bbb, ccc)
names(df) <- c("number", "beta", "p", "rsq")
write.csv(df, file = "../results/number.csv", row.names = F)
df <- read.csv("../results/number.csv")


