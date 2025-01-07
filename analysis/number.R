
# slope sig rsq

qwerty <- data.frame()
for (qw in (0:97)*0.01) {
  
  
  thrs <- qw
  
  dat <- read.csv("../data/data.csv")
  
  # combine raw contig results
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
  
  # write csv
  write.csv(final, "../results/parsed.csv", row.names = FALSE)
  # Zhaobo Hu
  # zhaobohu2002@gmail.com
  
  # Description: Calculates statistics to summarize repeat landscape
  # characteristics for each species. saves one file with parsed
  # contigs and another file for parsed contigs
  
  # calculate stats
  files <- list.files("../results/divsums")
  sp <- gsub("_", " ", gsub(".divsum$", "", files))
  
  dat <- read.csv("../results/parsed.csv")
  asmbsz <- dat[!duplicated(dat$species), ]
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
    total.rep.pct <- sum((rep.bp / asmbsz[i]) * 100)
    total.rep.median <- which(cumsum(rep.bp) > sum(rep.bp)/2)[1]
    
    for (k in classes) {
      assign(paste0(tolower(k), ".rep.pct"), sum(table[k] / asmbsz[i] * 100))
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
  df <- merge(dat, repstats, by = "species", all.x = TRUE)
  
  # reorganize and save results
  df <- df[, c(1:20, 25:36, 21:24)]
  write.csv(df,
            "../results/parsed.csv", 
            row.names = FALSE)
  
  
  
  
  qwert <- list()
  classes <- c("total", "line", "sine", "ltr", "dna", "rc")
  for (qwe in classes) {
    dat <- read.csv("../results/parsed.csv")
    d <- dat
    dat <- dat[!duplicated(dat$species), ]
    dat <- dat[!is.na(dat$chromnum.1n), ]
    dat$median.trans <- 1 - (dat[[paste0(qwe, ".rep.median")]]/70)
    dat <- na.omit(dat[, c("species", "rsq", "clade", "median.trans")])
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
    signal <- phylosig(pruned.tree, res, method="lambda", test=TRUE)[4]
    if (signal < 0.05) {
      library(nlme)
      qwert[[paste0(qwe, ".slope")]] <- summary(gls(rsq ~ median.trans, data = dat1))$tTable[2, 1]
      qwert[[paste0(qwe, ".p")]] <- summary(gls(rsq ~ median.trans, data = dat1))$tTable[2, 4]
      library(piecewiseSEM)
      qwert[[paste0(qwe, ".r2")]] <- rsquared(gls(rsq ~ median.trans, data = dat1))[[5]]
    } else {
      qwert[[paste0(qwe, ".slope")]] <- summary(glm(rsq ~ median.trans, data = dat))$coefficients[2, 1]
      qwert[[paste0(qwe, ".p")]] <- summary(glm(rsq ~ median.trans, data = dat))$coefficients[2, 4]
      library(piecewiseSEM)
      qwert[[paste0(qwe, ".r2")]] <- rsquared(glm(rsq ~ median.trans, data = dat))[[5]]
    }
  }
  qwert <- as.data.frame(qwert)
  qwert$number <- qw
  qwerty <- rbind(qwerty, as.data.frame(qwert))
}


df <- data.frame(c((0:97)*0.01), aaa, bbb, ccc)
names(df) <- c("number", "beta", "p", "rsq")
write.csv(df, file = "../results/number.csv", row.names = F)
df <- read.csv("../results/number.csv")


