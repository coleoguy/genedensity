
# hypothesis: gene density varies WITHIN chromosomes. fission might break a 
# chromosome so that one fragment has high gene density variation and another
# has low gene density variation. 

# visualize gene start positions
dat <- read.csv("../results/parsed.csv")
sp <- gsub("_", " ", sub("\\..*", "", list.files("../data/annot")))
for (i in sp) { # for each species
  annot <- read.table(paste0("../data/annot/", gsub(" ", "_", i), ".gtf"), 
                      header = TRUE, 
                      sep = "\t")
  annot <- annot[annot[, 3] == "gene", ]
  annot[, c(4, 5)] <- annot[, c(4, 5)] / 1000000
  chromnum <- nrow(dat[dat$species == i, ])
  if (chromnum == 0) {
    next
  }
  max.chromsize <- max(dat[dat$species == i, ]$size.Mbp)
  svg(filename = paste0("genepos/", gsub(" ", "_", i), ".svg"), 
      width = 7, 
      height = 7)
  plot(1, 
       type = "n", 
       xlim = c(0, 1.1*max.chromsize), 
       ylim = c(0, chromnum), 
       xlab = "Size (Mb)", 
       ylab = "Chromosomes", 
       axes = FALSE, 
       main = i)
  axis(1)
  chromcount <- 1
  for (j in dat[dat$species == i, ]$name) { # for each chromosome
    chrom <- annot[annot[, 1] == j, ]
    end <- dat[dat$species == i & dat$name == j, ]$size.Mbp
    # rect(ybottom = chromcount - 0.4, 
         # xleft = 0, 
         # ytop = chromcount + 0.4, 
         # xright = end, 
         # col = "white", 
         # border = "black")
    if (chromcount %% 2 == 0) {
      text(-0.01*max.chromsize, chromcount, j, cex = 0.6, xpd = NA, adj = c(1, 0.5))
    } else {
      text(-0.03*max.chromsize, chromcount, j, cex = 0.6, xpd = NA, adj = c(1, 0.5))
    }
    mid <- chrom[, 4]+((chrom[, 5]-chrom[, 4])/2)
    b <- rep(chromcount - 0.4, nrow(chrom))
    l <- mid - 0.00004*end
    t <- rep(chromcount + 0.4, nrow(chrom))
    r <- mid + 0.00004*end
    x <- c(rbind(l, r, r, l, NA))
    y <- c(rbind(b, b, t, t, NA))
    polygon(x, y, col = "black", border = NA)
    
    chromcount <- chromcount + 1
  }
  dev.off()
}





















# calculate CVs for gene positions on the 0.5n largest chromosomes and calculate 
# the median of those CVs

dat <- read.csv("../results/parsed.csv")
sp <- gsub("_", " ", sub("\\..*", "", list.files("../data/annot")))
df <- data.frame()
for (i in sp) { # for each species
  df1 <- data.frame()
  annot <- read.table(paste0("../data/annot/", gsub(" ", "_", i), ".gtf"), 
                      header = TRUE, 
                      sep = "\t")
  annot <- annot[annot[, 3] == "gene", ]
  annot[, c(4, 5)] <- annot[, c(4, 5)] / 1000000
  chromnum <- nrow(dat[dat$species == i, ])
  if (chromnum == 0) {
    next
  }
  top <- round(0.5 * chromnum, 0)
  for (j in dat[dat$species == i, ]$name[1:top]) { # for each chromosome
    chrom <- annot[annot[, 1] == j, ]
    end <- dat[dat$species == i & dat$name == j, ]$size.Mbp
    mid <- chrom[, 4]+((chrom[, 5]-chrom[, 4])/2)
    cv <- sd(mid) / mean(mid)
    df1 <- rbind(df1, data.frame(i, cv))
  }
  med <- median(df1$cv)
  df <- rbind(df, data.frame(i, med))
}

dat <- dat[!duplicated(dat$species), ]
dat <- merge(dat, df, by.y = "i", by.x = "species")
dat <- dat[, c("species", "med", "rsq")]
plot(dat$med, dat$rsq)
model <- glm(rsq ~ med, data = dat)
summary(model)

library(phytools)
library(caper)
tree <- read.tree("../data/formatted_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
int <- intersect(tree$tip.label, dat$species)
pruned.tree <- keep.tip(tree, int)
res <- setNames(resid(model), dat$species)
phylosig(pruned.tree, 
         res, 
         method = "lambda", 
         test = TRUE, 
         nsim = 10000)
cd <- comparative.data(pruned.tree, 
                       dat, 
                       names.col = "species", 
                       vcv = TRUE)
summary(pgls(rsq ~ med, data = cd))

