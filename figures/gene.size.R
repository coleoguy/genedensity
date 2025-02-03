

# visualize gaps between genes
dat <- read.csv("../results/parsed.csv")
sp <- gsub("_", " ", sub("\\..*", "", list.files("../data/annot")))
for (i in sp) { # for each species
  annot <- read.table(paste0("../data/annot/", gsub(" ", "_", i), ".gtf"), 
                      header = TRUE, 
                      sep = "\t")
  annot <- annot[annot[, 3] == "gene", ]
  annot[, c(4, 5)] <- annot[, c(4, 5)] / 1000000
  chromnum <- nrow(dat[dat$species == i, ])
  max.chromsize <- max(dat[dat$species == i, ]$size.Mbp)
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
    b <- rep(chromcount - 0.4, nrow(chrom))
    l <- chrom[, 4]
    t <- rep(chromcount + 0.4, nrow(chrom))
    r <- chrom[, 5]
    x <- c(rbind(l, r, r, l, NA))
    y <- c(rbind(b, b, t, t, NA))
    polygon(x, y, col = "black", border = NA)
    
    chromcount <- chromcount + 1
  }
}

sub <- dat[dat$species == "Suricata suricatta", ]
plot(sub$size.Mbp, sub$genecount)
