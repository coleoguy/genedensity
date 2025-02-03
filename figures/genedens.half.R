
# hypothesis: gene density varies WITHIN chromosomes. fission might break a 
# chromosome so that one fragment has high gene density variation and another
# has low gene density variation. 

library(viridis)
# visualize gene density variation on either half of chromosomes
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
  df <- data.frame()
  for (j in dat[dat$species == i, ]$name) { # for each chromosome
    chrom <- annot[annot[, 1] == j, ]
    end <- dat[dat$species == i & dat$name == j, ]$size.Mbp
    left.dens <- nrow(chrom[chrom[, 5] < end/2, ]) / end/2
    right.dens <- nrow(chrom[chrom[, 4] > end/2, ]) / end/2
    df <- rbind(df, data.frame(j, left.dens, right.dens, end))
  }
  mindens <- min(df[, c(2, 3)])
  maxdens <- max(df[, c(2, 3)])
  
  res <- 500
  palette <- viridis(res)
  left.norm <- (df$left.dens - mindens) / (maxdens - mindens)
  left.idx <- round(left.norm * (res-1)) + 1
  left.cols <- palette[left.idx]
  right.norm <- (df$right.dens - mindens) / (maxdens - mindens)
  right.idx <- round(right.norm * (res-1)) + 1
  right.cols <- palette[right.idx]
  
  svg(filename = paste0("half.dens/", gsub(" ", "_", i), ".svg"), 
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
  
  for (k in 1:nrow(df)) {
    rect(ybottom = k - 0.4, 
         xleft = 0, 
         ytop = k + 0.4, 
         xright = df$end[k]/2, 
         col = left.cols[k], 
         border = NA)
    rect(ybottom = k - 0.4, 
         xleft = df$end[k]/2, 
         ytop = k + 0.4, 
         xright = df$end[k], 
         col = right.cols[k], 
         border = NA)
    text(-0.05*df$end[1], k, round(df$left.dens[k], 3), cex = 0.6, xpd = NA, adj = c(0.5, 0.5))
    text(df$end[k] + 0.05*df$end[1], k, round(df$right.dens[k], 3), cex = 0.6, xpd = NA, adj = c(0.5, 0.5))
    #if (k %% 2 == 0) {
    #  text(-0.01*max.chromsize, k, df$j[k], cex = 0.6, xpd = NA, adj = c(1, 0.5))
    #} else {
    #  text(-0.03*max.chromsize, k, df$j[k], cex = 0.6, xpd = NA, adj = c(1, 0.5))
    #}
  }
  
  #rasterImage(as.raster(viridis(res)), 
  #            ybottom = 0.9*k, 
  #            xleft = 0.9*max.chromsize, 
  #            ytop = 0.4*k, 
  #            xright = 1*max.chromsize)
  #text(0.95*max.chromsize, 
  #     k, 
  #     paste0(round(maxdens, 3), " genes/Mb"), 
  #     cex = 0.8, 
  #     xpd = NA, 
  #     adj = c(0.5, 1))
  #text(0.95*max.chromsize, 
  #     0.3*k, 
  #     paste0(round(mindens, 3), " genes.Mb"), 
  #     cex = 0.8, 
  #     xpd = NA, 
  #     adj = c(0.5, 0))
  dev.off()
}
