

library(data.table)
source("../scripts/functions.R")

verbose <- T

chrom <- read.csv("../data/chrom.and.gsz.csv")

# rsq scatter

sp <- c("Ornithorhynchus_anatinus", "Theropithecus_gelada")
for (i in 1:length(sp)) {
  fasta.path <- paste0("../data/genomes/", sp[i], ".fa")
  annot.path <- paste0("../data/genomes/", sp[i], ".gtf")
  max.contig <- 60
  fasta.data <- dataFromFasta(fasta.path = fasta.path, 
                              max.contig = max.contig,
                              verbose = TRUE)
  fasta.data <- fasta.data[fasta.data$size.Mb >= 10, ]
  name <- fasta.data$name
  size.Mb <- fasta.data$size.Mb
  asmblysize.Mb <- fasta.data$asmblysize.Mb[1]
  rm(fasta.data)
  gc()
  genecount <- dataFromGtf(annot.path, name, verbose) / 1000
  assign(paste0("sp", i), 
         data.frame(sp[i], size.Mb, genecount))
}


par(mfrow = c(1, 2))
par(mar = c(5.1, 4.1, 1, 0.9))
color <- adjustcolor("#405070", alpha.f = 0.7)
plot(sp1$size.Mb, 
     sp1$genecount, 
     pch = 16, 
     cex = 0.8, 
     col = color, 
     ylim = c(0.25, 3), 
     xlim = c(0, 250), 
     xlab = NA, 
     ylab = "Gene count (thousands)")
abline(glm(genecount ~ size.Mb, data = sp1), lwd = 1.5, col = color)
par(mar = c(5.1, 2.9, 1, 2.1))
color <- adjustcolor("#405070", alpha.f = 0.7)
plot(sp2$size.Mb, 
     sp2$genecount, 
     pch = 16, 
     cex = 0.8, 
     col = color, 
     ylim = c(0.25, 3), 
     xlim = c(0, 250), 
     xlab = NA, 
     ylab = NA)
abline(glm(genecount ~ size.Mb, data = sp2), lwd = 1.5, col = color)

#text(-40, -0.01, "Chromosome size (Mb)", xpd = NA, adj = 0)
par(mfrow = c(1, 1))
par(mar = c(5.1, 4.1, 4.1, 2.1))
