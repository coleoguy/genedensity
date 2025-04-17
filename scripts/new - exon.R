

# get a csv with number of exons per gene for every gene in every species
source("functions.R")
all.species <- unique(gsub("\\..*$", "", list.files("../data/genomes")))
all.species <- all.species[all.species != "readme"]
df <- data.frame()
for (i in 1:length(all.species)) {
  sp <- all.species[i]
  gtf <- paste0("../data/genomes/", sp, ".gtf")
  gtf <- read.table(gtf, header = FALSE, sep = "\t")
  
  exon <- gtf[gtf[, 3] == "exon", ]
  exon <- exon[order(exon$V5), ]
  exon <- exon[order(exon$V4), ]
  exon <- unique(exon)
  
  gene <- gtf[gtf[, 3] == "gene", ]
  gene <- gene[order(gene$V5), ]
  gene <- gene[order(gene$V4), ]
  gene <- unique(gene)
  
  exons <- data.frame(sp = sp, 
                     exons = apply(X = gene, MARGIN = 1, FUN = howmany, exons = exon))
  df <- rbind(df, exons)
  
  rm(gtf)
  rm(exons)
  gc()
}
write.csv(df, "exon.num.csv", row.names = F)









# for each species, find the exon number at the 80% quantile

# read in repeats.csv and removes species where the exon number at the 80% 
# quantile is greater than 1000; for reference, this number for the human t2t 
# annotation is 404 and the drosophila reference genome annotation is 50
# a <- read.csv("../results/rsq.csv")
# b <- read.csv("../results/repeats.csv")
# int <- intersect(a$species, b$species)
df <- read.csv("../results/exon.num.csv")
all.sp <- unique(df$sp)
result <- c()
for (i in 1:length(all.sp)) {
  sp <- all.sp[i]
  spdat <- df[df$sp == sp, ]
  metric <- quantile(spdat$exons, probs = 0.8, na.rm = TRUE)
  result <- c(result, setNames(metric, sp))
}
result <- result[result <= 1000]
# int <- intersect(int, names(result))
# write.csv(b[b$species %in% int, ], "../results/repeats.csv")


# overall gene density (total # of genes / assembly size) vs assembly size
# larger genomes are less gene dense; exceptions with tetraploids and 
# well-annotated species
rsq <- read.csv("../results/rsq.csv")
num.gene <- c()
for (j in unique(df$sp)) {
  num.gene <- c(num.gene, setNames(nrow(df[df$sp == j, ]), j))
}
int <- intersect(rsq$species, names(num.gene))
rsq <- rsq[rsq$species %in% int, ]
num.gene <- num.gene[names(num.gene) %in% int]
genedens <- num.gene / rsq$assem.sz
plot(rsq$assem.sz, genedens)



