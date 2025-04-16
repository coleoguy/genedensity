

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
a <- read.csv("../results/rsq.csv")
b <- read.csv("../results/repeats.csv")
int <- intersect(a$species, b$species)
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
int <- intersect(int, names(result))
write.csv(b[b$species %in% int, ], "../results/repeats.csv")

library(phytools)
library(caper)
tree <- read.tree("../data/formatted.tree.nwk")
int2 <- intersect(tree$tip.label, int)
pruned.tree <- keep.tip(tree, int2)

y <- a[a$species %in% int2, ]$rsq
x <- result[names(result) %in% int2]
x <- x[order(names(x))]
res <- setNames(resid(glm(y ~ x)), sort(int2))

l <- phylosig(pruned.tree, res, method = "lambda", niter = 10000, test = T)

dat <- data.frame(sp = names(x), 
                  rsq = y, 
                  metric = x)
cd <- comparative.data(pruned.tree, 
                       dat, 
                       names.col = "sp", 
                       vcv = T)
model <- pgls(rsq ~ metric, data = cd)
summary(model)
