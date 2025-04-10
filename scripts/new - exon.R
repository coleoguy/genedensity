

library(data.table)
all.species <- unique(gsub("\\..*$", "", list.files("../data/genomes")))
all.species <- all.species[all.species != "readme"]
df <- data.frame()
for (i in 1:length(all.species)) {
  sp <- all.species[i]
  gtf <- paste0("../data/genomes/", sp, ".gtf")
  gtf <- read.table(gtf, header = FALSE, sep = "\t")
  exon <- gtf[gtf[, 3] == "exon", ]
  exonlength <- exon[, 5] - exon[, 4]
  df <- rbind(
    df,
    data.frame(
      sp,
      genenum = nrow(gtf[gtf[, 3] == "gene", ]),
      exonnum = length(exonlength),
      exontotlen = sum(exonlength),
      exonavglen = mean(exonlength),
      exonmedlen = median(exonlength)
    )
  )
}
write.csv(df, "../results/exons.csv", row.names = FALSE)





library(data.table)
library(parallel)
all.species <- unique(gsub("\\..*$", "", list.files("../data/genomes")))
all.species <- all.species[all.species != "readme"]
do <- function(sp) {
  gtf <- paste0("../data/genomes/", sp, ".gtf")
  gtf <- read.table(gtf, header = FALSE, sep = "\t")
  exon <- gtf[gtf[, 3] == "exon", ]
  exonlength <- exon[, 5] - exon[, 4]
  data.frame(
    sp = sp,
    genenum = nrow(gtf[gtf[, 3] == "gene", ]),
    exonnum = length(exonlength),
    exontotlen = sum(exonlength),
    exonavglen = mean(exonlength),
    exonmedlen = median(exonlength)
  )
}

ncores <- 8
cl <- makeCluster(ncores)
clusterEvalQ(cl, library(data.table))
clusterExport(cl, 
              varlist = "all.species",
              envir = environment())
results <- parLapply(cl, all.species, do)
results <- rbindlist(results)














exon <- read.csv("../results/exons.csv")
rsq <- read.csv("../results/rsq.csv")
merg <- merge(rsq, exon, by.x = "species", by.y = "sp")
merg$assem.sz <- merg$assem.sz * 1000000
merg$vec <- merg$totlength/merg$assem.sz
merg$vec2 <- merg$number/merg$assem.sz

