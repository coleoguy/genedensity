species <- read.csv("../../../data/chrom-and-gsz.csv")
r2 <- read.csv("../../../results/rsq.csv")
landscape <- read.csv("../../../results/repeat-results.csv")
busco <- read.csv("../../../data/busco.csv")

clades <- species$class
clades[clades %in% c("Reptilia", "Aves")] <- "Sauropsida"
clades[!clades %in% c("Mammalia", "Actinopterygii", "Sauropsida")] <- "Other"

species$has.landscape <- NA
species$has.landscape[match(landscape$species, species$species)] <- "Yes"
species$has.landscape[which(is.na(species$has.landscape))] <- "No"
species <- merge(species, r2, by = "species", all = T)
species <- merge(species, busco, by = "species", all = T)
species <- species[, c("species", "clade", "rsq", "assem.sz", "has.landscape", "single")]
species$rsq <- signif(species$rsq, 3)
species$assem.sz <- signif(species$assem.sz, 4)
species$clade <- clades
species$species <- gsub("_", " ", species$species)
colnames(species) <- c("Species", "Clade", "R²", "Assembly size (Mb)", "Repeats analyzed?", "BUSCO (Single)")
write.csv(species, "table-s1.csv", row.names = F)
