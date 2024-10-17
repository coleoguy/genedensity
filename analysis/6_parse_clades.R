

# Zhaobo Hu
# haobohu2002@gmail.com

# Description: removes species without genome size estimates and then 
# pulls the vertebrate clades Actinopteryii, Mammalia, and Sauria from
# class names

major.classes <- c("Actinopterygii", "Aves", "Mammalia", "Reptilia")
clades.gnsz <- read.csv("../data/vertebrates/clades_gnsz.csv")
if ("clade" %in% colnames(clades.gnsz)) {
  clades.gnsz <- subset(clades.gnsz, select = -clade)
}
clades.gnsz <- clades.gnsz[!is.na(clades.gnsz$genome.size.est_bp), ]
clade <- c()
clade[which(clades.gnsz$class == major.classes[1])] <- "Actinopterygii"
clade[which(clades.gnsz$class == major.classes[2])] <- "Sauria"
clade[which(clades.gnsz$class == major.classes[3])] <- "Mammalia"
clade[which(clades.gnsz$class == major.classes[4])] <- "Sauria"
clade[which(!(clades.gnsz$class %in% major.classes))] <- "Others"

clades.gnsz <- cbind(clades.gnsz[1], clade, clades.gnsz[, 2:ncol(clades.gnsz)])
write.csv(clades.gnsz, "../data/vertebrates/clades_gnsz.csv", row.names = FALSE)




