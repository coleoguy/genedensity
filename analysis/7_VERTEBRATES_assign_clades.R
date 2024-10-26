

# Zhaobo Hu
# haobohu2002@gmail.com

# Description: pulls the vertebrate clades Actinopteryii, Mammalia, and Sauria
# from class names

# list of classes that will be covered by clade assignments
major.classes <- c("Actinopterygii", "Aves", "Mammalia", "Reptilia")
# read in data file
taxo.gnsz <- read.csv("../data/vertebrates/taxo_gnsz.csv")
# read in all results
final.results <- read.csv("../results/vertebrates/final_results.csv")
# subset of data file that contains non-NA classes
dat <- taxo.gnsz[!is.na(taxo.gnsz$class), ]
# subset containing relevant information
dat <- dat[, c("species", "class")]
# create new column
dat$clade <- NA
# species in class Actinopterygii are assigned to clade Actinopterygii
dat[which(dat$class == major.classes[1]), ]$clade <- "Actinopterygii"
# species in class Aves are assigned to clade Sauria
dat[which(dat$class == major.classes[2]), ]$clade <- "Sauria"
# species in class Mammalia are assigned to clade Mammalia
dat[which(dat$class == major.classes[3]), ]$clade <- "Mammalia"
# species in class Reptilia are assigned to clade Sauria
dat[which(dat$class == major.classes[4]), ]$clade <- "Sauria"
# species not belonging to aforementioned classes are assigned to Others
dat[which(!(dat$class %in% major.classes)), ]$clade <- "Others"
# remove classes from dataframe
dat <- dat[, c("species", "clade")]
# write csv
write.csv(dat, "../results/vertebrates/clades.csv", row.names = FALSE)




