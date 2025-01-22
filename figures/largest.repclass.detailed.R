





# what is the largest repeat class in each species? not like LINEs 
# SINEs etc but more specific like TC/mariner hAT L2 Gypsy etc



name <- cl <- r2 <- largest <- c()
dat <- read.csv("../data/data.csv")

# combine raw contig results
library(data.table)
dir <- "../results/individual_species_results"
files <- paste0(dir, "/",  list.files(dir))
contigs <- lapply(files, fread)
contigs <- as.data.frame(rbindlist((contigs), fill = TRUE))

# parse by contig size
contigs <- contigs[contigs$size.Mbp >= 10, ]

# remove species with less than 3 contigs
rm <- names(table(contigs$species)[table(contigs$species) < 3])
contigs <- contigs[!(contigs$species %in% rm), ]

df <- data.frame()
parsed <- data.frame()
for (z in unique(contigs$species)) {
  sub <- contigs[contigs$species == z, ]
  cont <- sum(sub$size.Mbp)
  total <- contigs[contigs$species == z, ]$asmblysize[1]
  if (cont/total >= 0) {
    parsed <- rbind(parsed, sub)
  }
}

# calculate stats based on parsed results
sp <- unique(parsed$species)
for (species in sp) {
  sub <- parsed[which(parsed$species == species), ]
  # fit <- summary(glm(sub$genecount ~ sub$size.Mbp))
  # beta <- fit$coefficients[2, 1]
  # pval.beta <- fit$coefficients[2, 4]
  rsq <- summary(lm(sub$genecount ~ sub$size.Mbp))$r.squared
  # weightmean <- sum(sub$genedens * sub$size.Mbp) / sum(sub$size.Mbp)
  # weightsd <- sqrt(sum(sub$size.Mbp * (sub$genedens - weightmean)^2) / sum(sub$size.Mbp))
  # weightcv <- weightsd / weightmean
  stats <- data.frame(species, rsq, 0)
  sub <- merge(sub, stats, by = "species", all = TRUE)
  sub <- merge(dat[dat$species == species, ], sub, by = "species", all = TRUE)
  df <- rbind(df, sub)
}

#assign clades
df$clade <- df$class
df[df$clade %in% "Aves", ]$clade <- "Sauria"
df[df$clade %in% "Reptilia", ]$clade <- "Sauria"
df[!(df$clade %in% c("Actinopterygii", "Mammalia", "Sauria")), ]$clade <- "Others"












files <- list.files("../results/divsums")
sp <- gsub("_", " ", gsub(".divsum$", "", files))
asmbsz <- df[!duplicated(df$species), ]
asmbsz <- asmbsz[asmbsz$species %in% sp, ]
asmbsz <- setNames(asmbsz$asmblysize.Mbp*1000000, asmbsz$species)
repstats <- data.frame()

for (i in 1:length(sp)) {
  species <- sp[i]
  # read text file into lines
  divsum <- readLines(paste0("../results/divsums/", files[i]))
  # look for the start of relevant information
  phrase <- "Coverage for each repeat class and divergence (Kimura)"
  start.index <- match(phrase, divsum) + 1
  # condense relevant lines into a table
  divsum <- divsum[start.index:length(divsum)]
  divsum <- read.table(textConnection(divsum), 
                       sep = " ", 
                       header = TRUE)
  # drop columns with all NA
  divsum <- divsum[-c(which(sapply(divsum, function(col) all(is.na(col)))))]
  divsum <- colSums(divsum)
  name <- c(name, species)
  cl <- c(cl, unique(df[df$species == species, ]$clade))
  r2 <- c(r2, unique(df[df$species == species, ]$rsq))
  largest <- c(largest, names(tail(sort(divsum), 1)))
}
df <- data.frame(name, cl, r2, largest)



library(beeswarm)
library(viridis)
cols <- viridis(4)
levels <- factor(as.factor(df$cl), levels = c("Others", "Actinopterygii", "Sauria", "Mammalia"))
par(mar = c(5, 4, 4, 7)+0.1)
beeswarm(r2 ~ largest, 
         data = df, 
         xlab = "largest repeat class", 
         ylab = "R2", 
         pch = 16, 
         pwcol = cols[levels])
text(8.02, (0.25+((0.8-0.25))), "Mammalia", xpd = NA, adj = 0)
text(8.02, (0.25+(2*(0.8-0.25)/3)), "Actinopterygii", xpd = NA, adj = 0)
text(8.02, (0.25+((0.8-0.25)/3)), "Sauria", xpd = NA, adj = 0)
text(8.02, 0.25, "Others", xpd = NA, adj = 0)
points(7.92, (0.25+((0.8-0.25))), pch = 16, col = "#FDE725FF", xpd = NA)
points(7.92, (0.25+(2*(0.8-0.25)/3)), pch = 16, col = "#31688EFF", xpd = NA)
points(7.92, (0.25+((0.8-0.25)/3)), pch = 16, col = "#35B779FF", xpd = NA)
points(7.92, 0.25, pch = 16, col = "#440154FF", xpd = NA)
par(mar = c(5.1, 4.1, 4.1, 2.1))


