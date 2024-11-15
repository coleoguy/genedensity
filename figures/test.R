
vert.invert <- "vertebrates"

source("../analysis/functions.R")
dat <- read.csv("../results/vertebrates/final_results.csv")
rep <- read.csv("../results/vertebrates/all_repeat_landscapes.csv")

sp <- unique(rep$species)

files <- list.files("../results/vertebrates/repeat_landscape_divsums")
asmblysz <- read.csv("../results/vertebrates/assembly_sizes.csv")
asmblysz <- asmblysz[asmblysz$species %in% sp, ]
asmblysz <- asmblysz[order(asmblysz$species == sp), ]

dat <- data.frame(asmblysz, files)


calcRepLandscStats <- function(i, file, asmblysz.Mbp, vert.invert) {
  # read text file into lines
  divsum.vector <- readLines(
    paste0("../results/", vert.invert, "/repeat_landscape_divsums/", file))
  # look for the start of relevant information
  phrase <- "Coverage for each repeat class and divergence (Kimura)"
  start.index <- match(phrase, divsum.vector) + 1
  # condense relevant lines into a table
  divsum.vector <- divsum.vector[start.index:length(divsum.vector)]
  divsum.table <- read.table(textConnection(divsum.vector), 
                             sep = " ", 
                             header = TRUE)
  # drop NA columns
  divsum.table <- divsum.table[
    -c(which(sapply(divsum.table, function(col) all(is.na(col)))))]
  # vector of divergence scores
  divergence <- divsum.table$Div
  # vector of the frequencies of each divergence score
  frequency <- rowSums(divsum.table[, !names(divsum.table) == "Div"])
  # repeat content in Mbp
  repcontent.Mbp <- sum(frequency) / 1000000
  # the bin that contains the median
  median.bin <- which(cumsum(frequency) > sum(frequency)/2)[1]
  # frequency of the previous bin
  lower <- cumsum(frequency)[median.bin-1]
  # frequency of the next bin
  upper <- cumsum(frequency)[median.bin+1]
  # median frequency
  mid <- sum(frequency)/2
  # median bin
  k2p.median <- median.bin + (mid-lower)/(upper-lower)
  # calculate mean
  k2p.mean <- sum(divergence*frequency)/sum(frequency)
  # repeat content in percent coverage
  repcontent.percentcoverage <- (repcontent.Mbp / asmblysz.Mbp) * 100
  
  percentfreq <- ((frequency / 1000000) / asmblysz.Mbp) * 100
  
  smoothfreq <- smooth.spline(1:length(percentfreq), percentfreq, spar = 0.6)$y
  
  peakpos <- findpeaks(smoothfreq, minpeakdistance = 5)[, 2]
  if (smoothfreq[1] > smoothfreq[2]) {
    peakpos <- c(1, peakpos)
  }
  peakpos <- min(peakpos)
  deriv <- diff(smoothfreq)
  upper <- sd(deriv) * 0.1
  lower <- -sd(deriv) * 0.1
  peakstart <- which(deriv > upper)[1]
  if (peakpos == 1) {
    peakstart <- 1
  }
  end.cand <- which(deriv < lower)[which(deriv < lower) > peakpos]
  for (k in seq_along(end.cand)[-length(end.cand)]) {
    if (end.cand[k + 1] - end.cand[k] > 1) {
      peakend <- end.cand[k]
      break
    }
  }
  peakend <- end.cand[k]
  width <- peakend - peakstart
  peakarea <- sum(smoothfreq[peakstart:peakend])
  peakweightedarea <- sum(smoothfreq[peakstart:peakend] * seq(width, 0, by = -1))
  
  
  
  
  
  # build dataframe
  df <- data.frame(i, 
                   repcontent.Mbp,
                   repcontent.percentcoverage, 
                   k2p.mean, 
                   k2p.median,
                   peakpos,
                   peakstart,
                   peakend,
                   peakarea,
                   peakweightedarea)
  return(df)
}

replandsc.stats <- data.frame()

for (i in sp) {
  file <- dat[dat$species == i, ]$files
  asmblysz.Mbp <- dat[dat$species == i, ]$asmbly.size.Mbp
  result <- lapply(i, 
                   calcRepLandscStats, 
                   file = file,
                   asmblysz.Mbp = asmblysz.Mbp,
                   vert.invert = vert.invert)[[1]]
  replandsc.stats <- rbind(replandsc.stats, result)
}
colnames(replandsc.stats)[colnames(replandsc.stats) == "i"] <- "species"


chromnums <- read.csv(paste0("../data/", vert.invert, "/chromnums.csv"))
assembly.sizes <- read.csv(paste0("../results/", vert.invert, "/assembly_sizes.csv"))
contig.stats <- read.csv(paste0("../results/", vert.invert, "/contig_stats.csv"))
clades <- read.csv(paste0("../results/", vert.invert, "/clades.csv"))
taxo.gnsz <- read.csv(paste0("../data/", vert.invert, "/taxo_gnsz.csv"))

a <- merge(chromnums, assembly.sizes, by = "species", all.x = TRUE)
b <- merge(a, contig.stats, by = "species", all.x = TRUE)
c <- merge(b, clades, by = "species", all.x = TRUE)
d <- merge(c, taxo.gnsz, by = "species", all.x = TRUE)
final.results <- merge(d, replandsc.stats, by = "species", all.x = TRUE)


































packages <- c("ape", "ggplot2", "nlme")
lapply(packages, library, character.only = TRUE)
tree <- read.tree(paste0("../data/", vert.invert, "/formatted_tree.nwk"))
tree$tip.label <- gsub("_", " ", tree$tip.label)

# gather and subset relevant results
final.result <- final.results[!is.na(final.results$chromnum.1n), ]
dat <- na.omit(final.results[, c("species", "weightcv", "peakarea", "clade")])
sp.intersect <- intersect(dat$species, tree$tip.label)
dat <- dat[dat$species %in% sp.intersect, ]

# prune tree
pruned.tree <- keep.tip(tree, sp.intersect)

# create PGLS object for trendline
pgls.model <- gls(weightcv ~ peakarea, 
                  data = dat, 
                  correlation = corBrownian(phy = pruned.tree, form = ~species),
                  method = "ML")
summary <- summary(pgls.model)
slope <- signif(summary$tTable[2, 1], 3)
intercept <- signif(summary$tTable[1, 1], 3)
slope.pval <- signif(summary$tTable[2, 4], 3)

# calculate PICs for permutation test of pearson correlation coefficient
y <- pic(setNames(dat$weightcv, dat$species), pruned.tree)
x <- pic(setNames(dat$peakarea, dat$species), pruned.tree)
perm.pval <- signif(permTest(x, y, 100000, "pearson"), 3)

# set factors for figure legend
dat$clade <- factor(dat$clade, levels = c("Mammalia", "Actinopterygii", "Sauria", "Others"))

# graph
ggplot(dat, aes(x = peakarea, y = weightcv, color = clade)) +
  geom_point(shape = 16, alpha = 0.4, size = 2.3) +
  scale_color_manual(labels = c(
    paste0("Mammals\n(n = ", sum(dat$clade == "Mammalia"), ")"),
    paste0("Ray-finned fish\n(n = ", sum(dat$clade == "Actinopterygii"), ")"), 
    paste0("Reptiles\n(n = ", sum(dat$clade == "Sauria"), ")"),
    paste0("Others\n(n = ", sum(dat$clade == "Others"), ")")
  ), values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3"))+
  theme(plot.title = element_text(hjust = 0.475),
        plot.subtitle = element_text(hjust = 0.475), 
        axis.line = element_line(color = "black"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "#f2f2f2", color = "black", linewidth = 0.5),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "black", linetype = "dotted", size = 0.25),
        legend.position = c(0.86, 0.69),
        legend.key.size = unit(21, "points"))+
  geom_abline(intercept = intercept, slope = slope, color = "black", linetype = "dashed", linewidth = 0.5) +
  labs(title = bquote("Weighted CV vs peakarea"), 
       subtitle = bquote(italic(β) * "-coefficient" == .(slope) * "," ~~ italic(β) ~ italic(p) * "-value" == .(slope.pval) * "," ~~ "permutation" ~ italic(p) * "-value" == .(perm.pval)),
       x = "peakarea", 
       y = bquote("Weighted CV"))


