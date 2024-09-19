


# constants
requiredPackages <- c("data.table", "ggplot2")

# paths
functionsPath <- "../../analysis/functions.R"
parsedResultsCsvPath <- "../../results/vertebrates/parsedResults.csv"
cladeDataCsvPath <- "../../data/vertebrates/cladeData.csv"
divsumDirPath <- "../../results/vertebrates/repeatLandscape"







# load stuff in
source(functionsPath)
loadPackages(requiredPackages)
parsedResults <- fread(parsedResultsCsvPath)
cladeData <- fread(cladeDataCsvPath)



getK2pMean <- function(files) {
  # read text file into lines
  divsumVector <- readLines(paste0(divsumDirPath, "/", files))
  # look for the start of useful information
  phrase <- "Coverage for each repeat class and divergence (Kimura)"
  startIndex <- match(phrase, divsumVector) + 1
  # condense the useful lines into a table
  divsumVector2 <- divsumVector[startIndex:length(divsumVector)]
  divsumTable <- read.table(textConnection(divsumVector2), sep = " ", header = TRUE)
  # drop NA columns
  divsumTable2 <- divsumTable[-c(which(sapply(divsumTable, 
                                              function(col) all(is.na(col)))))]
  divergence <- divsumTable$Div
  frequency <- rowSums(divsumTable2[, !names(divsumTable2) == "Div"])
  k2pMean <- sum(divergence*frequency)/sum(frequency)
  return(k2pMean)
}
getRsq <- function(species) {
  filterNames <- parsedResults$species == species
  if (any(filterNames)) {
    rsq <- unique(parsedResults$speciesRsquared[filterNames])
    return(rsq)
  } else {
    return(NA)
  }
}
getClass <- function(species) {
  filterNames <- cladeData$species == species
  if (any(filterNames)) {
    class <- cladeData$class[filterNames]
    return(class)
  } else {
    return(NA)
  }
}
files <- list.files(divsumDirPath)
k2pMean <- sapply(files, getK2pMean)
speciesLower <- gsub("_", " ", gsub("_summary\\.divsum$", "", files))
species <- gsub("^(\\w)(.*)", "\\U\\1\\L\\2", speciesLower, perl = TRUE)
rsq <- sapply(species, getRsq)
class <- sapply(species, getClass)
customClade <- class
customClade[customClade == "Actinopterygii"] <- "Ray-finned fish"
customClade[customClade == "Aves"] <- "Reptiles"
customClade[customClade == "Mammalia"] <- "Mammals"
customClade[customClade == "Reptilia"] <- "Reptiles"
otherClades <- !customClade %in% c("Ray-finned fish", "Reptiles", "Mammals")
customClade[otherClades] <- "Others"
customClade <- factor(customClade, levels = c("Mammals", "Ray-finned fish", "Reptiles", "Others"))
rsqVsK2pMean <- na.omit(data.frame(species, k2pMean, rsq, class, customClade))
fit <- lm(rsqVsK2pMean$rsq ~ rsqVsK2pMean$k2pMean)
slope <- signif(summary(fit)$coefficients[2, 1], 3)
intercept <- signif(summary(fit)$coefficients[1, 1], 3)
slopePvalue <- signif(summary(fit)$coefficients[2, 4], 3)
slopeRsquared <- signif(summary(fit)$adj.r.squared, 3)

ggplot(rsqVsK2pMean, aes(x = k2pMean, y = rsq, color = customClade)) +
  geom_point(shape = 16, alpha = 0.4, size = 2.3) +
  scale_color_manual(labels = c(
    paste0("Mammals\n(n = ", sum(customClade == "Mammals"), ")"),
    paste0("Ray-finned fish\n(n = ", sum(customClade == "Ray-finned fish"), ")"), 
    paste0("Reptiles\n(n = ", sum(customClade == "Reptiles"), ")"),
    paste0("Others\n(n = ", sum(customClade == "Others"), ")")
  ), values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3"))+
  ggtitle(bquote(italic(r)^2~"vs Estimated Genome Size"))+
  theme(plot.title = element_text(hjust = 0.475),
        axis.line = element_line(color = "black"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "#f2f2f2", color = "black", linewidth = 0.5),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "black", linetype = "dotted", size = 0.25),
        legend.position = c(0.86, 0.69),
        legend.key.size = unit(21, "points"))+
  xlim(c(10, 26)) +
  ylim(c(0, 1))+
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5, fullrange = TRUE)+
  labs(x = "Mean K2P Distance (subs per site)", y = bquote(italic(r)^2))
ggsave(filename = "rsqVsK2pMean.jpg", 
       plot = last_plot(), 
       width = 7680, 
       height = 4320, 
       units = "px", 
       dpi = 1100)

