


dat <- read.csv("../results/rsq.csv")
files <- list.files("../results/divsums")
sp <- gsub(".divsum$", "", files)

for (i in 1:length(sp)) {
  species <- sp[i]
  asmbsz <- dat[dat$species == species, ]$assem.sz * 1000000
  if (length(asmbsz) == 0L) {
    next
  }
  
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
  l <- as.data.frame(divsum$Unknown)
  m <- divsum[, !(names(divsum) %in% c("Unknown", "Div"))]
  l$others <- rowSums(m)
  l <- l[order(l$`divsum$Unknown`), ]
  vec <- l$`divsum$Unknown`/(l$others + l$`divsum$Unknown`)
  vec[which(is.nan(vec))] <- 0
  plot(0:70, vec, 
       xlab = "div", ylab = "proportion of unknowns", 
       #ylim = c(0, 1), 
       pch = 16, cex = 0.8, 
       main = species)
  if (summary(glm(vec ~ c(0:70)))$coefficients[2, 4] < 0.05) {
    print(round(summary(glm(vec ~ c(0:70)))$coefficients[2, 1], 3))
  }
  readline("Press [Enter] to continue...")
}
