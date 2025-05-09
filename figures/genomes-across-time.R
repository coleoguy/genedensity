

dat <- read.delim("../data/assembly.dates.tsv")
results <- read.csv("../results/masked.size.csv")
results$accession <- gsub("(^\\./|/.*$)", "", results$file)
results$assembly <- sub("^(?:[^_]*_){3}", "", results$file, perl = TRUE)
results$assembly <- sub("_genomic\\.fna$", "", results$assembly)
results <- results[, c("accession", "assembly", "lower", "N", "assem_sz")]
foo <- merge(results, dat, by.x = "accession", by.y = "Assembly.Accession")
foo$Assembly.Release.Date <- as.Date(foo$Assembly.Release.Date, format = "%Y-%m-%d")
foo$Organism.Name <- sub("^(([^ ]+ [^ ]+)).*$", "\\1", foo$Organism.Name)
foo <- foo[order(foo$Assembly.Release.Date, decreasing = TRUE), ]
foo <- foo[order(foo$Organism.Name), ]

for (i in 1:length(unique(foo$Organism.Name))) {
  sub <- foo[foo$Organism.Name == unique(foo$Organism.Name)[i], ]
  y <- sub$lower / sub$assem_sz
  plot(sub$Assembly.Release.Date, y, 
       xlab = "date", ylab = "masked/total", 
       main = unique(sub$Organism.Name))
  fit <- glm(y ~ sub$Assembly.Release.Date)
  abline(fit$coefficients)
  print(summary(fit))
  readline(prompt = "Press [Enter] to continue...")
}

