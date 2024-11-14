


# Zhaobo Hu
# zhaobohu2002@gmail.com

# Description: for each species, fit a linear model for parsed contigs and 
# calculate slope, p-value of slope, r-squared, and coefficient of varience

# "vertebrates" or "invertebrates"
vert.invert <- "vertebrates"

# read results
input <- read.csv(paste0("../results/", 
                           vert.invert, 
                           "/all_contigs_results.csv"))
all.species <- unique(input$species)

contig.stats <- data.frame()
for (species in all.species) {
  sub <- input[input$species == species, ]
  fit <- summary(glm(sub$contig.gene.count ~ sub$contig.size.Mbp))
  beta <- fit$coefficients[2, 1]
  pval.beta <- fit$coefficients[2, 4]
  rsq <- summary(lm(sub$contig.gene.count ~ sub$contig.size.Mbp))$r.squared
  cv <- sd(sub$contig.genedens.geneperMbp) / mean(sub$contig.genedens.geneperMbp)
  weightmean <- sum(sub$contig.genedens.geneperMbp * sub$contig.size.Mbp) / sum(sub$contig.size.Mbp)
  weightsd <- sqrt(sum(sub$contig.size.Mbp * (sub$contig.genedens.geneperMbp - weightmean)^2) / sum(sub$contig.size.Mbp))
  weightcv <- weightsd / weightmean
  contig.stats <- rbind(contig.stats, data.frame(species, beta, pval.beta, rsq, cv, weightcv))
}

# write csv
write.csv(contig.stats, 
          paste0("../results/", vert.invert, "/contig_stats.csv"), 
          row.names = FALSE)



