
# Zhaobo Hu
# zhaobohu2002@gmail.com

# Description: Parses results and calculates additional statistics
# to summarize contigs for each species

input <- read.csv("../results/vertebrates/unparsed.csv")

# parsing
parsed <- input[input$size.Mbp >= 10, ]
sp.lessthanthree <- names(which(table(parsed$species) < 3))
parsed <- parsed[!(parsed$species %in% sp.lessthanthree), ]

# calculate stats
sp <- unique(input$species)
final <- data.frame()
for (species in sp) {
  sub <- input[input$species == species, ]
  sub2 <- parsed[parsed$species == species, ]
  if (nrow(sub2) > 0){
    fit <- summary(glm(sub2$genecount ~ sub2$size.Mbp))
    beta <- fit$coefficients[2, 1]
    pval.beta <- fit$coefficients[2, 4]
    rsq <- summary(lm(sub2$genecount ~ sub2$size.Mbp))$r.squared
    cv <- sd(sub2$genedens) / mean(sub2$genedens)
    weightmean <- sum(sub2$genedens * sub2$size.Mbp) / sum(sub2$size.Mbp)
    weightsd <- sqrt(sum(sub2$size.Mbp * (sub2$genedens - weightmean)^2) / sum(sub2$size.Mbp))
    weightcv <- weightsd / weightmean
    contig.stats <- data.frame(species, beta, pval.beta, rsq, cv, weightcv)
  } else {
    beta <- pval.beta <- rsq <- cv <- weightcv <- NA
    contig.stats <- data.frame(species, beta, pval.beta, rsq, cv, weightcv)
  }
  final <- rbind(final, merge(sub, contig.stats, by = "species"))
}

# write csv
write.csv(final, "../results/vertebrates/unparsed.csv", row.names = FALSE)
