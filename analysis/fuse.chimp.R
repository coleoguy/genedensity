
# make a "human" genome by fusing chimpanzee chromosomes 12 & 13 and
# compare gene density variation between chimp and human

# r squared became more like human but not by much

dat <- read.csv("../results/parsed.csv")
c <- dat[dat$species == "Pan troglodytes", ]
h <- dat[dat$species == "Homo sapiens", ]
sub <- c[c$name %in% c("NC_072410.2", "NC_072411.2"), ]

size <- sum(sub$size.Mbp)
genecount <- sum(sub$genecount)
genedens <- size / genecount
new <- sub[1, ]
new$size.Mbp <- size
new$genecount <- genecount
new$genedens <- genedens

c <- c[!(c$name %in% c("NC_072410.2", "NC_072411.2")), ]
c <- rbind(c, new)
c$rsq <- summary(lm(c$genecount ~ c$size.Mbp))$r.squared


