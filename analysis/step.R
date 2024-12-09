
# model
dat <- read.csv("../results/vertebrates/parsed.csv")
dat <- dat[!duplicated(dat$species), ]

d <- na.omit(dat[, c("species", "rsq", "clade", "median", "totalrep.pct")])
summary(step(glm(d$rsq ~ d$median * d$totalrep.pct)))

m <- d[d$clade == "Mammalia", ]
summary(step(glm(m$rsq ~ m$median * m$totalrep.pct)))

f <- d[d$clade == "Actinopterygii", ]
summary(step(glm(f$rsq ~ f$median * f$totalrep.pct)))

r <- d[d$clade == "Sauria", ]
summary(step(glm(r$rsq ~ r$median * r$totalrep.pct)))

# test for phylogenetic signal
library(phytools)
tree <- read.tree("../data/vertebrates/formatted_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
int <- intersect(d$species, tree$tip.label)
pruned.tree <- keep.tip(tree, int)
res <- setNames(resid(step(glm(d$rsq ~ d$median:d$totalrep.pct))), d$species)
phylosig(pruned.tree, res, method = "lambda", test = TRUE)

# PGLS model
library(caper)
cd <- comparative.data(tree, d, species)
pgls.model <- pgls(rsq ~ median + totalrep.pct + median : totalrep.pct, data = cd)
summary(pgls.model)

# compare AICs
d <- d[d$species %in% int, ]
summary(gls(rsq ~ median + totalrep.pct + median:totalrep.pct, 
    data = d, 
    correlation = corBrownian(phy = pruned.tree, form = ~species),
    method = "ML")) # 21.78158
summary(gls(rsq ~ median, 
            data = d, 
            correlation = corBrownian(phy = pruned.tree, form = ~species),
            method = "ML")) # 23.40565
summary(gls(rsq ~ totalrep.pct, 
            data = d, 
            correlation = corBrownian(phy = pruned.tree, form = ~species),
            method = "ML")) # 25.3864

