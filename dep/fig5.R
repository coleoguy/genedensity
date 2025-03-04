

terms <- c(
  "intercept", 
  "median.trans", 
  "rep.prop", 
  "median.trans:rep.prop"
)
dat <- read.csv("../results/models.csv")
names(dat) <- c("clade", "thrs", "repeat", "phylosig", "n", "stat", terms, "f2")


# models where every term is either significant or NA
df <- dat[apply(dat[, terms[-1]], 1, function(x) all(is.na(x) | x < 0.05)), ]

# remove models with all NA terms
df <- df[rowSums(is.na(df[, terms[-1]])) < length(terms[-1]), ]

# filter for significant clade-threshold-repeat combinations
df <- df[df$stat == "p", ]
hit <- c("clade", "thrs", "repeat")
df <- df[, hit]

# get beta coefficients for significant combinations
df <- merge(dat, df, by = intersect(names(dat), names(df)))
df <- df[df$stat == "beta", ]
df <- df[!df$`repeat` %in% "rc", ]
df <- df[order(df$thrs), ]
df <- df[order(df$`repeat`), ]
df <- df[order(df$clade), ]
rownames(df) <- c(1:nrow(df))
df <- df[df$thrs >= 0.8 & df$thrs <= 0.95, ]


#combinations of clades and repeats
comb <- unique(df[, c("clade", "repeat"), ])
beta <- c()
for (i in 1:nrow(comb)) {
  cl <- comb[i, ][[1]]
  rep <- comb[i, ][[2]]
  sub <- df[df$clade == cl, ]
  sub <- sub[sub$`repeat` == rep, ]
  sub <- sub[, c(8, 9)]
  beta <- c(beta, mean(as.matrix(sub)[which(!is.na(sub))]))
}
comb$beta <- beta
comb$abs <- abs(beta)

# define tick marks
n <- 5  
log.max <- log10(max(comb$abs) + 1)  
ticks <- seq(0, log.max, length.out = n)
labels <- 10^ticks - 1  

# color mapping
library(viridis)
dat.trans <- log10(comb$abs + 1)
norm <- dat.trans / log.max  
res <- 1000  
palette <- viridis(res)  
comb$cols <- palette[round(norm * (res - 1)) + 1]  
