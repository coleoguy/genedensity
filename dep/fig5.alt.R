


terms <- c(
  "intercept", 
  "median.trans", 
  "rep.prop", 
  "median.trans:rep.prop"
)
dat <- read.csv("../results/models1.csv")
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
df <- df[df$thrs >= 0.8, ]

#combinations of clades and repeats
comb <- unique(df[, c("clade", "repeat"), ])
beta <- c()
for (i in 1:nrow(comb)) {
  cl <- comb[i, ][[1]]
  rep <- comb[i, ][[2]]
  sub <- df[df$clade == cl, ]
  sub <- sub[sub$`repeat` == rep, ]
  sub <- sub[, terms[-1]]
  to.keep <- c()
  for (j in 1:nrow(sub)) {
    if (sum(is.na(sub[j, terms[-1]])) == length(terms[-1])-1) {
      to.keep <- c(to.keep, j)
    }
  }
  sub <- sub[to.keep, ]
  if (nrow(sub) == 0) {
    beta <- c(beta, NA)
    next
  }
  sub <- sub[, terms[-c(1, 4)]]
  beta <- c(beta, median(as.matrix(sub)[which(!is.na(sub))]))
}
comb$beta <- beta
comb$abs <- abs(beta)
comb <- na.omit(comb)

# define tick marks
n <- 5
seq(from = 0, to = max(comb$abs), length.out = 5)

# color mapping
library(viridis)
res <- 1000  
palette <- viridis(res)  
comb$cols <- palette[round((comb$abs / (max(comb$abs))) * (res-1)) + 1]





















# define tick marks
n <- 5
seq(from = min(comb$beta), to = max(comb$beta), length.out = 5)

# color mapping
library(viridis)
library(RColorBrewer)
res <- 1000  
palette <- colorRampPalette(c(brewer.pal(8, "RdYlBu")))(res) 
# palette <- brewer.pal(11, "PiYG")
comb$cols <- palette[round(((comb$beta - min(comb$beta)) / diff(range(comb$beta))) * (res-1)) + 1]


# Create a matrix representing the color bar
z <- matrix(1:res, nrow = 1)
# Plot the color bar without axes or labels
image(z, col = palette, axes = FALSE, xlab = "", ylab = "")






