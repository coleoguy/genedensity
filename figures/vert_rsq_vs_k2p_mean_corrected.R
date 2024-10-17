


packages <- c("ape", "ggplot2", "nlme")
lapply(packages, library, character.only = TRUE)
source("../analysis/functions.R")
vert.invert <- "vertebrates"
tree <- read.tree("../data/vertebrates/chordates_species.nwk")
parsed.results <- read.csv("../results/vertebrates/parsed_results.csv")
clades.gnsz <- read.csv("../data/vertebrates/clades_gnsz.csv")
files <- list.files("../results/vertebrates/repeat_landscape")
k2p.mean <- sapply(files, getK2pMean)
species.lower <- gsub("_", " ", gsub("_summary\\.divsum$", "", files))
species <- gsub("^(\\w)(.*)", "\\U\\1\\L\\2", species.lower, perl = TRUE)
rsq <- sapply(species, getRsq)
class <- sapply(species, getClass)
custom.clade <- class
custom.clade[custom.clade == "Actinopterygii"] <- "Ray-finned fish"
custom.clade[custom.clade == "Aves"] <- "Reptiles"
custom.clade[custom.clade == "Mammalia"] <- "Mammals"
custom.clade[custom.clade == "Reptilia"] <- "Reptiles"
other.clades <- !custom.clade %in% c("Ray-finned fish", "Reptiles", "Mammals")
custom.clade[other.clades] <- "Others"
custom.clade <- factor(custom.clade, levels = c("Mammals", "Ray-finned fish", "Reptiles", "Others"))
dat <- na.omit(data.frame(species, k2p.mean, rsq, custom.clade))


# prune and format tree
sp <- dat$species
spf <- sub("^([^_]*_[^_]*)_.*", "\\1", gsub(" ", "_", sp))
sp.intersect <- intersect(tree$tip.label, spf)
pruned.tree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% sp.intersect)])
spn <- sp[match(spf[match(pruned.tree$tip.label, spf)], spf)]
pruned.tree$tip.label <- spn

# subset and reorder dataframe based on pruned tree data
dat <- dat[dat$species %in% spn, ]
dat <- dat[match(spn, dat$species), ]

pgls.model <- gls(rsq ~ k2p.mean, 
                       data = dat[1:3], 
                       correlation = corBrownian(phy = pruned.tree, form = ~species),
                       method = "ML")
summary <- summary(pgls.model)


intercept <- signif(summary$tTable[1, 1], 3)
slope <- signif(summary$tTable[2, 1], 3)
pval <- signif(summary$tTable[2, 4], 3)

ggplot(dat, aes(x = k2p.mean, y = rsq, color = custom.clade)) +
  geom_point(shape = 16, alpha = 0.4, size = 2.3) +
  scale_color_manual(labels = c(
    paste0("Mammals\n(n = ", sum(dat$custom.clade == "Mammals"), ")"),
    paste0("Ray-finned fish\n(n = ", sum(dat$custom.clade == "Ray-finned fish"), ")"), 
    paste0("Reptiles\n(n = ", sum(dat$custom.clade == "Reptiles"), ")"),
    paste0("Amphibians\n(n = ", sum(dat$custom.clade == "Others"), ")")
  ), values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3"))+
  theme(plot.title = element_text(hjust = 0.475), 
        plot.subtitle = element_text(hjust = 0.475), 
        axis.line = element_line(color = "black"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "#f2f2f2", color = "black", linewidth = 0.5),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "black", linetype = "dotted", size = 0.25),
        legend.position = c(0.86, 0.65),
        legend.key.size = unit(21, "points"))+
  geom_abline(intercept = intercept, slope = slope, color = "black", linetype = "dashed", linewidth = 0.5) +
  xlim(c(8, 27)) +
  labs(title = bquote(italic(r)^2 ~ "vs Estimated Genome Size"), 
       subtitle = bquote(italic(β) * "-coefficient" == .(slope) * "," ~~ italic(p) * "-value" == .(pval)),
       x = "Estimated Genome Size (Gbp)", 
       y = bquote(italic(r)^2))
ggsave(filename = "vert_rsq_vs_k2p_mean_corrected.jpg", 
       plot = last_plot(), 
       width = 7680, 
       height = 4320, 
       units = "px", 
       dpi = 1100)

























rsq <- replicate(3, sample(dat$rsq, nrow(dat), replace = TRUE))
rsq.means <- apply(rsq, 2, mean)
rsq.vars <- apply(rsq, 2, var)
rsq.mean.int <- t.test(rsq.means, conf.level = 0.99)$conf.int
rsq.var.int <- t.test(rsq.vars, conf.level = 0.99)$conf.int

k2p.mean <- replicate(3, sample(dat$k2p.mean, nrow(dat), replace = TRUE))
k2p.mean.means <- apply(k2p.mean, 2, mean)
k2p.mean.vars <- apply(k2p.mean, 2, var)
k2p.mean.mean.int <- t.test(k2p.mean.means, conf.level = 0.99)$conf.int
k2p.mean.var.int <- t.test(k2p.mean.vars, conf.level = 0.99)$conf.int




# Number of bootstrap iterations
n_boot <- 1000
x <- dat[, 2]
y <- dat[, 3]
calc_slope <- function(x, y) {
  fit <- lm(y ~ x)
  return(coef(fit)[2])
}
n <- length(x)
boot_indices <- matrix(sample(n, size = n * n_boot, replace = TRUE), nrow = n, ncol = n_boot)
boot_slopes <- apply(boot_indices, 2, function(idx) calc_slope(x[idx], y[idx]))

# Calculate statistics
mean_slope <- mean(boot_slopes)
se_slope <- sd(boot_slopes)
ci_slope <- quantile(boot_slopes, c(0.025, 0.975))

# Print results
cat("Mean bootstrapped slope:", mean_slope, "\n")
cat("Standard error of bootstrapped slope:", se_slope, "\n")
cat("95% CI of slope: [", ci_slope[1], ",", ci_slope[2], "]\n")

# Plot histogram of bootstrap slopes
hist(boot_slopes, main = "Bootstrap Distribution of Slope", xlab = "Slope")

# Add original slope to plot
abline(v = calc_slope(x, y), col = "red", lwd = 2)

# Add mean bootstrapped slope to plot
abline(v = mean_slope, col = "blue", lwd = 2, lty = 2)

# Add legend
legend("topright", legend = c("Original Slope", "Mean Bootstrap Slope"),
       col = c("red", "blue"), lwd = 2, lty = c(1, 2))





