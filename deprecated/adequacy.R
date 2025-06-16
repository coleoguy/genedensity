
source("../../scripts/functions.R")
library(viridis)

significant <- read.csv("../../results/model-averaging.csv")
null <- read.csv("../../results/permute.csv")

titles <- c("LTR age x LTR proportion", 
            "SINE age x SINE proportion", 
            "SINE age x SINE proportion", 
            "Other repeats proportion", 
            "LTR age", 
            "Unknown repeats age")

par(
  mfrow = c(2, 3), 
  mar   = c(3, 5, 2, 1), 
  oma   = c(3, 2, 1, 9) 
)

clade.cols <- setNames(c("black", "#d95f02", "#7570b3", "#1b9e77"), 
                       unique(significant$clade))

for (i in 1:nrow(significant)) {
  null.id.cols <- c("dataset", "stat", "variable")
  null.result.cols <- setdiff(colnames(null), null.id.cols)
  cur.variable <- significant[i, ]
  sub <- null[null$dataset == cur.variable$clade, ]
  sub <- sub[sub$variable == cur.variable$model, ]
  null.impo <- as.numeric(sub[sub$stat == "importance", null.result.cols])
  null.impo <- null.impo[!is.na(null.impo)]
  null.lower <- as.numeric(sub[sub$stat == "lower", null.result.cols])
  null.upper <- as.numeric(sub[sub$stat == "upper", null.result.cols])
  null.beta <- (null.upper + null.lower) / 2
  null.beta <- null.beta[!is.na(null.beta)]
  
  title <- titles[i]
  
  null.beta.dens <- density(null.beta)
  cur.beta <- mean(as.numeric(cur.variable[, c("lower", "upper")]))
  
  null.impo.dens <- density(null.impo)
  cur.impo <- cur.variable$importance
  
  beta.ci <- quantile(null.beta, probs = c(0.025, 0.975))
  beta.to.shade <- null.beta.dens$x >= beta.ci[1] & null.beta.dens$x <= beta.ci[2]
  
  impo.ci <- quantile(null.impo, probs = c(0.025, 0.975))
  impo.to.shade <- null.impo.dens$x >= impo.ci[1] & null.impo.dens$x <= impo.ci[2]
  
  x.range <- range(cur.beta, # modeled beta
                   beta.ci, # null beta
                   cur.impo, # modeled importance
                   impo.ci) # null importance
  x.range <- ifelse(x.range == max(x.range), 
                    x.range + abs(0.1 * x.range), 
                    x.range - abs(0.1 * x.range))
  y.range <- c(0, 1.1 * max(null.beta.dens$y))
 
  clade.col <- clade.cols[match(cur.variable$clade, names(clade.cols))]
  
  # beta
  plot(null.beta.dens, 
       xlim = x.range, 
       ylim = y.range, 
       xlab = NA, 
       ylab = NA, 
       main = title, 
       col.main = clade.col, 
       cex.main = 1, 
       col = "steelblue", 
       yaxs     = "i" 
       )
  box(col = clade.col, lwd = 0.8)
  polygon(c(null.beta.dens$x[beta.to.shade], rev(null.beta.dens$x[beta.to.shade])),
          c(null.beta.dens$y[beta.to.shade], rep(0, sum(beta.to.shade))),
          col = co("steelblue", 0.4), border = NA)
  abline(v = cur.beta, col = "steelblue", lty = 3, lwd = 2)
  
  # importance
  lines(null.impo.dens$x, null.impo.dens$y, col = "firebrick", lwd = 0.5)
  polygon(c(null.impo.dens$x[impo.to.shade], rev(null.impo.dens$x[impo.to.shade])),
          c(null.impo.dens$y[impo.to.shade], rep(0, sum(impo.to.shade))),
          col = co("firebrick", 0.4), border = NA)
  abline(v = cur.variable$importance, col = "firebrick", lty = 3, lwd = 2)
}

par(xpd = NA)

# axes
text(-2.65, -1.4, "Value", cex = 1.8, adj = c(0.5,0.5))
text(-7.3, 4.1, "Likelihood", srt = 90, cex = 1.8, adj = c(0.5,0.5))

y <- 8
x <- 2.1

# legend
text(x, y, "All", cex = 1, adj = c(0,0.5))
text(x, y-0.3, "species", cex = 1, adj = c(0,0.5))
segments(x-0.55, y-0.15, x-0.2, y-0.15, col = clade.cols[1], lwd = 3,)

text(x, y-0.9, "Mammals", cex = 1, adj = c(0,0.5))
segments(x-0.55, y-0.9, x-0.2, y-0.9, col = clade.cols[2], lwd = 3,)

text(x, y-1.5, "Ray-finned", cex = 1, adj = c(0,0.5))
text(x, y-1.8, "fish", cex = 1, adj = c(0,0.5))
segments(x-0.55, y-1.65, x-0.2, y-1.65, col = clade.cols[3], lwd = 3,)

text(x, y-2.4, "Reptiles", cex = 1, adj = c(0,0.5))
segments(x-0.55, y-2.4, x-0.2, y-2.4, col = clade.cols[4], lwd = 3,)

text(x, y-4, "Null", cex = 1, adj = c(0,0.5))
text(x, y-4.3, "beta", cex = 1, adj = c(0,0.5))
segments(x-0.55, y-4.15, x-0.2, y-4.15, col = "steelblue", lwd = 3,)

text(x, y-4.9, "Observed", cex = 1, adj = c(0,0.5))
text(x, y-5.2, "beta", cex = 1, adj = c(0,0.5))
segments(x-0.55, y-5.05, x-0.2, y-5.05, col = "steelblue", lty = 3, lwd = 3)

text(x, y-5.8, "Null", cex = 1, adj = c(0,0.5))
text(x, y-6.1, "variable", cex = 1, adj = c(0,0.5))
text(x, y-6.4, "importance", cex = 1, adj = c(0,0.5))
segments(x-0.55, y-6.1, x-0.2, y-6.1, col = "firebrick", lwd = 3,)

text(x, y-7, "Observed", cex = 1, adj = c(0,0.5))
text(x, y-7.3, "variable", cex = 1, adj = c(0,0.5))
text(x, y-7.6, "importance", cex = 1, adj = c(0,0.5))
segments(x-0.55, y-7.3, x-0.2, y-7.3, col = "firebrick", lty = 3, lwd = 3)

par(xpd = FALSE)
