
source("../../scripts/functions.R")
library(viridis)

significant <- read.csv("../../results/model-averaging.csv")
null <- read.csv("../../results/permute.csv")

titles <- c("All species", 
            "Mammals", 
            "Ray-finned fish", 
            "Reptiles")

par(
  mfrow = c(2, 2), 
  mar = c(2, 1, 2, 3), 
  oma = c(2, 4, 1, 7) 
)

var.cols <- c("steelblue", "firebrick")
datasets <- unique(null$dataset)

for (i in 1:length(datasets)) {
  
  # subset data
  null.id.cols <- c("dataset", "stat", "variable")
  null.result.cols <- setdiff(colnames(null), null.id.cols)
  cur.ds <- datasets[i]
  sub <- null[null$dataset == cur.ds, ]
  cur.vars <- significant[significant$clade == cur.ds, ]$model
  sub <- sub[sub$variable %in% cur.vars, ]
  
  title <- titles[i]
  
  # get data for each variable
  for (j in 1:length(cur.vars)) {
    
    # get kde of null beta
    cur.var <- cur.vars[j]
    cur.null <- sub[sub$variable == cur.var, ]
    upper <- unlist(cur.null[cur.null$stat == "upper", null.result.cols])
    lower <- unlist(cur.null[cur.null$stat == "lower", null.result.cols])
    null.beta <- (upper + lower)/2
    null.beta <- null.beta[!is.na(null.beta)]
    assign(paste0("null.beta.dens.", j), 
           density(null.beta))
    
    # modeled value
    modeled.variable <- significant[significant$model == cur.var 
                                    & significant$clade == cur.ds, ]
    assign(paste0("obs.beta.", j), 
           mean(as.numeric(modeled.variable[, c("lower", "upper")])))
    
    # confidence interval
    ci.temp <- quantile(null.beta, probs = c(0.025, 0.975))
    assign(paste0("beta.ci.", j), 
           ci.temp)
    assign(paste0("beta.to.shade.", j), 
           density(null.beta)$x >= ci.temp[1] & density(null.beta)$x <= ci.temp[2])
    
  }
  obs.beta <- ls(pattern = "^obs.beta\\.\\d+$")
  null.beta.dens <- ls(pattern = "^null.beta.dens\\.\\d+$")
  beta.ci <- ls(pattern = "^beta.ci\\.\\d+$")
  beta.to.shade <- ls(pattern = "^beta.to.shade\\.\\d+$")
  
  
  x.range <- range(unlist(lapply(obs.beta, get)), 
                   unlist(lapply(beta.ci, get)))
  x.range <- ifelse(x.range == max(x.range), 
                    x.range + abs(0.1 * x.range), 
                    x.range - abs(0.1 * x.range))
  y.max <- max(unlist(lapply(null.beta.dens, function(var) get(var)$y)))
  y.range <- c(0, 1.2 * y.max)
  
  plot(NA, 
       xlim = x.range, 
       ylim = y.range,
       xlab = NA, 
       ylab = NA, 
       yaxs = "i", 
       main = titles[i], 
       xaxt = "n", 
       yaxt = "n")
  axis(side = 1, 
       at = pretty(x.range, n = 4), 
       cex = 1.2, 
       tck = -0.05, 
       mgp = c(1, 0.6, 0))
  axis(side = 2, 
       at = pretty(y.range, n = 4), 
       cex = 1.2, 
       tck = -0.05, 
       mgp = c(1, 0.6, 0))
  
  for (k in 1:length(cur.vars)) {
    cur.obs <- get(paste0("obs.beta.", k))
    cur.dens <- get(paste0("null.beta.dens.", k))
    cur.ci <- get(paste0("beta.to.shade.", k))
    lines(cur.dens$x, cur.dens$y, col = var.cols[k], lwd = 0.5)
    polygon(c(cur.dens$x[cur.ci], rev(cur.dens$x[cur.ci])),
            c(cur.dens$y[cur.ci], rep(0, sum(cur.ci))),
            col = co(var.cols[k], 0.3), border = NA)
    abline(v = cur.obs, col = var.cols[k], lty = 3, lwd = 2)
  }
  rm(list = c(obs.beta, 
              null.beta.dens, 
              beta.ci, 
              beta.to.shade), envir = .GlobalEnv)
}

par(xpd = NA)



# axes
text(cx(0.46), cy(0.06), "Variable β", cex = 1.3, adj = c(0.5,0.5))
text(cx(0.035), cy(0.51), "Density", srt = 90, cex = 1.3, adj = c(0.5,0.5))


x <- 0.9
y1 <- 0.7
y2 <- 0.58
y3 <- 0.35

polygon(cx(c(x-0.05, x-0.05, x-0.036, x-0.036)), 
        cy(c(y3-0.023, y3+0.022, y3+0.022, y3-0.023)), 
        border = "steelblue", col = co("steelblue", 0.3), lwd = 1)
polygon(cx(c(x-0.034, x-0.034, x-0.02, x-0.02)), 
        cy(c(y3-0.023, y3+0.022, y3+0.022, y3-0.023)), 
        border = "firebrick", col = co("firebrick", 0.3), lwd = 1)


text(cx(x), cy(y3), "Variables", cex = 1.2, adj = c(0,0.5))
text(cx(x), cy(y1+0.02), "Null", cex = 1.2, adj = c(0,0.5))
text(cx(x), cy(y1-0.02), "expectation", cex = 1.2, adj = c(0,0.5))
text(cx(x), cy(y2+0.02), "Observed", cex = 1.2, adj = c(0,0.5))
text(cx(x), cy(y2-0.02), "value", cex = 1.2, adj = c(0,0.5))

segments(cx(x-0.02), cy(y1), cx(x-0.05), cy(y1), 
         col = "black", lwd = 2, lty = 1)

segments(cx(x-0.02), cy(y2), cx(x-0.05), cy(y2), 
         col = "black", lwd = 2, lty = 3)


