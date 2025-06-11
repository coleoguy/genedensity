library(rstan)
library(posterior)
source("functions.R")

set.seed(1234)
num.cores <- round(0.8 * parallel::detectCores())
options(mc.cores = num.cores)

rsq <- read.csv("../results/rsq.csv")
repeats <- read.csv("../results/repeat-results.csv")
dat <- merge(rsq, repeats, by.x = "species", by.y = "species", all.x = T, all.y = T)
results <- read.csv("../results/model-averaging.csv")

# for each clade...
for (i in 1:length(unique(results$clade))) {
  
  # subset data by clade
  clade <- unique(results$clade)[i]
  if (clade == "All") {
    sub <- dat
  } else {
    sub <- dat[dat$clade == clade, ]
  }
  dat.vars <- colnames(sub)[grep("^(prop|age)\\.", colnames(sub))]
  sub <- na.omit(sub[, c("species", "clade", "rsq", dat.vars)])
  
  # subset data from significant variables
  est.vars <- results[results$clade == clade, ]$model
  dat.new <- as.data.frame(matrix(NA, nrow(sub), 0))
  for (k in 1:length(est.vars)) {
    est.var <- est.vars[k]
    if (grepl(":", est.var)) {
      strsplit(est.var, ":")
      dat.cols <- sub[, unlist(strsplit(est.var, ":"))]
      dat.new <- cbind(dat.new, apply(dat.cols, 1, prod))
      names(dat.new)[k] <- est.var
    } else {
      dat.cols <- sub[, est.var]
      dat.new <- cbind(dat.new, dat.cols)
      names(dat.new)[k] <- est.var
    }
  }
  row.names(dat.new) <- sub$species
  
  # rescale data
  for (j in 1:length(dat.new)) {
    col <- names(dat.new)[j]
    dat.new[[col]] <- (dat.new[[col]]-min(dat.new[[col]])) / diff(range(dat.new[[col]]))
  }
  
  # simulate R2 using estimated predictors
  coefs <- results[results$clade == clade, ]$estimate
  dat.new <- as.matrix(dat.new)
  r2.new <- dat.new %*% coefs
  r2.new <- setNames(as.vector(r2.new), row.names(r2.new))
  
  # match simulated data with observed data
  r2.observed <- na.omit(setNames(dat$rsq, dat$species))
  int <- intersect(names(r2.new), names(r2.observed))
  r2.new <- r2.new[names(r2.new) %in% int]
  r2.observed <- r2.observed[names(r2.observed) %in% int]
  
  # kde for observed r2
  obs.dens <- density(r2.observed)
  
  # standardize predicted and observed r2
  mean.x <- mean(r2.new); sd.x <- sd(r2.new)
  mean.y <- mean(r2.observed); sd.y <- sd(r2.observed)
  r2.new      <- (r2.new - mean.x) / sd.x
  r2.observed <- (r2.observed - mean.y) / sd.y
  
  # assemble dataframe for stan
  df <- list(
    N = length(r2.new),
    x = r2.new,
    y = r2.observed
  )
  
  # fit
  fit.bayes <- stan(
    file = "mix.stan",
    data = df,
    chains = 4, 
    cores = num.cores, 
    iter = 8000,
    warmup = 3000,
    control = list(adapt_delta = 0.995, max_treedepth = 15),
    seed = 1234
  )
  
  # posterior predictive backtransform
  post.std <- extract(fit.bayes)$y_rep
  post.bt <- sweep(sweep(post.std, 2, sd.y, "*"), 2, mean.y, "+")
  xcommon <- seq(
    min(c(r2.observed, post.bt)),
    max(c(r2.observed, post.bt)),
    length.out = 512
  )
  
  # posterior predictive mean
  pp.mat <- sapply(1:nrow(post.bt), function(i)
    approx(density(post.bt[i,])$x, density(post.bt[i,])$y, xout = xcommon, rule = 2)$y
  )
  pp.mean <- rowMeans(pp.mat) # find posterior predictive mean
  
  # posterior predictive check
  plot(0,0, type = "n", xlim = range(obs.dens$x),
       ylim = c(0, 1.2 * max(obs.dens$y, pp.mean)),
       xlab = "observed R2", ylab = "density",
       main = clade)
  for (l in sample(nrow(post.bt), 30)) lines(density(post.bt[l,]), col = co("wheat", 0.6))
  lines(obs.dens, col = "firebrick", lwd = 2)
  lines(xcommon, pp.mean, col = "steelblue", lwd = 2, lty = 2)
  
  # 5) posterior of beta1 & beta2
  post.df  <- as_draws_df(fit.bayes, variable = c("beta1", "beta2"))
  dens1    <- density(post.df$beta1)
  dens2    <- density(post.df$beta2)
  ci1      <- quantile(post.df$beta1, probs = c(0.025, 0.975))
  ci2      <- quantile(post.df$beta2, probs = c(0.025, 0.975))
  xlim.all <- c(-1, 2)
  ylim.all <- c(0, 1.2 * max(dens1$y, dens2$y))
  plot(0, 0, type = "n", xlim = xlim.all, ylim = ylim.all,
       xlab = "slope", ylab = "density",
       main = clade)
  lines(dens1, col = co("firebrick", 0.5), lwd = 2)
  idx1 <- which(dens1$x >= ci1[1] & dens1$x <= ci1[2])
  polygon(
    c(dens1$x[idx1], rev(dens1$x[idx1])),
    c(dens1$y[idx1], rep(0, length(idx1))),
    col = co("firebrick", 0.2), border = NA
  )
  lines(dens2, col = co("steelblue", 0.5), lwd = 2)
  idx2 <- which(dens2$x >= ci2[1] & dens2$x <= ci2[2])
  polygon(
    c(dens2$x[idx2], rev(dens2$x[idx2])),
    c(dens2$y[idx2], rep(0, length(idx2))),
    col = co("steelblue", 0.2), border = NA
  )
  abline(v = 0, col = "black", lty = "73", lwd = 2)
  abline(v = 1, col = "black", lty = "73", lwd = 2)
}






















############################## brms snippet #################################

# center
r2.observed <- r2.observed - median(r2.observed)
r2.new <- r2.new - median(r2.new)

df <- data.frame(
  r2.new,
  r2.observed, 
  obs = seq_along(r2.new)
)

bf <- bf(
  r2.observed ~ 0 + beta * r2.new,
  beta ~ 1 + (1 | obs),
  nl = TRUE
)

priors <- c(
  prior(normal(1, 0.5), class="b", nlpar="beta"),    # slope mean
  prior(student_t(3, 0, 1), class = "sigma"),   # residual SD
  prior(student_t(3, 0, 1), class = "sd", nlpar = "beta")  # slope SD
)

fit0 <- brm(
  bf, 
  data = df, 
  prior = priors,
  chains = 6, 
  iter = 5000, 
  warmup = 2000, 
  family = gaussian(link = "identity"),
  cores = 15, 
  control = list(adapt_delta = 0.99), 
  seed = 1234
)

# posterior predictive check
yrep <- posterior_predict(fit0, draws = 1000)
dens <- apply(yrep, 1, density)
pp.x <- seq(-1, 1, length.out = 512)
pp.matrix <- sapply(dens, function(d) approx(d$x, d$y, xout = pp.x, rule = 2)$y)
pp.mean <- rowMeans(pp.matrix)
plot(0, 0, type = "n",
     xlim = c(-1, 1), 
     ylim = c(0, max(c(pp.mean, density(df$r2.observed)$y))),
     main = "Posterior predictive check",
     xlab = "observed R2", ylab = "density")
for (l in sample(1:nrow(yrep), 30)) {
  lines(density(yrep[l, ]), col = co("wheat", 0.8))}
lines(density(df$r2.observed), col = "firebrick", lwd = 2)
lines(pp.x, pp.mean, col = "steelblue", lwd = 2, lty = 2)

# posterior distribution of slope
post <- as_draws_df(fit0, variable = "b_beta_Intercept")
dens <- density(post$b_beta_Intercept)
ci <- quantile(post$b_beta_Intercept, probs = c(0.025, 0.975))
plot(0, 0, type = "n", xlim = range(dens$x), ylim = c(0, max(dens$y)),
     xlab = "beta (observed ~ simulated R2)", ylab = "density", main = "posterior distribution")
lines(dens, col = "firebrick", lwd = 2)
idx <- which(dens$x >= ci[1] & dens$x <= ci[2])
polygon(c(dens$x[idx], rev(dens$x[idx])),
        c(dens$y[idx], rep(0, length(idx))),
        col = co("firebrick", 0.3), border = NA)
abline(v = 1, col = "steelblue", lty = 3, lwd = 2)
abline(v = 0, col = "steelblue", lty = 3, lwd = 2)

########################### end brms snippet ##############################












############################### eb snippet ##################################

# Prepare data
x_obs <- r2.new
y_obs <- r2.observed
df <- data.frame(
  R2_obs       = x_obs,
  R2_pred_base = y_obs
)

# Fast loss function for prior predictive fit
loss_for_fast <- function(int_mean, beta_mean, sigma_scale,
                          x_obs, y_obs, n_draws = 1000) {
  alpha_draws <- rnorm(n_draws, int_mean, 0.15)
  beta_draws  <- rnorm(n_draws, beta_mean, 0.15)
  sigma_draws <- pmax(abs(rt(n_draws, df = 3) * sigma_scale / sqrt(3)), 1e-4)
  
  ysim <- vapply(seq_len(n_draws), function(j) {
    alpha_draws[j] + beta_draws[j] * x_obs +
      rnorm(length(x_obs), 0, sigma_draws[j])
  }, numeric(length(x_obs)))
  
  if (!all(is.finite(ysim))) return(1e9)
  
  # For each x_obs[i], compare simulated y values to y_obs[i]
  loss_vec <- sapply(seq_along(y_obs), function(i) {
    yi_sim <- ysim[, i]
    
    mean_sq_error     <- (mean(yi_sim) - y_obs[i])^2
    quant_error       <- mean((quantile(yi_sim, probs = c(0.1, 0.5, 0.9)) - y_obs[i])^2)
    wasserstein_error <- mean(abs(sort(yi_sim) - y_obs[i]))
    
    # Weighted total loss
    1.0 * mean_sq_error +
      0.5 * quant_error +
      0.5 * wasserstein_error
  })
  
  mean(loss_vec)
}



scoringFunction <- function(mu_int, beta_mean, sigma_scale) {
  tryCatch({
    score <- -loss_for_fast(mu_int, beta_mean, sigma_scale, x_obs, y_obs, n_draws = 500)
    list(Score = max(score, -1e9))
  }, error = function(e) {
    # Always return the same structure, even on error
    list(Score = -1e9)  # Keep just 1 named element
  })
}


# Optimization bounds (no more conservative clipping)
bounds_list <- list(
  mu_int      = c(0, 1),
  beta_mean   = c(0, 1.5),
  sigma_scale = c(0, 1.5)
)

# Run ParBayesianOptimization

# Register a cluster of 15 workers
cl <- makeCluster(15)
clusterExport(cl, varlist = c("loss_for_fast", "x_obs", "y_obs"))
registerDoParallel(cl)
opt_res <- bayesOpt(
  FUN = scoringFunction,
  bounds = bounds_list,
  initPoints = 20,
  iters.n = 100,
  gsPoints = 1, 
  acq = "ei", 
  parallel = TRUE, 
  verbose = 1
)
stopCluster(cl)
registerDoSEQ()

# Extract best parameters
best <- getBestPars(opt_res)
mu_int      <- best["mu_int"]
beta_mean   <- best["beta_mean"]
sigma_scale <- best["sigma_scale"]

# Define priors
priors_auto <- c(
  set_prior(paste0("normal(", mu_int, ",0.15)"), class = "Intercept"),
  set_prior(paste0("normal(", beta_mean, ",0.15)"), class = "b", coef = "R2_obs"),
  set_prior(paste0("normal(0,", sigma_scale, ")"), class="sigma")
  
)

# Fit brms model (prior-only)
fit_clade <- brm(
  R2_pred_base ~ 1 + R2_obs,
  data         = df,
  prior        = priors_auto,
  sample_prior = "only",
  chains       = 2,
  iter         = 2000,
  cores        = 15,
  seed         = 1234
)

# Prior-predictive plot
ppd <- posterior_predict(fit_clade, ndraws = 500)
x <- x_obs
cols <- rgb(0, 0, 1, 0.05)

plot(range(x), range(ppd, y_obs),
     type = "n", xlab = "Simulated R²", ylab = "Observed R² (prior)",
     main = "Prior-Predictive (BO-tuned priors)")
for (l in seq_len(nrow(ppd))) {
  points(x, ppd[l, ], pch = 16, col = cols)
}
abline(0, 1, lty = 2, col = "gray")
points(x, y_obs, pch = 19, col = "red")

# Hypothesis tests
print(hypothesis(fit_clade, "Intercept = 0"))
print(hypothesis(fit_clade, "R2_obs = 1"))

# Plot posterior densities
post <- as_draws_df(fit_clade, variable = c("Intercept", "b_R2_obs"))
dens_int   <- density(post$Intercept)
dens_slope <- density(post$b_R2_obs)

ci_int   <- quantile(post$Intercept, probs = c(0.025, 0.975))
ci_slope <- quantile(post$b_R2_obs, probs = c(0.025, 0.975))

xlim_all <- range(dens_int$x, dens_slope$x)
ylim_all <- range(dens_int$y, dens_slope$y)

plot(0, 0, type = "n", xlim = xlim_all, ylim = ylim_all,
     xlab = "Value", ylab = "Density", main = "Posterior Distributions")
lines(dens_int, col = "steelblue", lwd = 2)
lines(dens_slope, col = "firebrick", lwd = 2, lty = 2)

ix_int <- which(dens_int$x >= ci_int[1] & dens_int$x <= ci_int[2])
polygon(c(dens_int$x[ix_int], rev(dens_int$x[ix_int])),
        c(dens_int$y[ix_int], rep(0, length(ix_int))),
        col = rgb(70,130,180, max = 255, alpha = 80), border = NA)

ix_slope <- which(dens_slope$x >= ci_slope[1] & dens_slope$x <= ci_slope[2])
polygon(c(dens_slope$x[ix_slope], rev(dens_slope$x[ix_slope])),
        c(dens_slope$y[ix_slope], rep(0, length(ix_slope))),
        col = rgb(178,34,34, max = 255, alpha = 80), border = NA)

abline(v = 0, col = "steelblue", lty = 3)
abline(v = 1, col = "firebrick", lty = 3)

legend("topright",
       legend = c("Intercept", "Slope (b_R2_obs)"),
       col = c("steelblue", "firebrick"),
       lty = c(1, 2), lwd = 2,
       bty = "n")

############################### end eb snippet ################################





