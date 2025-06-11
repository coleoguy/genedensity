data {
  int<lower=1> N;
  vector[N]    x;  // standardized predictor
  vector[N]    y;  // standardized outcome
}
parameters {
  real          beta0;     // pop-level slope on std scale
  vector[N]     u_raw;     // unit-normal deviations
  real<lower=0> sigma_u;   // SD of those deviations
  real<lower=0> sigma;     // residual SD on std scale
}
transformed parameters {
  vector[N] u = u_raw * sigma_u;  // random slopes on std scale
}
model {
  // weakly informative priors on std scale
  beta0 ~ normal(0, 2);
  sigma_u ~ normal(0, 2) T[0, ];
  u_raw   ~ normal(0, 1);
  sigma ~ normal(0, 5) T[0, ];

  // likelihood on std scale
  y ~ normal((beta0 + u) .* x, sigma);
}
generated quantities {
  vector[N] y_rep;
  for (n in 1:N)
    y_rep[n] = normal_rng((beta0 + u[n]) * x[n], sigma);
}
