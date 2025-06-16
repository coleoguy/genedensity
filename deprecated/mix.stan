data {
  int<lower=1> N; // number of observations
  vector[N] x; // predictor (standardized)
  vector[N] y; // outcome (standardized)
}

parameters {
  ordered[2] alpha; // ordered intercepts: alpha[1] < alpha[2]
  real beta1; // slope for component 1
  real beta2; // slope for component 2
  real<lower=0> sigma1; // residual SD for component 1
  real<lower=0> sigma2; // residual SD for component 2
  real theta_logit; // mixture weight for component 1
}

transformed parameters {
  real<lower=0,upper=1> theta = inv_logit(theta_logit);
}

model {
  // priors
  alpha ~ normal(0, 0.5);
  beta1 ~ normal(1, 0.2);
  beta2 ~ normal(1, 0.2);
  sigma1 ~ student_t(3, 0, 1);
  sigma2 ~ student_t(3, 0, 1);
  theta_logit ~ normal(0, 1); // same as uniform on theta

  // log likelihood for mix
  for (n in 1:N) {
    target += log_mix(theta,
      normal_lpdf(y[n] | alpha[1] + beta1 * x[n], sigma1),
      normal_lpdf(y[n] | alpha[2] + beta2 * x[n], sigma2)
    );
  }
}

generated quantities {
  vector[N] y_rep;
  for (n in 1:N) { // avoid numerical instability when mixing
    if (bernoulli_rng(theta))
      y_rep[n] = normal_rng(alpha[1] + beta1 * x[n], sigma1);
    else
      y_rep[n] = normal_rng(alpha[2] + beta2 * x[n], sigma2);
  }
}


