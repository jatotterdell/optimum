data {
  int<lower=0> N; // sample size
  int<lower=1> K; // parameters
  int y[N];    // responses
  int n[N];    // trials (all 1 if bernoulli)
  matrix[N, K] X;
}
parameters {
  vector[K] beta;
}
model {
  vector[N] eta = X * beta;
  beta ~ normal(0, 1);
  for(i in 1:N) {
    y[i] ~ binomial(n[i], inv_logit(eta[i]));
  }
}
