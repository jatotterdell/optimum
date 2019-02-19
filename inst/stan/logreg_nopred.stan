data {
  int N;
  int P;
  int y[N];
  int n[N];
  matrix[N, P] X;
}

parameters {
  vector[P] beta;
}

transformed parameters {
  real p[N];
  for (i in 1:N) {
    p[i] = inv_logit(X[i, ] * beta);
  }
}

model {
  beta ~ normal(0, 1);
  y ~ binomial(n, p);
}
