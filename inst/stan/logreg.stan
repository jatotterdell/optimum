data {
  int N;
  int P;
  int y[N];
  int n[N];
  matrix[N, P] X;
  
  int Npred;
  int npred[Npred];
  matrix[Npred, P] Xpred;
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
  beta ~ normal(0, 10);
  y ~ binomial(n, p);
}

generated quantities {
  int ypred[Npred];
  for(i in 1:Npred) {
    ypred[i] = binomial_rng(npred[i], inv_logit(Xpred[i, ] * beta));
  }
}
