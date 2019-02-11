// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rmath.h>

using namespace Rcpp;

//' Evaluate standard normal cdf for matrix of variates
//' 
//' @param m A matrix of variates 
// [[Rcpp::export]]
arma::mat pnorm_mat(arma::mat& m) {
  int p = m.n_cols;
  int n = m.n_rows;
  arma::mat out(n, p);
  
  for (int i = 0; i < n; i++) {
    for(int j = 0; j < p; j++) {
      out(i, j) = R::pnorm(m(i, j), 0.0, 1.0, 1, 0);
    }
  }
  return out;
}

//' Evaluate standard normal density for matrix of variates
//' 
//' @param m A matrix of variates 
// [[Rcpp::export]]
arma::mat dnorm_mat(arma::mat& m) {
  int p = m.n_cols;
  int n = m.n_rows;
  arma::mat out(n, p);
  
  for (int i = 0; i < n; i++) {
    for(int j = 0; j < p; j++) {
      out(i, j) = R::dnorm(m(i, j), 0.0, 1.0, 0);
    }
  }
  return out;
}

//' Perform Jaakkola-Jordan update of variational parameters
//' 
//' @param X The design matrix
//' @param y The response vector
//' @param eta1 The current value of 1st natural parameter
//' @param eta2 The current value of the 2nd natural parameter
//' @param eta1_p The prior value of the 1st natural parameter
//' @param eta2_p The prior value of the 2nd natural parameter
// [[Rcpp::export]]
double jaakkola_jordan(
  const arma::mat& X, const arma::vec& y, 
  arma::vec& eta1, arma::vec& eta2, 
  const arma::vec& eta1_p, const arma::vec& eta2_p
) {
  
  int p = X.n_cols;
  arma::mat Sigma0 = -0.5*inv(reshape(eta2_p, p, p));
  arma::vec mu0 = Sigma0 * eta1_p;
  arma::mat tmp1 = inv(reshape(eta2 + eta2_p, p, p));
  arma::mat tmp2 = (eta1 + eta1_p) * trans(eta1 + eta1_p);
  arma::mat Xi = 0.25*tmp1*(tmp2*tmp1 - 2*arma::eye<arma::mat>(p,p));
  arma::vec xi = sqrt(diagvec(X*Xi*trans(X)));
  arma::mat K = trans(X) * diagmat( tanh(xi / 2) / (4 * xi) ) * X;
  eta1 = trans(X) * (y - 0.5);
  eta2 = -vectorise(K);
  arma::mat Sigma = -0.5*inv(reshape(eta2 + eta2_p, p, p));
  arma::vec mu = Sigma * (eta1 + eta1_p);

  double l = real(0.5*log_det(Sigma) - 0.5*log_det(Sigma0)) +
    as_scalar(0.5*trans(mu)*inv(Sigma)*mu - 0.5*trans(mu0)*inv(Sigma0)*mu0) +
    sum(xi / 2 - log(1 + exp(xi)) + (xi / 4) % tanh(xi / 2));
  return l;
};

//' Perform Saul-Jordan update of variational parameters
//' 
//' @param X The design matrix
//' @param y The response vector
//' @param eta1 The current value of 1st natural parameter
//' @param eta2 The current value of the 2nd natural parameter
//' @param eta1_p The prior value of the 1st natural parameter
//' @param eta2_p The prior value of the 2nd natural parameter
//' @param omega1 The current value of the Omega1 variational parameter
// [[Rcpp::export]]
double saul_jordan(
    const arma::mat& X, const arma::vec& y, 
    arma::vec& eta1, arma::vec& eta2, 
    const arma::vec& eta1_p, const arma::vec& eta2_p,
    arma::vec& omega1
) {
  
  int p = X.n_cols;
  arma::mat Sigma0 = -0.5*inv(reshape(eta2_p, p, p));
  arma::vec mu0 = Sigma0 * eta1_p;
  arma::mat tmp1 = -0.5*X*inv(reshape(eta2 +eta2_p, p, p));
  arma::vec mu = tmp1 * (eta1 + eta1_p);
  arma::vec sigma2 = diagvec(tmp1 * trans(X));
  arma::vec omega0 = mu + 0.5 * (1 - 2 * omega1) % sigma2;
  omega1 = 1/(1 + exp(-omega0));
  arma::vec omega2 = 1 / (2.0 * (1 + cosh(omega0)));
  eta1 = trans(X) * (y - omega1 + omega2 % mu);
  eta2 = -0.5*vectorise(trans(X) * diagmat(omega2) * X);
  
  arma::mat Sigma = -0.5*inv(reshape(eta2 + eta2_p, p, p));
  mu = Sigma * (eta1 + eta1_p);
  
  double l = real(0.5*log_det(Sigma) - 0.5*log_det(Sigma0)) + 0.5*p -
    0.5*trace( inv(Sigma0) * (Sigma + (mu - mu0) * trans(mu - mu0)) ) +
    dot(y, X * mu) - 0.5 * dot(omega1 % omega1, sigma2) -
    sum(log(1 + exp(X*mu + 0.5 * (1 - 2 * omega1) % sigma2)));
  return l;
}

//' Perform Knowles-Minka-Wand update of variational parameters
//' 
//' @param X The design matrix
//' @param y The response vector
//' @param eta1 The current value of 1st natural parameter
//' @param eta2 The current value of the 2nd natural parameter
//' @param eta1_p The prior value of the 1st natural parameter
//' @param eta2_p The prior value of the 2nd natural parameter
//' @param MS_p The 
//' @param MS_s The
// [[Rcpp::export]]
double knowles_minka_wand(
    const arma::mat& X, const arma::vec& y, 
    arma::vec& eta1, arma::vec& eta2, 
    const arma::vec& eta1_p, const arma::vec& eta2_p,
    const arma::vec& MS_p, const arma::vec& MS_s
) {
  
  int p = X.n_cols;
  arma::mat Sigma0 = -0.5*inv(reshape(eta2_p, p, p));
  arma::vec mu0 = Sigma0 * eta1_p;
  
  arma::mat tmp1 = -0.5*X*inv(reshape(eta2 + eta2_p, p, p));
  arma::vec m = tmp1 * (eta1 + eta1_p);
  arma::vec sigma2 = diagvec(tmp1 * trans(X));
  arma::mat Omega = sqrt(1 + sigma2 * trans(MS_s % MS_s));
  arma::mat tmp2 = (m * trans(MS_s)) / Omega;
  arma::vec omega3 = pnorm_mat(tmp2) * MS_p;
  arma::vec omega4 = (dnorm_mat(tmp2) / Omega) * (MS_p % MS_s);
  eta1 = trans(X) * (y - omega3 + omega4 % m);
  eta2 = -0.5 * vectorise(trans(X) * diagmat(omega4) * X);
  
  arma::mat Sigma = -0.5*inv(reshape(eta2 + eta2_p, p, p));
  arma::vec mu = Sigma * (eta1 + eta1_p);
  
  double l = real(0.5*log_det(Sigma) - 0.5*log_det(Sigma0)) + 0.5*p -
    0.5*trace( inv(Sigma0) * (Sigma + (mu - mu0) * trans(mu - mu0)) ) +
    dot(y, X * mu) - sum(omega3);
  return l;
}

//' Perform variational inference for logistic regression model
//' 
//' @param X The design matrix
//' @param y The response vector
//' @param mu0 The prior mean for beta paramter
//' @param Sigma0 The prior variance for beta parameter
//' @param tol The tolerance level to assess convergence
//' @param maxiter The maximum number of iterations
//' @param maxiter_jj The maximum number of Jaakkola-Jordan iterations to initialise estimation
//' @param string_alg The algorithm used for final estimation of variational parameters. 
//' Must be one of "jj", "sj", "kmw".
// [[Rcpp::export]]
List vb_logistic(
  const arma::mat& X, const arma::vec& y,
  const arma::vec& mu0, const arma::mat& Sigma0,
  double tol = 1e-8, int maxiter = 1000, int maxiter_jj = 25,
  std::string alg = "jj"
) {
  
  int p = X.n_cols;
  int n = X.n_rows;
  
  bool converged = 0;
  bool jj_converged = 0;
  int iterations = 0;
  arma::mat trace(maxiter, p);
  arma::vec elbo(maxiter);
  
  arma::vec eta2_p = -0.5*arma::vectorise(inv(Sigma0));
  arma::vec eta1_p = solve(Sigma0, mu0);
  arma::vec eta1 = arma::zeros(p);
  arma::vec eta2 = -0.5*arma::vectorise(arma::diagmat(arma::ones(p)));
  
  arma::vec omega1 = 0.5 + arma::zeros(n);
  
  
  arma::vec MS_p;
  MS_p << 0.003246343272134 << arma::endr
       << 0.051517477033972 << arma::endr
       << 0.195077912673858 << arma::endr
       << 0.315569823632818 << arma::endr
       << 0.274149576158423 << arma::endr
       << 0.131076880695470 << arma::endr
       << 0.027912418727972 << arma::endr
       << 0.001449567805354 << arma::endr;
  arma::vec MS_s;
  MS_s << 1.365340806296348 << arma::endr
       << 1.059523971016916 << arma::endr
       << 0.830791313765644 << arma::endr
       << 0.650732166639391 << arma::endr
       << 0.508135425366489 << arma::endr
       << 0.396313345166341 << arma::endr
       << 0.308904252267995 << arma::endr
       << 0.238212616409306 << arma::endr;
  
  for(int i = 0; i < maxiter && !converged; i++) {
    if(!jj_converged && i < maxiter_jj) {
      elbo(i) = jaakkola_jordan(X, y, eta1, eta2, eta1_p, eta2_p);
    } 
    else if(alg == "jj") {
      elbo(i) = jaakkola_jordan(X, y, eta1, eta2, eta1_p, eta2_p);
    } else if(alg == "sj") {
      elbo(i) = saul_jordan(X, y, eta1, eta2, eta1_p, eta2_p, omega1);
    } else if (alg == "kmw") {
      elbo(i) = knowles_minka_wand(X, y, eta1, eta2, eta1_p, eta2_p, MS_p, MS_s);
    }

    // check for convergence
    if(i > 0 && !jj_converged && i < maxiter_jj && fabs(elbo(i) - elbo(i - 1)) < tol) {
      jj_converged = 1;
    } else if(i > 0 && fabs(elbo(i) - elbo(i - 1)) < tol) {
      converged = 1;
      iterations = i;
    }
  }
  
  arma::mat Sigma = -0.5*inv(reshape(eta2 + eta2_p, p, p));
  arma::vec mu = Sigma * (eta1 + eta1_p);
  
  return List::create(Named("converged") = converged,
                      Named("jj_converged") = jj_converged,
                      Named("elbo") = elbo.subvec(0, iterations),
                      Named("mu") = mu,
                      Named("Sigma") = Sigma);
}
