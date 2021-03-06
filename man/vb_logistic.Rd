% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{vb_logistic}
\alias{vb_logistic}
\title{Perform variational inference for logistic regression model}
\usage{
vb_logistic(X, y, mu0, Sigma0, tol = 1e-08, maxiter = 1000L,
  maxiter_jj = 25L, alg = "jj")
}
\arguments{
\item{X}{The design matrix}

\item{y}{The response vector}

\item{mu0}{The prior mean for beta paramter}

\item{Sigma0}{The prior variance for beta parameter}

\item{tol}{The tolerance level to assess convergence}

\item{maxiter}{The maximum number of iterations}

\item{maxiter_jj}{The maximum number of Jaakkola-Jordan iterations to initialise estimation}

\item{alg}{The algorithm used for final estimation of variational parameters. 
Must be one of "jj", "sj", "kmw".}
}
\description{
Perform variational inference for logistic regression model
}
