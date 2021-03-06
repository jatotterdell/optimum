## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
	echo = FALSE,
	message = FALSE,
	warning = FALSE,
	comment = "#>",
  eval = TRUE # Change to TRUE if want to build vignette
)

## ------------------------------------------------------------------------
library(brms)
library(optimum)
library(mvtnorm)
library(ggplot2)
library(gridExtra)

## ------------------------------------------------------------------------
lpdf_pi <- function(p, mu, sigma) {
  -(log(p / (1 - p)) - mu)^2 / (2*sigma^2) - 0.5*log(2*pi) - log(p*(1 - p)) - log(sigma)
}
pdf_pi <- function(p, mu, sigma)
  exp(lpdf_pi(p, mu, sigma))

## ------------------------------------------------------------------------
parseq <- seq(0.001, 0.999, 0.001)
s <- c(0.5, 1, 1.78)
pdat <- data.frame(theta = rep(parseq, 3),
                   p = exp(c(sapply(seq(0.001, 0.999, 0.001), lpdf_pi, mu = 0, s = s[1]),
                     sapply(seq(0.001, 0.999, 0.001), lpdf_pi, mu = 0, s = s[2]),
                     sapply(seq(0.001, 0.999, 0.001), lpdf_pi, mu = 0, s = s[3]))),
           mu = rep(0, length(parseq)*3),
           s = rep(s, each = length(parseq)))
ggplot(pdat, aes(theta, p, colour = factor(s))) +
  geom_line() +
  scale_colour_discrete(expression(sigma))

## ------------------------------------------------------------------------
ggplot(data.frame(x = c(0, 1)), aes(x)) +
  stat_function(fun = pdf_pi, geom = "line", args = list(mu = 0, sigma = 1),
                aes(colour = "p0")) +
  stat_function(fun = pdf_pi, geom = "line", args = list(mu = 0, sigma = sqrt(1 + 1)),
                aes(colour = "p1")) +
  scale_color_discrete("Prior on", labels = c(expression(p[0]), expression(p[1]))) +
  labs(x = expression(theta), y = expression(f(theta))) +
  xlim(0,1) + ylim(0, 3)

## ------------------------------------------------------------------------
x <- c(0, 1)
X <- cbind(Intercept = 1, x = x)
n <- c(25, 25)
y <- c(12, 13)

dat <- data.frame(x = x, y = y, n = n)
fit_mcmc <- brm(y | trials(n) ~ x, data = dat, family = binomial, iter = 50000,
                prior = c(prior(normal(0,1.75), class = "b"), prior(normal(0,sqrt(10)), class = "Intercept")))
fit_vb <- vb_logistic_n(X, y, n, c(0, 0), diag(c(1.75^2, 10)), alg = "sj")

hist(posterior_samples(fit_mcmc)[, 1], freq = F, main = "VB vs. HMC posterior", breaks = 100,
     xlab = expression(beta[0]))
curve(dnorm(x, fit_vb$mu[1], fit_vb$Sigma[1,1]^(1/2)), add = T)

hist(posterior_samples(fit_mcmc)[, 2], freq = F, main = "VB vs. HMC posterior", breaks = 100,
     xlab = expression(beta[1]))
curve(dnorm(x, fit_vb$mu[2], fit_vb$Sigma[2,2]^(1/2)), add = T)

ypred_mcmc <- posterior_predict(fit_mcmc)

B <- rmvnorm(1e5, fit_vb$mu, fit_vb$Sigma)
ypred_vb <- cbind(rbinom(1e5, 25, plogis(B[, 1])), rbinom(1e5, 25, plogis(B[, 1] + B[, 2])))
ypred_m_vb <- plogis((X%*%fit_vb$mu) * 1 / sqrt(1 + pi*diag(fit_vb$Sigma)/8))


p_y_mcmc <- prop.table(table(factor(ypred_mcmc[, 1], levels = 0:25)))
p_y_vb <- prop.table(table(factor(ypred_vb[, 1], levels = 0:25)))
plot(p_y_vb - p_y_mcmc, main = "VB-PP minus HMC-PP")
barplot(p_y_mcmc, ylim = c(0, 0.15), ylab = "p_y_mcmc", main = "PP - HMC")
barplot(p_y_vb, ylim = c(0, 0.15), ylab = "p_y_vb", main = "PP - VB")

## ---- fig.height = 8-----------------------------------------------------
X <- cbind(Int = 1, x = c(0, 1))
n_foll <- seq(200, 3000, 400)/2
n_enro <- pmin(3000/2, n_foll + 1500/2)
n_miss <- n_enro - n_foll

# Beta-Binomial model
sim_dat <- sim_trial(p1tru = 0.1, p2tru = 0.1, n1int = n_foll, n2int = n_foll, 
                     a1 = 1, b1 = 1, a2 = 1, b2 = 1)
sim_dat <- calc_scenario_ppos(sim_dat, ppos_name = "final", post_method = "exact",
                              m1 = n_miss, m2 = n_miss,
                              l1 = max(n_foll) - n_foll, l2 = max(n_foll) - n_foll)

# Logistic model using VB
ptail_vb <- lapply(1:nrow(sim_dat), function(i) {
  fit <- vb_logistic_n(X, sim_dat[i, c(y1, y2)], sim_dat[i, c(n1, n2)], 
                       mu0 = c(0, 0), Sigma0 = diag(c(1.75, 1.75)), alg = "sj")
  beta <- rmvnorm(1e3, fit$mu, fit$Sigma)
  ypred <- data.table(y1 = rbinom(1e3, sim_dat[i, l1], plogis(beta[, 1])),
                      y2 = rbinom(1e3, sim_dat[i, l2], plogis(beta[, 1] + beta[, 2])))
  ypred <- ypred[, .N, keyby = .(y1, y2)]
  for(j in 1:nrow(ypred)) {
    pred_fit <- vb_logistic_n(X, 
                              sim_dat[i, c(y1, y2)] + ypred[j, c(y1, y2)],
                              sim_dat[i, c(n1, n2)] + sim_dat[i, c(l1, l2)],
                              mu0 = c(0, 0), Sigma0 = diag(c(10, 10)), alg = "sj")
    ypred[j, P := pnorm(0, pred_fit$mu[2], sqrt(pred_fit$Sigma[2,2]))]
  }
  ppos <- ypred[, sum(N*(P > 0.9)) / sum(N)]
  ptail <- pnorm(0, fit$mu[2], sqrt(fit$Sigma[2,2]))
  list(ptail, ppos, fit$mu, fit$Sigma)
})

m <- sapply(ptail_vb, `[[`, 3)
s <- sapply(ptail_vb, `[[`, 4)

## ---- fig.cap="At first analysis."---------------------------------------
hist(qlogis(rbeta(1e5, sim_dat[1, a1], sim_dat[1, b1])), freq = F, breaks = 100,
     ylim = c(0, 2.5), main = "Beta-binomial posterior", xlab = "Log-odds")
curve(dnorm(x, m[1, 1], sqrt(s[1, 1])), add = TRUE)
legend("topleft", lty = 1, legend = "VB approx", bty = 'n')

hist(qlogis(rbeta(1e5, sim_dat[1, a2], sim_dat[1, b2])) -
     qlogis(rbeta(1e5, sim_dat[1, a1], sim_dat[1, b1])), 
     freq = F, breaks = 100, ylim = c(0, 2.5), main = 'Beta-binomial posterior',
     xlab = "Difference in log-odds")
curve(dnorm(x, m[2, 1], sqrt(s[4, 1])), add = TRUE)
legend("topleft", lty = 1, legend = "VB approx", bty = 'n')

B <- rmvnorm(1e5, m[, 1], matrix(s[, 1], 2, 2))
hist(plogis(B[, 1]), freq = F, breaks = 100, 
     main = 'VB approx logistic model', xlab = bquote(theta[1]))
curve(dbeta(x, sim_dat[1, a1], sim_dat[1, b1]), add = TRUE)
legend("topright", lty = 1, legend = "Beta-binomial", bty = 'n')

hist(plogis(B[, 1] + B[, 2]), freq = F, breaks = 100, 
     main = 'VB approx logistic model', xlab = bquote(theta[2]))
curve(dbeta(x, sim_dat[1, a2], sim_dat[1, b2]), add = TRUE)
legend("topright", lty = 1, legend = "Beta-binomial posterior", bty = 'n')

## ---- fig.cap="At final analysis."---------------------------------------
hist(qlogis(rbeta(1e5, sim_dat[8, a1], sim_dat[8, b1])), freq = F, breaks = 100,
     main = "Beta-binomial posterior", xlab = "Log-odds")
curve(dnorm(x, m[1, 8], sqrt(s[1, 8])), add = TRUE)
legend("topleft", lty = 1, legend = "VB approx logistic model", bty = 'n')

hist(qlogis(rbeta(1e5, sim_dat[8, a2], sim_dat[8, b2])) -
     qlogis(rbeta(1e5, sim_dat[8, a1], sim_dat[8, b1])), 
     freq = F, breaks = 100, main = 'Beta-binomial posterior',
     xlab = "Difference in log-odds")
curve(dnorm(x, m[2, 8], sqrt(s[4, 8])), add = TRUE)
legend("topleft", lty = 1, legend = "VB approx logistic model", bty = 'n')

B <- rmvnorm(1e5, m[, 8], matrix(s[, 8], 2, 2))
hist(plogis(B[, 1]), freq = F, breaks = 100, 
     main = 'VB approx logistic model', xlab = bquote(theta[1]))
curve(dbeta(x, sim_dat[8, a1], sim_dat[8, b1]), add = TRUE)
legend("topright", lty = 1, legend = "Beta-binomial posterior", bty = 'n')

hist(plogis(B[, 1] + B[, 2]), freq = F, breaks = 100, 
     main = 'VB approx logistic model', xlab = bquote(theta[1]))
curve(dbeta(x, sim_dat[8, a2], sim_dat[8, b2]), add = TRUE)
legend("topright", lty = 1, legend = "Beta-binomial posterior", bty = 'n')

