# library(mystan)
# library(rstan)
# library(brms)
# library(data.table)
# library(optimum)
# 
# n <- 50
# X <- cbind(x0 = 1, x1 = c(rep(0, n/2), rep(1, n / 2)))
# y <- sample(c(0, 1), n, replace = T)
# d1 <- cbind(X, y)
# d2 <- data.table(d1)[, .(y = sum(y), n = .N), by = .(x0, x1)]
# mod <- stan_model("inst/stan/logreg_nopred.stan")
# mod_pred <- stan_model("inst/stan/logreg.stan")
# m <- sampling(mod, data = list(X = as.matrix(d2[, 1:2]), y = d2[,y], n = d2[, n],
#                                N = 2, P = 2))
# 
# eta1 <- c(0,0)
# eta2 <- c(-0.5,0,0,-0.5)
# optimum:::jaakkola_jordan(X, y, eta1, eta2, c(0,0), c(-0.5,0,0,-0.5))
# optimum:::jaakkola_jordan_n(as.matrix(d2[,1:2]), as.matrix(d2[, 3]), as.matrix(d2[, 4]),
#                             eta1, eta2, c(0,0), c(-0.5,0,0,-0.5))
# 
# omega1 <- rep(0.5, length(y))
# optimum:::saul_jordan(X, y, eta1, eta2, c(0,0), c(-0.5,0,0,-0.5), omega1)
# omega1 <- rep(0.5, nrow(d2))
# optimum:::saul_jordan_n(as.matrix(d2[,1:2]), as.matrix(d2[, 3]), as.matrix(d2[, 4]),
#                             eta1, eta2, c(0,0), c(-0.5,0,0,-0.5), omega1)
# 
# o1 <- optimum:::vb_logistic(X, y, c(0,0), diag(2), alg = "jj")
# o2 <- optimum:::vb_logistic_n(as.matrix(d2[,1:2]), as.matrix(d2[,3]), as.matrix(d2[,4]), c(0,0), diag(2), alg = "jj")
# 
# o1 <- optimum:::vb_logistic(X, y, c(0,0), diag(2), alg = "sj")
# o2 <- optimum:::vb_logistic_n(as.matrix(d2[,1:2]), as.matrix(d2[,3]), as.matrix(d2[,4]), c(0,0), diag(2), alg = "sj")
# 
# hist(extract(m)$beta[, 1], freq = F)
# curve(dnorm(x, o2$mu[1], sqrt(o2$Sigma[1,1])), add = T)
# 
# hist(extract(m)$beta[, 2], freq = F)
# curve(dnorm(x, o2$mu[2], sqrt(o2$Sigma[2,2])), add = T)
# 
# 
# 
# 
# lpdf_pi <- function(p, mu, sigma) {
#   -(-mu + log(p / (1 - p)))^2 / (2*sigma^2) - 0.5*log(2*pi) - log(p - p^2) - log(sigma)
# }
# curve(exp(lpdf_pi(x, 0, 10)), xlim = c(0, 1), ylim = c(0, 2))
# curve(exp(lpdf_pi(x, 0, 1.75)), xlim = c(0, 1), ylim = c(0, 2))
# curve(dbeta(x, 1, 1), add = T, col = "red")
# 
# sims <- 100
# sim_dat <- optimum:::sim_scenario(sims, p1tru = 0.23, p2tru = 0.2, n1int = c(500, 550), n2int = c(500, 550), post_method = "exact")
# sim_ppos <- optimum:::calc_scenario_ppos(sim_dat, post_method = "exact")
# 
# P <- outer(0:50, 0:50, function(x, y)
#   dbetabinom(x, 50, sim_dat[1, ]$a1, sim_dat[1, ]$b1) * dbetabinom(y, 50, sim_dat[1, ]$a2, sim_dat[1, ]$b2))
# Post <- outer(0:50, 0:50, function(x,y)
#   Vectorize(beta_ineq)(sim_dat[1, ]$a1 + x, sim_dat[1, ]$b1 + 50 - x,
#                        sim_dat[1, ]$a2 + y, sim_dat[1, ]$b2 + 50 - y))
# 
# vb <- vb_logistic_n(cbind(1, c(0, 1)), c(sim_dat[1, y1], sim_dat[1, y2]), c(500, 500),mu0 = c(0, 0), Sigma = diag(2), alg = "sj")
# m <- logreg_stan_pred(cbind(1, c(0, 1)), c(sim_dat[1, y1], sim_dat[1, y2]), c(500, 500), X, c(50, 50), chains = 1, iter = 1000)
# m <- sampling(mod_pred, dat = list(X = cbind(1, c(0, 1)), y = c(sim_dat[1, y1], sim_dat[1, y2]), n = c(500, 500), 
#                                    Xpred = cbind(1, c(0, 1)), npred = c(50, 50), Npred = 2, N = 2, P = 2), chains = 1, iter = 1000)
# 
# Pvb <- pnorm(0, vb$mu[2], sqrt(vb$Sigma[2,2]))
# Ptail <- mean(extract(m)$beta[, 2] < 0)
# ypred <- data.table(extract(m)$ypred)[, .N, by = .(V1, V2)]
# for(i in 1:nrow(ypred)) {
#   tmp <- logreg_stan(X,
#                      sim_dat[sim_id == 1 & stage == 1, c(y1, y2)] + ypred[i, c(V1, V2)],
#                      c(550, 550),
#                      chains = 4, iter = 1000, refresh = 0)
#   ypred[i, P := mean(extract(tmp)$beta[, 2] < 0)]
# }
# c(sum((Post > 0.9) * P), ypred[, sum(N*(P > 0.9)) / sum(N)])
# 
# mean(plogis(extract(m)$beta[, 1]))
# mean(plogis(apply(extract(m)$beta[, 1:2], 1, sum)))
# mean(extract(m)$beta[, 2] < 0)
# hist(plogis(extract(m)$beta[, 1]), freq = FALSE, breaks = 50)
# curve(dbeta(x, sim_dat[1, a1], sim_dat[1, b1]), xlim = c(0.1, 0.3), add = T)
# 
# hist(plogis(extract(m)$beta[, 1] + extract(m)$beta[, 2]), freq = FALSE, breaks = 50)
# curve(dbeta(x, sim_dat[1, a2], sim_dat[1, b2]), xlim = c(0.1, 0.3), add = T)
# 
# PPoS <- rep(0, sims)
# Ptail <- rep(0, sims)
# 
# pt <- proc.time()
# for(s in 1:sims) {
#   mod1 <- logreg_stan_pred(X,
#                            sim_dat[sim_id == s & stage == 1, c(y1, y2)],
#                            sim_dat[sim_id == s & stage == 1, c(n1, n2)],
#                            X,
#                            sim_ppos[sim_id == s & stage == 1, c(ppos_m1, ppos_m2)],
#                            chains = 1, iter = 2000, refresh = 0)
#   Ptail[s] <- mean(extract(mod1)$beta[, 2] < 0)
#   ypred <- data.table(extract(mod1)$ypred)[, .N, by = .(V1, V2)]
#   for(i in 1:nrow(ypred)) {
#     tmp <- logreg_stan(X,
#                        sim_dat[sim_id == s & stage == 1, c(y1, y2)] + ypred[i, c(V1, V2)],
#                        sim_dat[sim_id == s & stage == 1, c(n1, n2)] + sim_ppos[sim_id == s & stage == 1, c(ppos_m1, ppos_m2)],
#                        chains = 1, iter = 2000, refresh = 0)
#     ypred[i, P := mean(extract(tmp)$beta[, 2] < 0)]
#   }
#   PPoS[s] <- ypred[, sum(N*(P > 0.9)) / sum(N)]
# }
# proc.time() - pt
# 
# library(doParallel)
# registerDoParallel(cores = 4)
# pt <- proc.time()
# Ppar <- foreach(s = 1:sims, .packages = c("mystan", "optimum", "rstan", "data.table"), .combine = rbind) %dopar% {
#   mod1 <- logreg_stan_pred(X,
#                            sim_dat[sim_id == s & stage == 1, c(y1, y2)],
#                            sim_dat[sim_id == s & stage == 1, c(n1, n2)],
#                            X,
#                            sim_ppos[sim_id == s & stage == 1, c(ppos_m1, ppos_m2)],
#                            chains = 2, iter = 2000, refresh = 0)
#   Ptail <- mean(extract(mod1)$beta[, 2] < 0)
#   ypred <- data.table(extract(mod1)$ypred)[, .N, by = .(V1, V2)]
#   for(i in 1:nrow(ypred)) {
#     tmp <- logreg_stan(X,
#                        sim_dat[sim_id == s & stage == 1, c(y1, y2)] + ypred[i, c(V1, V2)],
#                        sim_dat[sim_id == s & stage == 1, c(n1, n2)] + sim_ppos[sim_id == s & stage == 1, c(ppos_m1, ppos_m2)],
#                        chains = 2, iter = 2000, refresh = 0)
#     ypred[i, P := mean(extract(tmp)$beta[, 2] < 0)]
#   }
#   c(ypred[, sum(N*(P > 0.9)) / sum(N)], Ptail)
# }
# proc.time() - pt
# 
# pp <- cbind(sim_ppos[stage == 1, .(ptail, ppos)], ptail_mcmc = Ppar[, 2], ppos_mcmc = Ppar[, 1])
# plot(pp$ptail, pp$ptail_mcmc); abline(0,1)
# plot(pp$ppos, pp$ppos_mcmc); abline(0,1)
# plot(pp$ptail - pp$ptail_mcmc)
# plot(pp$ppos - pp$ppos_mcmc)
