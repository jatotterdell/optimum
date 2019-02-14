# library(mystan)
# library(rstan)
# library(data.table)
# 
# lpdf_pi <- function(p, mu, sigma) {
#   -(-mu + log(p / (1 - p)))^2 / (2*sigma^2) - 0.5*log(2*pi) - log(p - p^2) - log(sigma)
# }
# curve(exp(lpdf_pi(x, 0, 10)), xlim = c(0, 1), ylim = c(0, 2))
# curve(exp(lpdf_pi(x, 0, 1.75)), xlim = c(0, 1), ylim = c(0, 2))
# curve(dbeta(x, 1, 1), add = T, col = "red")
# 
# sims <- 10
# sim_dat <- sim_scenario(sims, p1tru = 0.23, p2tru = 0.2, n1int = c(500, 550), n2int = c(500, 550), post_method = "exact")
# sim_ppos <- calc_scenario_ppos(sim_dat, post_method = "exact")
# 
# P <- outer(0:50, 0:50, function(x, y) 
#   dbetabinom(x, 50, sim_dat[1, ]$a1, sim_dat[1, ]$b1) * dbetabinom(y, 50, sim_dat[1, ]$a2, sim_dat[1, ]$b2))
# Post <- outer(0:50, 0:50, function(x,y) 
#   Vectorize(beta_ineq)(sim_dat[1, ]$a1 + x, sim_dat[1, ]$b1 + 50 - x,
#                        sim_dat[1, ]$a2 + y, sim_dat[1, ]$b2 + 50 - y))
# 
# m <- logreg_stan_pred(cbind(1, c(0, 1)), c(sim_dat[1, y1], sim_dat[1, y2]), c(500, 500), X, c(50, 50), chains = 1, iter = 1000)
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
# X <- cbind(1, c(0, 1))
# PPoS <- rep(0, sims)
# Ptail <- rep(0, sims)
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
# Ppar <- foreach(s = 1:sims, .packages = c("mystan", "optimum", "rstan", "data.table"), .combine = c) %dopar% {
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
#   ypred[, sum(N*(P > 0.9)) / sum(N)]
# }
# stopImplicitCluster(cl)
# proc.time() - pt
# 
# pp <- cbind(sim_ppos[stage == 1, .(ptail, ppos)], Ptail, PPoS)
# plot(pp$ptail, pp$Ptail); abline(0,1)
# plot(pp$ppos, pp$PPoS); abline(0,1)
