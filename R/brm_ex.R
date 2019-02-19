# library(brms)
# library(simr)
# 
# trt <- c(0, 1)
# site <- 1:4
# hist <- c(0, 1)
# 
# var <- expand.grid(Intercept = 1, trt = trt, site = site, hist = hist, n = 1000)
# var$y <- rbinom(nrow(var), var$n, 0.5)
# 
# Xdes <- model.matrix( ~ trt*hist, data = var)
# 
# beta <- c(qlogis(0.1), -0.2, 0.1, 0.05)
# gama <- rnorm(4)
# linpred <- Xdes %*% beta + gama[var$site]
# p <- plogis(linpred)
# var$y <- rbinom(nrow(var), var$n, p)
# 
# mod <- glmer(cbind(y, n - y) ~ trt*hist + (1 | site), data = var, family = binomial)
# 
# priors <- c(prior(normal(0,0.75), class = "Intercept"),
#             prior(normal(0, 1), class = "b"))
# m <- brm(y | trials(n) ~ factor(trt)*factor(hist) + (1 | site), data = var, 
#          family = "binomial")
# 
# cbind(fixef(mod), fixef(m)[, 1])
# cbind(ranef(mod)$site, ranef(m)$site[, 1, 1])
# cbind(coef(mod)$site, coef(m)$site[,1,])
# 
