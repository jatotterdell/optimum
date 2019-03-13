## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
	echo = FALSE,
	message = FALSE,
	warning = FALSE,
	collapse = TRUE,
	comment = "#>"
)

## ------------------------------------------------------------------------
library(knitr)
library(dplyr)
library(optimum)
library(bench)

## ---- fig.cap="Example Normal approximation to Beta densities.", fig.height=3, fig.width=4----
par(mfrow = c(2,2), oma = c(1,1,1,1), mar = c(1,1,1,1), cex = 0.5)
plot_beta_norm(2, 100, xlim = c(0, 0.2), xaxt = 'n', yaxt = 'n', main = "Beta(2,100)")
plot_beta_norm(5, 100, xlim = c(0, 0.2), xaxt = 'n', yaxt = 'n', main = "Beta(5,100)")
plot_beta_norm(10, 100, xlim = c(0, 0.2), xaxt = 'n', yaxt = 'n', main = "Beta(10,100)")
plot_beta_norm(15, 100, xlim = c(0, 0.2), xaxt = 'n', yaxt = 'n', main = "Beta(15,100)")

## ------------------------------------------------------------------------
a <- 5
b <- 100
c <- 10
d <- 100
m1 <- a / (a + b)
m2 <- c / (c + d)
v1 <- a*b / ( (a + b)^2 * (a + b + 1))
v2 <- c*d / ( (c + d)^2 * (c + d + 1))
y1 <- rbeta(1e5, a, b)
y2 <- rbeta(1e5, c, d)

par(oma = c(2,2,1,1), mar = c(1,1,1,1), cex = 0.75)
hist(y1 - y2, freq = F, breaks = 100, main = "X - Y; a = 5, b = 100, c = 10, d = 100")
x <- seq(min(y1- y2), max(y1 - y2), length.out = 1e3)
lines(x, dnorm(x, m1 - m2, sqrt(v1 + v2)))
legend("topleft", legend = "Normal approx", lty = 1, bty = 'n')

a <- 2
b <- 50
c <- 10
d <- 50
m1 <- a / (a + b)
m2 <- c / (c + d)
v1 <- a*b / ( (a + b)^2 * (a + b + 1))
v2 <- c*d / ( (c + d)^2 * (c + d + 1))
y1 <- rbeta(1e5, a, b)
y2 <- rbeta(1e5, c, d)

par(oma = c(2,2,1,1), mar = c(1,1,1,1), cex = 0.75)
hist(y1 - y2, freq = F, breaks = 100, main = "X - Y; a = 2, b = 50, c = 10, d = 50")
x <- seq(min(y1- y2), max(y1 - y2), length.out = 1e3)
lines(x, dnorm(x, m1 - m2, sqrt(v1 + v2)))
legend("topleft", legend = "Normal approx", lty = 1, bty = 'n')

## ------------------------------------------------------------------------
mark(
  beta_ineq(3, 100, 13, 90),
  beta_ineq_approx(3, 100, 13, 90),
  beta_ineq_sim(3, 100, 13, 90, sims = 1000),
  check = F,
  iterations = 1000
) %>%
  arrange(median) %>%
  select(expression, min, mean, median, max, `itr/sec`) %>%
  kable("html", row.names = F, booktabs = TRUE, digits = 2) %>%
  kableExtra::kable_styling(latex_options = "hold_position")

## ---- fig.cap="Deviation from exact value (adaptive quadrature) of $\\mathbb P(X>Y+\\delta)$."----
P_exact <- outer(0:50, 0:50, function(x, y) Vectorize(beta_ineq)(1+x, 1+50-x, 1+y, 1+50-y))
P_approx <- outer(0:50, 0:50, function(x, y) Vectorize(beta_ineq_approx)(1+x, 1+50-x, 1+y, 1+50-y))
P_sim <- outer(0:50, 0:50, function(x, y) Vectorize(beta_ineq_sim)(1+x, 1+50-x, 1+y, 1+50-y, sims = 1000))

par(mfrow = c(1, 2), mar = c(4,1,1,1), oma = c(0,1,1,1), mgp = c(2, 1, 0), cex = 0.7)
matplot(P_approx - P_exact, type = 'l', lty = 1, col = "grey50", ylim = c(-0.08, 0.08), main = "Approx", xlab = "a")
matplot(P_sim - P_exact, type = 'l', lty = 1, col = "grey50", ylim = c(-0.08, 0.08), main = "Sim (N = 1,000)", xlab = "a")

