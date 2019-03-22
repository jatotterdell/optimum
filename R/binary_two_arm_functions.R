# Define global variables to avoid
# notes from R CMD check
P <- NULL
N <- NULL
ptail <- NULL
sim_id <- NULL
a1 <- b1 <- a2 <- b2 <- NULL
parm <- stage <- variable <- value <- grp <- NULL


#' Calculate density of beta-binomial distribution
#' 
#' @useDynLib optimum, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' 
#' @param x The value at which to evaluate the density
#' @param n The sample size
#' @param a First parameter
#' @param b Second parameter
#' 
#' @return Value of beta-binomial(n,a,b) evaluated at x
#' 
#' @examples 
#' dbetabinom(5, 10, 2, 3)
#' 
#' @export
dbetabinom <- function(x, n, a = 1, b = 1){
  if(!(all(c(a, b) > 0))) stop("a and b must be > 0")
  if(any(n < 1)) stop("n must be > 0")
  if(any(x < 0)) stop("x must be >= 0")
  
  num <- lgamma(a + b) + lgamma(n + 1) + lgamma(x + a) + lgamma(n - x + b)
  den <- lgamma(a) + lgamma(b) + lgamma(x + 1) + lgamma(n - x + 1) + lgamma(n + a + b)
  prob <- exp(num - den)
  prob
}

#' Compare Normal approximation to Beta
#' 
#' @param a First Beta parameter
#' @param b Second Beta parameter
#' @param ... Other arguments to `curve()`
#' @return A plot of Beta density and Normal approximation
#' 
#' @export
plot_beta_norm <- function(a, b, ...) {
  curve(dbeta(x, a, b), ...)
  curve(dnorm(x, a/(a + b), sqrt( a*b / ((a+b)^2*(a+b+1)) )), add = TRUE, col = "red", ...)
}

#' Draw random variates from beta-binomial distribution
#' 
#' @import stats
#' 
#' @param n The number of random values to sample
#' @param m The sample size
#' @param a First parameter
#' @param b Second parameter
#' 
#' @examples
#' rbetabinom(2, 10, 2, 3)
#' 
#' @export
rbetabinom <- function(n, m, a = 1, b = 1) {
  if(!(all(c(a, b) > 0))) stop("a and b must be > 0")
  if(!(all(n > 0))) stop("n must be > 0")
  
  stats::rbinom(n, m, stats::rbeta(n, a, b))
}

#' Calculate the predicted probability of success
#' 
#' @import data.table
#' 
#' @param a First parameter of first beta random variable
#' @param b Second parameter of first beta random variable
#' @param c First paramter of second beta random variable
#' @param d Second parameter of second beta random variable
#' @param m1 Sample size to predict for first beta random variable
#' @param m2 Sample size to predict for second beta random variable
#' @param k_ppos The posterior probability cut-point to be assessed
#' 
#' @return The predicted probability of success
#' 
#' @export
calc_ppos <- function(a, b, c, d, m1, m2, k_ppos, post_method = "exact") {
  require(data.table)
  if(!(all(c(a, b, c, d, m1, m2) > 0))) stop("a, b, c, d, m1, m2 must be > 0")
  if(k_ppos < 0 | k_ppos > 1) stop("k_ppos must be in [0, 1]")
  
  calc_post <- switch(post_method,
                      "exact" = beta_ineq,
                      "approx" = beta_ineq_approx,
                      "sim" = beta_ineq_sim)
  
  y1pred <- rbetabinom(10000, m1, a, b)
  y2pred <- rbetabinom(10000, m2, c, d)
  ypred <- data.table(y1pred = y1pred, y2pred = y2pred)[, .N, keyby = list(y1pred, y2pred)]
  ypred[, `:=`(P = Vectorize(calc_post)(a + y1pred, 
                                    b + m1 - y1pred, 
                                    c + y2pred,
                                    d + m2 - y2pred))]
  ypred[, c(sum(N*(P > k_ppos)) / sum(N))]
}

beta_diff_dens <- function(x, a, b, c, d, ...) {
  require(appell)
  
  if(abs(x) >= 1) stop("x must be in [-1,1]")
  if(!(all(c(a, b, c, d) > 0))) stop("a, b, c, d must be > 0")
  
  if(x > 0) {
    return(exp(lbeta(c, b) - lbeta(a,b) - lbeta(c,d) + 
          (b + d - 1)*log(x) + (c + b - 1)*log(1 - x)) *
    as.numeric(appellf1(b, a + c + b + d - 2, 1 - a, b + c, 1 - x, 1 - x^2, ...)$val))
  } else if (x == 0) {
    return(exp(lbeta(a + c - 1, b + d - 1) - lbeta(a,b) - lbeta(c,d)))
  } else if ( x < 0) {
    return(exp(lbeta(a, d) - lbeta(a,b) - lbeta(c,d) +
        (b + d - 1)*log(-x) + (a + d - 1)*log(1 + x)) *
    as.numeric(appellf1(d, 1 - c, a + b + c + d - 2, a + d, 1 - x^2, 1 + x, ...)$val))
  }
}

# Estimate
# P(X > Y + delta)
# # where X ~ Beta(a, b), Y ~ Beta(c, d)
# beta_ineq <- function(a, b, c, d, delta = 0, ...) {
#   integrand <- function(x) { dbeta(x, a, b)*pbeta(x - delta, c, d) }
#   tryCatch(
#     integrate(integrand, delta, 1, ...)$value,
#     error = function(err) NA)
# }

#' Calculate Pr(X > Y + delta) where X and Y are independent Beta random variables
#' using numerical integration
#' 
#' @param a Parameter one of beta density for X
#' @param b Parameter two of beta density for X
#' @param c Parameter one of beta density for Y
#' @param d Parameter two of beta density for Y
#' @param delta The difference we wish to assess (i.e. X - Y > delta)
#' @param ... other arguments passed to integrate/quadgk function
#' 
#' @return The value of the integral
#' 
#' @examples 
#' beta_ineq(5, 5, 3, 7)
#' 
#' @export
beta_ineq <- function(a, b, c, d, delta = 0, ...) {
  
  if(!(all(c(a, b, c, d) > 0))) stop("a, b, c, d must be > 0")
  
  integrand <- function(x) { stats::dbeta(x, a, b)*stats::pbeta(x - delta, c, d) }
  tryCatch(
    pracma::quadgk(integrand, delta, 1, ...),
    error = function(err) NA)
}

#' Calculate Pr(X > Y + delta) where X and Y are independent Beta random variables
#' using Normal approximation.
#' 
#' @param a Parameter one of beta density for X
#' @param b Parameter two of beta density for X
#' @param c Parameter one of beta density for Y
#' @param d Parameter two of beta density for Y
#' @param delta The difference we wish to assess (i.e. X - Y > delta)
#' 
#' @return The value of the integral
#' 
#' @examples 
#' beta_ineq_approx(5, 5, 3, 7)
#' 
#' @export
beta_ineq_approx <- function(a, b, c, d, delta = 0) {
  if(!(all(c(a, b, c, d) > 0))) stop("a, b, c, d must be > 0")
  
  m1 <- a / (a + b)
  v1 <- a*b / ( (a + b)^2 * (a + b + 1))
  m2 <- c / (c + d)
  v2 <- c*d / ( (c + d)^2 * (c + d + 1))
  z <- (m1 - m2 - delta) / sqrt(v1 + v2)
  return(stats::pnorm(z))
}

#' Calculate Pr(X > Y + delta) where X and Y are independent Beta random variables
#' using Monte Carlo method.
#' 
#' @param a Parameter one of beta density for X
#' @param b Parameter two of beta density for X
#' @param c Parameter one of beta density for Y
#' @param d Parameter two of beta density for Y
#' @param delta The difference we wish to assess (i.e. X - Y > delta)
#' @param sims The number of Monte Carlo variates to generate for estimation
#' 
#' @return The value of the integral
#' 
#' @examples 
#' beta_ineq_sim(5, 5, 3, 7)
#' 
#' @export
beta_ineq_sim <- function(a, b, c, d, delta = 0, sims = 10000) {
  if(!(all(c(a, b, c, d) > 0))) stop("a, b, c, d must be > 0")
  
  lens <- unlist(lapply(list(a, b, c, d), length))
  if(any(max(lens)- min(lens) != 0)) stop("a, b, c, d must be same len")
  
  X <- lapply(1:length(a), function(x) stats::rbeta(sims, a[x], b[x]))
  Y <- lapply(1:length(a), function(x) stats::rbeta(sims, c[x], d[x]))
  
  means <- lapply(1:length(a), function(x) mean(X[[x]] > Y[[x]] + delta))
  unlist(means)
}

#' Simulate trial data using Poisson process for accrual rate
#'
#' @param sim_id The id for the simulated trial data
#' @param p1tru True response rate under control
#' @param p2tru True response rate under treatment
#' @param nmax The maximum total sample size
#' @param enro_rate The baseline acrrual rate
#' @param enro_intensity Function for changing acrrual intensity
#' @param resp_delay Function returning response time value
#' @export
sim_trial_dat <- function(
  sim_id = 1,
  p1tru  = 0.1,
  p2tru  = 0.1,
  nmax = 100,
  enro_rate = 1,
  enro_intensity = function(t) 1,
  resp_delay = function(n) runif(n)
) {
  library(poisson)
  require(data.table)
  require(randomizr)
  
  if(!(all(c(p1tru, p2tru) > 0) & all(c(p1tru, p2tru) < 1))) 
    stop("p1tru and p2tru must be in (0,1)")
  
  enro_t <- nhpp.sim(rate = enro_rate, num.events = nmax, enro_intensity, prepend.t0 = F)
  resp_t <- enro_t + resp_delay(nmax)
  
  d <- data.table(p1tru = p1tru,
                  p2tru = p2tru,
                  pid = 1:nmax,
                  enro_t = enro_t,
                  resp_t = resp_t,
                  x = as.numeric(complete_ra(nmax, num_arms = 2)))
  
  d[order(resp_t), resp_o := 1:.N]
  d[, y := rbinom(nmax, 1, prob = ifelse(x == 1, p1tru, p2tru))]
  return(d)
}

#' Aggregate trial data over treatment
#' 
#' @param d The raw trial data.
#' @param stage_n The number of responses required to trigger each interim.
#' @param min_rem The minimum number of subjects still to be enrolled for an interim to trigger.
#' @return An aggregated trial dataset.
#' @export
agg_trial_dat <- function(d, stage_n, min_rem = 10) {
  resp_cut <- sort(d[, resp_t])[stage_n]
  
  # for all cutpoints, by group
  # n participants with resp_time <= cut point
  # y participants with events with resp_time <= cut point
  # m participants enrolled but without resp_time available at cut point
  # w participants with events that are enrolled but without resp_time available at cut point 
  # l participants with enrolment occurring after cut point
  # z participants with events with enrolment occurring after cut point
  
  dd <- d[,
          {
            list(resp = stage_n,
                 n = sapply(resp_cut, function(a) sum(resp_t <= a)),
                 y = sapply(resp_cut, function(a) sum((resp_t <= a)*y)),
                 m = sapply(resp_cut, function(a) sum(resp_t > a & enro_t <= a)),
                 w = sapply(resp_cut, function(a) sum((resp_t > a & enro_t <= a)*y)),
                 l = sapply(resp_cut, function(a) sum(enro_t > a)),
                 z = sapply(resp_cut, function(a) sum((enro_t > a)*y))
            )
          }, by = .(p1tru, p2tru, x)]
  
  # spread to wide.
  dcast(dd, resp + p1tru + p2tru ~ x, 
        value.var = c("n", "y", "m", "w", "l", "z"), sep = "")[(l1 + l2) > min_rem | resp == max(stage_n)]
}

#' Calculate trial probabilities (posterior and predictive)
#' 
#' @param d The Trial data
#' @param ppos_q The vector of cut-points to consider
#' @param ppos_sim Number of simulations from posterior predictive to use in PPoS(q) calculation
#' @param post_method Method used to calculate P(X > Y + delta), one of `exact` (integration), `approx` (normal), `sim` (monte carlo)
#' @param a1 Prior shape for treatment 1
#' @param b1 Prior scale for treatment 1
#' @param a2 Prior shape for treatment 2
#' @param b2 Prior scale for treatment 2
#' @return Updates `d` inplace but also returns the updated `d`. 
#' `post` - current posterior probability, 
#' `post_int` - posterior probability when follow-up enrolled individuals,
#' `post_fin` - posterior probability when follow-up to maximum sample size,
#' `ppos_intq` - predictive probability of success given q if follow-up enrolled individuals,
#' `ppos_finq` - predictive probability of success given q if follow-up to maximum sample size,
#' 
#' @export
est_trial_prob <- function(
  d, 
  ppos_q = 0.95, 
  ppos_sim = 1000, 
  post_method = "approx",
  a1 = 1, b1 = 1, a2 = 1, b2 = 1
) {
  if(!(all(c(a1, b1, a2, b2) > 0))) 
    stop("a1, b1, a2, b2 must be > 0")
  if(!post_method %in% c("exact", "approx", "sim")) 
    stop("post_method must be 'exact', 'approx', or 'sim'.")
  
  # What method used to estimate posterior probability?
  calc_post <- switch(post_method,
                      "exact" = beta_ineq,
                      "approx" = beta_ineq_approx,
                      "sim" = beta_ineq_sim)
  # Setup storage for prediction
  ypred <- data.table(y1 = rep(0, ppos_sim),
                      y2 = rep(0, ppos_sim))
  # Determine posterior parameters
  d[, `:=`(a1 = a1 + y1,
           b1 = b1 + n1 - y1,
           a2 = a2 + y2,
           b2 = a2 + n2 - y2)]
  # Calculate posterior probabilities given current data and assuming follow-up of currently enrolled
  d[, `:=`(post = calc_post(a1[1:.N], b1[1:.N], a2[1:.N], b2[1:.N]),
           post_int = calc_post(a1 + w1, b1 + m1 - w1, a2 + w2, b2 + m2 - w2))]
  d[, post_fin := post[.N]]
  
  # Calculate PPoS at each interim for currently enrolled and complete enrolment
  if(!is.null(ppos_q)) {
    for(i in 1:(nrow(d) - 1)) {
      # Do interim PPoS calculation
      ypred[, `:=`(y1 = rbetabinom(ppos_sim, d[i, m1], d[i, a1], d[i, b1]),
                   y2 = rbetabinom(ppos_sim, d[i, m2], d[i, a2], d[i, b2]))]
      ypred_agg <- ypred[, .N, by = .(y1, y2)]
      ypred_agg[, P := Vectorize(calc_post)(d[i, a1] + y1, 
                                            d[i, b1 + m1] - y1, 
                                            d[i, a2] + y2,
                                            d[i, b2 + m2] - y2)]
      d[i, paste0('ppos_int', ppos_q) := ypred_agg[, lapply(ppos_q, function(q) sum(N*(P > q)) / sum(N))]]
      
      # Do final PPoS calculation
      ypred[, `:=`(y1 = rbetabinom(ppos_sim, d[i, m1 + l1], d[i, a1], d[i, b1]),
                   y2 = rbetabinom(ppos_sim, d[i, m2 + l2], d[i, a2], d[i, b2]))]
      ypred_agg <- ypred[, .N, by = .(y1, y2)]
      ypred_agg[, P := Vectorize(calc_post)(d[i, a1] + y1, 
                                            d[i, b1 + m1 + l1] - y1, 
                                            d[i, a2] + y2,
                                            d[i, b2 + m2 + l2] - y2)]
      d[i, paste0('ppos_fin', ppos_q) := ypred_agg[, lapply(ppos_q, function(q) sum(N*(P > q)) / sum(N))]]
    } 
  }
  return(d)
}

#' Apply decision rule to trial data
#' 
#' @param trial Aggregated trial data with estimated probabilities (after call to est_trial_prob).
#' @param fut_k Futility boundary vector. Must has as many elements as the maximum number of interims in `trial`.
#' @param suc_k Success boundary vector.
#' @param inf_k Inferiority cut-off.
#' @param sup_k Superiority boundary vector.
#' @return A new `data.table` given the decision made and boundaries used
#' @export
dec_trial <- function(
  trial,
  fut_k = 0.1,
  suc_k = 0.9,
  inf_k = 0.05,
  sup_k = 0.95, ...) {
  
  # Use the ppos column which matches sup_k
  ppos_cols <- grep(paste0(sup_k, "$"), names(trial), value = T)
  if(length(ppos_cols) == 0) stop("sup_k did not match any columns in trial.")
  max_stage <- trial[, .N, by = sim_id][, max(N)]
  
  if (length(fut_k) == 1)
    fut_k <- rep(fut_k, max_stage - 1)
  if (length(suc_k) == 1)
    suc_k <- rep(suc_k, max_stage - 1)
  if (length(fut_k) < (max_stage - 1))
    stop("fut_k has too few elements")
  if (length(suc_k) < (max_stage - 1))
    stop("suc_k has too few elements")
  

  trial[,
        {
          fut <- match(TRUE, get(ppos_cols[2])[-.N] < fut_k[1:(.N-1)])
          suc <- match(TRUE, get(ppos_cols[1])[-.N] > suc_k[1:(.N-1)])
          if(any(!is.na(c(fut, suc)))) {
            res <- which.min(c(fut, suc))
            int <- min(c(fut, suc), na.rm = TRUE)
            list(p1tru = p1tru[1],
                 p2tru = p2tru[1],
                 res = switch(res, "futile", "expect success"),
                 fin = ifelse(post_int[int] > sup_k, "success", "failure"),
                 stage = int,
                 resp = n1[int] + n2[int],
                 enro = n1[int] + n2[int] + m1[int] + m2[int],
                 post = post[int],
                 post_int = post_int[int],
                 post_fin = post_fin[int],
                 ppos_int = get(ppos_cols[1])[int],
                 ppos_fin = get(ppos_cols[2])[int],
                 fut_k = paste(unique(fut_k), collapse = ","),
                 suc_k = paste(unique(suc_k), collapse = ","),
                 inf_k = inf_k,
                 sup_k = sup_k)
          } else {
            inf <- tail(post, 1) < tail(inf_k, 1)
            sup <- tail(post, 1) > tail(sup_k, 1)
            list(p1tru = p1tru[1],
                 p2tru = p2tru[1],
                 res = ifelse(inf, "inferior", ifelse(sup, "superior", "inconclusive")),
                 fin = ifelse(post[.N] > sup_k, "success", "failure"),
                 stage = .N,
                 resp = n1[.N] + n2[.N],
                 enro = n1[.N] + n2[.N] + m1[.N] + m2[.N],
                 post = post[.N],
                 post_int = post_int[.N],
                 post_fin = post_fin[.N],
                 ppos_int = get(ppos_cols[1])[.N],
                 ppos_fin = get(ppos_cols[2])[.N],
                 fut_k = paste(unique(fut_k), collapse = ","),
                 suc_k = paste(unique(suc_k), collapse = ","),
                 inf_k = inf_k,
                 sup_k = sup_k)
          }
        }, by = sim_id]
}


dec_trial_post <- function(
  trial,
  inf_k = 0.05,
  sup_k = 0.95  
) {
  
  max_stage <- trial[, .N, by = sim_id][, max(N)]
  
  if (length(inf_k) == 1)
    inf_k <- rep(inf_k, max_stage)
  if (length(sup_k) == 1)
    sup_k <- rep(sup_k, max_stage)
  if (length(inf_k) < max_stage)
    stop("inf_k has too few elements")
  if (length(sup_k) < max_stage)
    stop("sup_k has too few elements")
  
  trial[,
        {
          inf <- match(TRUE, post[1:(.N-1)] < inf_k[1:(.N-1)])
          sup <- match(TRUE, post[1:(.N-1)] > sup_k[1:(.N-1)])
          if(any(!is.na(c(inf, sup)))) {
            res <- which.min(c(inf, sup))
            int <- min(c(inf, sup), na.rm = TRUE)
            list(p1tru = p1tru[1],
                 p2tru = p2tru[1],
                 res = switch(res, "early failure", "early success"),
                 fin = ifelse(post_fin[.N] > sup_k[.N], "success", "failure"),
                 stage = int,
                 n1 = n1[int],
                 n2 = n2[int],
                 y1 = y1[int],
                 y2 = y2[int],
                 enro = n1[int] + n2[int] + m1[int] + m2[int],
                 post = post[int],
                 post_int = post_int[int],
                 post_fin = post_fin[int],
                 inf_k = paste(unique(inf_k), collapse = ","), 
                 sup_k = paste(unique(sup_k), collapse = ","))
          } else {
            inf <- tail(post, 1) < tail(inf_k, 1)
            sup <- tail(post, 1) > tail(sup_k, 1)
            list(p1tru = p1tru[1],
                 p2tru = p2tru[1],
                 res = ifelse(inf, "inferior", ifelse(sup, "superior", "inconclusive")),
                 fin = ifelse(post[.N] > sup_k[.N], "success", "failure"),
                 stage = .N,
                 n1 = n1[.N],
                 n2 = n2[.N],
                 y1 = y1[.N],
                 y2 = y2[.N],
                 enro = n1[.N] + n2[.N] + m1[.N] + m2[.N],
                 post = post[.N],
                 post_int = post_int[.N],
                 post_fin = post_fin[.N],
                 inf_k = paste(unique(inf_k), collapse = ","), 
                 sup_k = paste(unique(sup_k), collapse = ","))
          }
        }, by = sim_id]
}


#=======================#
# OLD FUNCTIONS BELOW...#
#=======================#



#' Simulate Bayesian two-arm trial 
#' with early termination for futility/success
#' Model:
#'  p_1      ~ beta(a, b)
#'  p_2      ~ beta(c, d)
#'  x_i|p_i  ~ binomial(n_i, p_i), k = 1...K, i = 1,2
#'  p_i|x_i  ~ beta(a + sum(x_i), b + n_i - sum(x_i))
#'
#' Hypothesis:
#'   H_0: p_1 > p_2 + d
#'   H_1: p_1 < p_2 + d
#'
#' Terimnal decision rule:
#'   Pr(p_1 < p_2 + d | x_i) > k_hi => H_1
#'   Pr(p_1 < p_2 + d | x_i) < k_lo => H_0
#'   Otherwise => inconclusive
#'
#' Interim decision rule:
#'   Pr(p_1 < p_2 + d | x_i) > k_hi => H_1
#'   Pr(p_1 < p_2 + d | x_i) < k_lo => H_0
#'
#'            OR
#'
#'   PPoS > q_hi => expect success => H_1
#'   PPoS < q_lo => futile => H_0
#'
#'   Otherwise => continue to k + 1
#'   
#' @param sim_id The simulation ID reference number
#' @param p1tru True value for response rate arm 1
#' @param p2tru True value for response rate arm 2
#' @param delta Relative difference of interest (X - Y > delta)
#' @param n1int Individuals with follow-up at each interim arm 1
#' @param n2int Individuals with follow-up at each interim arm 2
#' @param a1 Prior parameter arm 1
#' @param b1 Prior parameter arm 1
#' @param a2 Prior parameter arm 2
#' @param b2 Prior parameter arm 2
#' @param post_method What method to use for estimating the Beta inequality. One of `exact`, `approx` or `sim`.
#' 
#' @return A data.table of the simulated trial containing 
#' one row for each interim analysis.
#' 
#' @import data.table
#' 
#' @export
sim_trial <- function(
  sim_id = 1,
  p1tru  = 0.1,
  p2tru  = 0.1,
  delta  = 0.0,
  n1int  = seq(100, 500, 50),
  n2int  = seq(100, 500, 50),
  a1      = 1,
  b1      = 1,
  a2      = 1,
  b2      = 1,
  post_method = "exact"
) {
  
  if(!post_method %in% c("exact", "approx", "sim")) 
    stop("post_method must be either 'exact', 'approx', or 'sim'.")
  if(length(n1int) != length(n2int)) 
    stop("n1int and n2int must be of same dimension.")
  calc_post <- switch(post_method,
                      "exact" = beta_ineq,
                      "approx" = beta_ineq_approx,
                      "sim" = beta_ineq_sim)
  
  n_analyses <- length(n1int)
  
  # Generate data
  x1 <- rbinom(length(n1int), diff(c(0, n1int)), p1tru)
  x2 <- rbinom(length(n2int), diff(c(0, n2int)), p2tru)
  y1 <- cumsum(x1)
  y2 <- cumsum(x2)
  
  # Simulate interim analyses
  sim_interim <- function(interim) {
    # Posterior parameters
    y1curr <- y1[interim]
    y2curr <- y2[interim]
    n1curr <- n1int[interim]
    n2curr <- n2int[interim]
    a1post <- a1 + y1curr
    b1post <- b1 + n1curr - y1curr
    a2post <- a2 + y2curr
    b2post <- b2 + n2curr - y2curr
    
    # Posterior probability of p_1 > p_2
    ptail <- calc_post(a1post, b1post, a2post, b2post, delta)
    
    return(c(n1 = n1curr, n2 = n2curr, 
             x1 = x1[interim], x2 = x2[interim],
             y1 = y1curr, y2 = y2curr,
             a1 = a1post, b1 = b1post, a2 = a2post, b2 = b2post, 
             ptail = ptail))
  }
  # Iterate through all the interims
  intr <- lapply(1:n_analyses, sim_interim)
  # Collect the results
  sim_trial <- cbind(sim_id = sim_id, stage = 1:n_analyses,
                          p1tru = p1tru, p2tru = p2tru, delta = delta,
                          as.data.frame(do.call(rbind, intr)))
  return(data.table::as.data.table(sim_trial))
}

#' Simulate a set of Bayesian two-arm trials
#' with early termination for futility/success
#' for a given scenario.
#' 
#' @param sims The number of simulations to undertake for the scenario
#' @param ... Other arguments to `sim_trial`` function
#' 
#' @return A data.table of multiple simulated trials using the same parameters
#' 
#' @export
sim_scenario <- function(sims, ...) {
  res <- lapply(1:sims, function(i, ...) sim_trial(i, ...), ...)
  return(rbindlist(res))
}

calc_trial_ppos <- function(
  trial, 
  m1 = NULL, # Enrolled with no follow-up
  m2 = NULL,
  l1 = NULL, # Maximum left to follow-up
  l2 = NULL,
  k_ppos = 0.9,
  ppos_method = "sim", 
  post_method = "approx",
  pp_sim = 1000,
  ppos_name = "ppos") {
  
  if(! post_method %in% c("exact", "approx")) stop("post_method must be either 'exact' or 'approx'.")
  if(! ppos_method %in% c("sim", "exact")) stop("ppos_method must be either 'sim' or 'exact'.")
  if(any(is.null(c(m1, m2, l1, l2)))) {
    stop("Must specify sample size numbers m1, m2, l1 and l2.")
  }
  if(length(m1) != length(m2)) stop("m1 and m2 must have same dimension.")
  if(length(l1) != length(l2)) stop("l1 and l2 must have same dimension.")
  if(length(trial$n1) != length(m1)) stop("n1 and m1 must have same dimension.")
  if(length(trial$n1) != length(l1)) stop("n1 and l1 must have same dimension.")
  
  calc_post <- ifelse(post_method == "exact", beta_ineq, beta_ineq_approx)
  
  n_analyses <- nrow(trial)
  n_max <- sum(trial[n_analyses, c(n1, n2)])
  ppos_int <- rep(NA, n_analyses)
  ppos_fin <- rep(NA, n_analyses)
  
  for(i in 1:n_analyses) {
    # If no missing data, or number enrolled equal to maximum sample size, no point
    # in doing PPoS check
    current_enrolled <- trial$n1[i] + trial$n2[i] + m1[i] + m2[i]

    # If current enrolled is more than maximum sample size,
    # no point in stopping enrollment...
    if(current_enrolled >= n_max) {
      ppos_fin[i] <- NA
      ppos_int[i] <- NA
    } else {
      # If no missing data to predict, just check ptail
      if ((l1[i] == 0 & l2[i] == 0)) {
        ppos_fin[i] <- as.numeric(trial$ptail[i] > k_ppos)
      } else {
        if (ppos_method == "sim") {
          y1pred <- rbetabinom(pp_sim, l1[i], trial$a1[i], trial$b1[i])
          y2pred <- rbetabinom(pp_sim, l2[i], trial$a2[i], trial$b2[i])
          
          # No point computing posterior for duplicate values
          # just do once and multiply by the frequency
          ypred <- data.table(y1pred, y2pred)[, .N, by = list(y1pred, y2pred)]
          ypred[, P := Vectorize(calc_post)(trial$a1[i] + y1pred, 
                                            trial$b1[i] + l1[i] - y1pred, 
                                            trial$a2[i] + y2pred,
                                            trial$b2[i] + l2[i] - y2pred)]
          ppos_fin[i] <- ypred[, c(sum(N*(P > k_ppos)) / sum(N))]
        }
      }
      if ((m1[i] == 0 & m2[i] == 0)) {
        ppos_int[i] <-  as.numeric(trial$ptail[i] > k_ppos)
      } else {
        if (ppos_method == "sim") {
          y1pred <- rbetabinom(pp_sim, m1[i], trial$a1[i], trial$b1[i])
          y2pred <- rbetabinom(pp_sim, m2[i], trial$a2[i], trial$b2[i])
          
          # No point computing posterior for duplicate values
          # just do once and multiply by the frequency
          ypred <- data.table(y1pred, y2pred)[, .N, by = list(y1pred, y2pred)]
          ypred[, P := Vectorize(calc_post)(trial$a1[i] + y1pred, 
                                            trial$b1[i] + m1[i] - y1pred, 
                                            trial$a2[i] + y2pred,
                                            trial$b2[i] + m2[i] - y2pred)]
          ppos_int[i] <- ypred[, c(sum(N*(P > k_ppos)) / sum(N))]
        }
      }
    }
  }
  trial[, `:=`(m1 = m1, m2 = m2, l1 = l1, l2 = l2)]
  trial[, `:=`(ppos_int = ppos_int, ppos_fin = ppos_fin, ppos_cut = k_ppos)]
  invisible(trial)
}

#' Calculate the Predictive Probability of Success for
#' a given trial scenario. Repeatedly calls calc_trial_ppos.
#' 
#' @param scenario A data.table returned from `sim_scenario`.
#' @param useParallel Use parallel processing for calculating PPoS.
#' @param ... Other arguements to `calc_trial_ppos`.
#' 
#' @return A data.table of scenario with the PPoS results incorporated.
#' 
#' @export
calc_scenario_ppos <- function(scenario, useParallel = FALSE, ...) {
  scenario_split <- split(scenario, scenario$sim_id)
  if(useParallel) {
    cl <- parallel::makeCluster(parallel::detectCores())
    parallel::clusterExport(cl, c("beta_ineq", "beta_ineq_approx", "rbetabinom"))
    res <- parallel::parLapply(cl, scenario_split, calc_trial_ppos, ...)
    parallel::stopCluster(cl)
     
  } else {
    res <- lapply(scenario_split, calc_trial_ppos, ...)
  }
  return(data.table::rbindlist(res))
}

#' Apply a decision rule to a trial from `sim_trial`
#' 
#' @import utils
#' 
#' @param trial The data.table fro the trial in question
#' @param fut_var The name of the variable used for determining futility at interim analyses.
#' @param suc_var The name of the variable used for determining success at interim analyses.
#' @param fut_k The probability cut-off for futility.
#' @param suc_k The probability cut-off for success
#' @param inf_k The probability cut-off at final analysis for declaring theta_1 < theta_2
#' @param sup_k The probability cut-off at final analysis for declaring theta_1 > theta_2
#' 
#' @return A data.table of the stage at which a decision was made.
#' 
#' @export
decide_trial <- function(
  trial,
  fut_k = 0.1,
  suc_k = 0.9,
  inf_k = 0.05,
  sup_k = 0.95) {
  
  trial[,
    {
      fut <- match(TRUE, ppos_fin[-.N] < fut_k)
      suc <- match(TRUE, ppos_int[-.N] > suc_k)
      if(any(!is.na(c(fut, suc)))) {
        res <- which.min(c(fut, suc))
        int <- min(c(fut, suc), na.rm = TRUE)
        list(res = switch(res, "futile", "expect success"),
             stage = int,
             fut_k = paste(fut_k, collapse = ","),
             suc_k = paste(suc_k, collapse = ","),
             inf_k = inf_k,
             sup_k = sup_k)
      } else {
        inf <- tail(ptail, 1) < tail(inf_k, 1)
        sup <- tail(ptail, 1) > tail(sup_k, 1)
        list(res = ifelse(inf, "inferior", ifelse(sup, "superior", "inconclusive")),
             stage = .N,
             fut_k = paste(fut_k, collapse = ","),
             suc_k = paste(suc_k, collapse = ","),
             inf_k = inf_k,
             sup_k = sup_k)
      }
    }, by = sim_id]
}

#' Plot a trial
#' 
#' @param trial A trial output from sim_trial
#' @param par_seq The sequence of parameter values at which to evaluate the densities
#' 
#' @export
plot_trial <- function(trial, par_seq) {
  if(requireNamespace('ggridges', quietly = TRUE)) {
    pal <- viridisLite::viridis(2)
    d1 <- data.table::melt(trial[, lapply(par_seq, function(x) dbeta(x, a1, b1)), by = list(stage, s1 = a1, s2 = b1)], 
               id.vars = c("stage", "s1", "s2"))[order(stage, variable), `:=`(grp = "1", parm = par_seq)]
    d2 <- data.table::melt(trial[, lapply(par_seq, function(x) dbeta(x, a2, b2)), by = list(stage, s1 = a2, s2 = b2)], 
               id.vars = c("stage", "s1", "s2"))[order(stage, variable), `:=`(grp = "2", parm = par_seq)]
    d <- data.table::rbindlist(list(d1, d2))
    ggplot2::ggplot(d,
           ggplot2::aes(x = parm, y = stage, height = value)) +
      ggridges::geom_density_ridges(ggplot2::aes(fill = paste(stage, grp)), stat = "identity",
                          alpha = 0.8, colour = "grey50") +
      ggridges::scale_fill_cyclical(
        labels = c(bquote(theta[a]~"="~.(unique(trial$p1tru))), 
                   bquote(theta[w]~"="~.(unique(trial$p2tru)))),
        # values = c("#ff8080", "#8080ff"),
        values = c(pal[1], pal[2]),
        name = "Parameter", guide = "legend"
      ) +
      ggplot2::labs(x = bquote(theta), y = "Stage")  
  }
}


# 
# set.seed(123)
# trial <- sim_trial(p1tru = 0.105,
#                    n1int = seq(500, 3000, 500),
#                    n2int = seq(500, 3000, 500))
# trial_ppos <- calc_trial_ppos(trial, post_method = "approx", pp_sim = 1e5)
# trial_ppos <- calc_trial_ppos(trial, m1int = c(rep(50, 5), 0), m2int = c(rep(50, 5), 0),
#                               post_method = "approx", pp_sim = 1e5)
# 
# PP <- outer(0:500, 0:500, function(x,y) Vectorize(beta_ineq)(244+x,2258+500-x,223+y,2279+500-y))
# pp <- outer(0:500, 0:500, function(x,y) dbetabinom(x, 500, 244, 2258)*dbetabinom(y, 500, 223, 2279))
# sum((PP > 0.9)*pp)
# 
# PP <- outer(0:50, 0:50, function(x,y) Vectorize(beta_ineq)(244+x,2258+50-x,223+y,2279+50-y))
# pp <- outer(0:50, 0:50, function(x,y) dbetabinom(x, 50, 244, 2258)*dbetabinom(y, 50, 223, 2279))
# sum((PP > 0.9)*pp)
# 
# 
# sce <- sim_scenario(100)
# sce <- calc_scenario_ppos(sce, k_ppos = 0.95)
# sce <- calc_scenario_ppos(sce, m1int = c(rep(10, 8), 0), m2int = c(rep(10, 8), 0), k_ppos = 0.95)
