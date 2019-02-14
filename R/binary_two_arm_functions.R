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
calc_ppos <- function(a, b, c, d, m1, m2, k_ppos) {
  y1pred <- rbetabinom(10000, m1, a, b)
  y2pred <- rbetabinom(10000, m2, c, d)
  ypred <- data.table(y1pred = y1pred, y2pred = y2pred)[, data.table::.N, keyby = list(y1pred, y2pred)]
  ypred[, `:=`(P = Vectorize(beta_ineq)(a + y1pred, 
                                    b + m1 - y1pred, 
                                    c + y2pred,
                                    d + m2 - y2pred))]
  ypred[, c(sum(N*(P > k_ppos)) / sum(N))]
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
  X <- stats::rbeta(sims, a, b)
  Y <- stats::rbeta(sims, c, d)
  mean(X > Y + delta)
}

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
  m1int = NULL,
  m2int = NULL,
  k_ppos = 0.9,
  ppos_method = "sim", 
  post_method = "approx",
  pp_sim = 1000,
  ppos_name = "ppos") {
  
  if(! post_method %in% c("exact", "approx")) stop("post_method must be either 'exact' or 'approx'.")
  if(! ppos_method %in% c("sim", "exact")) stop("ppos_method must be either 'sim' or 'exact'.")
  if(is.null(m1int) || is.null(m2int)) {
    m1int = max(trial$n1) - trial$n1
    m2int = max(trial$n2) - trial$n2
  }
  if(length(m1int) != length(m2int)) stop("m1int and m2int must have same dimension.")
  if(length(trial$n1) != length(m1int)) stop("n1 and m1int must have same dimension.")
  
  calc_post <- ifelse(post_method == "exact", beta_ineq, beta_ineq_approx)
  
  n_analyses <- nrow(trial)
  ppos <- rep(0, n_analyses)
  
  for(i in 1:n_analyses) {
    if ((m1int > 0) || (m2int > 0)) {
      
      if (ppos_method == "sim") {
        y1pred <- rbetabinom(pp_sim, m1int[i], trial$a1[i], trial$b1[i])
        y2pred <- rbetabinom(pp_sim, m2int[i], trial$a2[i], trial$b2[i])
        
        # No point computing posterior for duplicate values
        # just do once and multiply by the frequency
        ypred <- data.table(y1pred, y2pred)[, .N, by = list(y1pred, y2pred)]
        ypred[, P := Vectorize(calc_post)(trial$a1[i] + y1pred, 
                                          trial$b1[i] + m1int[i] - y1pred, 
                                          trial$a2[i] + y2pred,
                                          trial$b2[i] + m2int[i] - y2pred)]
        ppos[i] <- ypred[, c(sum(N*(P > k_ppos)) / sum(N))]
      }
      
      if (ppos_method == "exact") {
        stop("Exact not implemented.")
        P <- outer(0:m1int, 0:m2int, function(x, y) 
          dbetabinom(x, m1int, trial$a1[i], trial$b1[i]) * dbetabinom(y, m2int, trial$a2[i], trial$b2[i]))
        Post <- outer(0:m1int, 0:m1int, function(x,y) 
          Vectorize(calc_post)(trial$a1[i] + x, trial$b1[i] + m1int - x,
                               trial$a2[i] + y, trial$b2[i] + m2int - y))
      }
      
    }
  }
  trial[, (ppos_name) := ppos]
  trial[, (paste0(ppos_name, "_cut")) := k_ppos]
  trial[, (paste0(ppos_name, "_m1")) := m1int]
  trial[, (paste0(ppos_name, "_m2")) := m2int]
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
  fut_var = "ppos_final",
  suc_var = "ppos_interim",
  fut_k = 0.1,
  suc_k = 0.9,
  inf_k = 0.05,
  sup_k = 0.95) {
  
  trial[,
    {
      fut <- match(TRUE, get(fut_var)[-.N] < fut_k)
      suc <- match(TRUE, get(suc_var)[-.N] > suc_k)
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
