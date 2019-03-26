library(testthat)
library(optimum)
library(data.table)


context("post approx")

test_that("post probs", {
  
  simdat <- lapply(1, sim_trial_dat, p1tru = 0.2, p2tru = 0.1)
  aggdat <- lapply(simdat, agg_trial_dat, stage_n = c(25, 50, 75, 100))
  prbdat1 <- lapply(aggdat, est_trial_prob, ppos_q = c(0.95, 0.975), post_method = "approx")
  
  aggdat <- lapply(simdat, agg_trial_dat, stage_n = c(25, 50, 75, 100))
  prbdat2 <- lapply(aggdat, est_trial_prob, ppos_q = c(0.95, 0.975), post_method = "exact")
  
  aggdat <- lapply(simdat, agg_trial_dat, stage_n = c(25, 50, 75, 100))
  prbdat3 <- lapply(aggdat, est_trial_prob, ppos_q = c(0.95, 0.975), post_method = "sim")

  expect_equal(prbdat1[[1]]$post, prbdat2[[1]]$post, tolerance = 0.05)
  expect_equal(prbdat1[[1]]$post_int, prbdat2[[1]]$post_int, tolerance = 0.05)
  expect_equal(prbdat1[[1]]$post_fin, prbdat2[[1]]$post_fin, tolerance = 0.05)
  
  expect_equal(prbdat3[[1]]$post, prbdat2[[1]]$post, tolerance = 0.05)
  expect_equal(prbdat3[[1]]$post_int, prbdat2[[1]]$post_int, tolerance = 0.05)
  expect_equal(prbdat3[[1]]$post_fin, prbdat2[[1]]$post_fin, tolerance = 0.05)
})
