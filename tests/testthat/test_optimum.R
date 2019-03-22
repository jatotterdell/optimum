library(testthat)
library(optimum)
library(data.table)


context("post approx")




test_that("post probs approx similar", {
  
  simdat <- lapply(1, sim_trial_dat, p1tru = 0.2, p2tru = 0.1)
  aggdat <- lapply(simdat, agg_trial_dat, stage_n = c(25, 50, 75, 100))
  
  a1 = 1; b1 = 1; a2 = 1; b2 = 1
  
  d <- copy(aggdat[[1]])
  d1 <- copy(aggdat[[1]])
  d2 <- copy(aggdat[[1]])
  d3 <- copy(aggdat[[1]])
  
  calc_post <- beta_ineq_approx
  # Determine posterior parameters
  d1[, `:=`(a1 = a1 + y1,
           b1 = b1 + n1 - y1,
           a2 = a2 + y2,
           b2 = a2 + n2 - y2)]
  # Calculate posterior probabilities given current data and assuming follow-up of currently enrolled
  d1[, `:=`(post = calc_post(a1, b1, a2, b2),
           post_int = calc_post(a1 + w1, b1 + m1 - w1, a2 + w2, b2 + m2 - w2))]
  d1[, post_fin := post[.N]]
  d1
  
  calc_post <- beta_ineq
  # Determine posterior parameters
  d2[, `:=`(a1 = a1 + y1,
            b1 = b1 + n1 - y1,
            a2 = a2 + y2,
            b2 = a2 + n2 - y2)]
  # Calculate posterior probabilities given current data and assuming follow-up of currently enrolled
  d2[, `:=`(post = calc_post(a1, b1, a2, b2),
            post_int = calc_post(a1 + w1, b1 + m1 - w1, a2 + w2, b2 + m2 - w2))]
  d2[, post_fin := post[.N]]
  d2
  
  
  calc_post <- beta_ineq_sim
  # Determine posterior parameters
  d3[, `:=`(a1 = a1 + y1,
            b1 = b1 + n1 - y1,
            a2 = a2 + y2,
            b2 = a2 + n2 - y2)]
  # Calculate posterior probabilities given current data and assuming follow-up of currently enrolled
  d3[, `:=`(post = calc_post(a1, b1, a2, b2),
            post_int = calc_post(a1 + w1, b1 + m1 - w1, a2 + w2, b2 + m2 - w2))]
  d3[, post_fin := post[.N]]
  d3
  
  #
  expect_equal(d1$post, d2$post, tolerance = 0.01)
  expect_equal(d1$post, d3$post, tolerance = 0.01)
  expect_equal(d2$post, d3$post, tolerance = 0.01)
  
  #
  expect_equal(d1$post_int, d2$post_int, tolerance = 0.01)
  expect_equal(d1$post_int, d3$post_int, tolerance = 0.01)
  expect_equal(d2$post_int, d3$post_int, tolerance = 0.01)
  
  #
  expect_equal(d1$post_fin, d2$post_fin, tolerance = 0.01)
  expect_equal(d1$post_fin, d3$post_fin, tolerance = 0.01)
  expect_equal(d2$post_fin, d3$post_fin, tolerance = 0.01)
})


