test_that("state occupancy proportions sum to one", {
  states <- c(1, 1, 2, 2, 2, 1)
  occ <- state_occupancy_vector(states, n_states = 2)
  expect_equal(sum(occ), 1)
})

test_that("residence diagnostics handle no-switch paths", {
  observed <- rep(1, 10)
  simulated <- matrix(rep(1, 30), nrow = 3)
  diag <- diagnostic_state_residence_time(observed, simulated)
  expect_true(any(grepl("low power", diag$warning, ignore.case = TRUE)))
})

test_that("transition-count diagnostics return finite rows", {
  observed <- c(1, 1, 2, 1, 2)
  simulated <- rbind(c(1, 2, 2, 1, 1), c(1, 1, 1, 2, 2))
  diag <- diagnostic_transition_counts(observed, simulated)
  expect_equal(diag$diagnostic, "transition_counts")
  expect_true(is.finite(diag$simulated_median))
})

test_that("Monte Carlo p-values treat numerical ties as ties", {
  observed <- 8.9e-16
  simulated <- rep(8.6e-16, 99)
  expect_equal(mc_rank_p_value(observed, simulated, alternative = "greater"), 1)
})

test_that("one-state occupancy and transition diagnostics are available", {
  observed <- rep(1L, 10)
  simulated <- matrix(1L, nrow = 9, ncol = 10)
  occupancy <- diagnostic_state_occupancy(observed, simulated)
  transition <- diagnostic_transition_counts(observed, simulated)

  expect_equal(occupancy$p_value, c(1, 1))
  expect_equal(transition$p_value, 1)
  expect_equal(occupancy$n_sims, c(9, 9))
  expect_equal(transition$n_sims, 9)
})
