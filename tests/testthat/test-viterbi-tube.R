test_that("Viterbi-tube paths are within epsilon", {
  fit <- example_hmmssf_list(n = 20, seed = 21)
  epsilon <- 3
  tube <- simulate_viterbi_tube_states(fit, epsilon = epsilon, n_sims = 10, seed = 22)

  expect_equal(ncol(tube$states), 20L)
  expect_true(all(tube$path_diagnostics$delta_from_viterbi <= epsilon + 1e-8))
  expect_true(all(tube$states >= 1 & tube$states <= 2))
})

test_that("Viterbi-tube paths work for a one-state model", {
  n <- 20
  fit <- list(
    initial = 1,
    transition = matrix(1, 1, 1),
    kernels = list(list(
      step = list(dist = "gamma", shape = 2, scale = 0.5),
      angle = list(dist = "wrapped_normal", mean = 0, sd = 1)
    )),
    observed_track = data.frame(x = seq_len(n), y = rep(0, n)),
    viterbi_path = rep(1L, n),
    posterior = matrix(1, nrow = n, ncol = 1),
    log_emission = matrix(0, nrow = n, ncol = 1)
  )

  tube <- simulate_viterbi_tube_states(fit, epsilon = 0, n_sims = 5, seed = 23)

  expect_equal(dim(tube$states), c(5L, n))
  expect_true(all(tube$states == 1L))
  expect_equal(dim(tube$occupancy), c(5L, 1L))
  expect_equal(tube$occupancy[, 1], rep(1, 5))
})
