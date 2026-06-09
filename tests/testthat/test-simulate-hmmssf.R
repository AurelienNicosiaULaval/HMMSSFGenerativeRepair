test_that("simulate_generative_hmmssf returns valid simulated tracks", {
  fit <- example_hmmssf_list(n = 30, seed = 11)
  sim <- simulate_generative_hmmssf(
    fit = fit,
    observed_track = fit$observed_track,
    n_sims = 5,
    method = "markov",
    parameter_uncertainty = FALSE,
    seed = 12
  )

  expect_s3_class(sim, "gmov_hmmssf_simulations")
  expect_length(sim$simulated_tracks, 5)
  expect_equal(dim(sim$simulated_states), c(5L, 30L))
  expect_true(all(vapply(sim$simulated_tracks, nrow, integer(1L)) == 30L))
})

test_that("Viterbi-fixed simulation preserves the exact state sequence", {
  fit <- example_hmmssf_list(n = 25, seed = 13)
  sim <- simulate_generative_hmmssf(
    fit = fit,
    observed_track = fit$observed_track,
    n_sims = 3,
    method = "viterbi",
    parameter_uncertainty = FALSE
  )
  expect_true(all(sim$simulated_states[1, ] == fit$viterbi_path))
  expect_true(all(sim$simulated_states[2, ] == fit$viterbi_path))
  expect_match(paste(sim$warnings, collapse = " "), "diagnostic only")
})
