test_that("diagnose_hmmssf returns a complete structured object", {
  fit <- example_hmmssf_list(n = 30, seed = 41)
  diag <- diagnose_hmmssf(
    fit = fit,
    observed_track = fit$observed_track,
    n_sims = 5,
    methods = c("markov", "viterbi", "viterbi_tube", "posterior"),
    diagnostics = c(
      "ud_wasserstein",
      "msd",
      "sinuosity",
      "state_occupancy",
      "state_residence_time",
      "switching_rate",
      "transition_counts",
      "state_conditioned_step_length",
      "state_conditioned_turning_angle"
    ),
    parameter_uncertainty = FALSE,
    epsilon = 5,
    seed = 42
  )
  expect_s3_class(diag, "gmov_hmmssf_diagnostics")
  expect_true(all(c(
    "method",
    "diagnostic",
    "observed_statistic",
    "observed_discrepancy",
    "simulated_median",
    "simulated_discrepancy_median",
    "p_value",
    "mc_p_value",
    "mc_p_value_resolution",
    "statistic_type",
    "sharpness_label",
    "interpretation_label"
  ) %in% names(diag$diagnostics)))
  expect_true(all(diag$diagnostics$statistic_type == "discrepancy"))
  informative <- diag$diagnostics$n_effective > 0
  expect_equal(
    diag$diagnostics$mc_p_value_resolution[informative],
    1 / (diag$diagnostics$n_effective[informative] + 1)
  )
  expect_true(all(c("markov", "viterbi", "viterbi_tube", "posterior") %in% names(diag$simulations)))
})

test_that("diagnose_hmmssf_simulations handles precomputed tracks", {
  fit <- example_hmmssf_list(n = 30, seed = 51)
  sim <- simulate_generative_hmmssf(
    fit = fit,
    observed_track = fit$observed_track,
    n_sims = 5,
    method = "markov",
    parameter_uncertainty = FALSE,
    seed = 52
  )

  diag <- diagnose_hmmssf_simulations(
    observed_track = sim$observed_track,
    simulated_tracks = sim$simulated_tracks,
    observed_states = fit$viterbi_path,
    simulated_states = sim$simulated_states,
    reference_transition = fit$transition,
    label = "markov_model",
    diagnostics = c("ud_wasserstein", "msd", "sinuosity", "state_occupancy", "state_residence_geometric")
  )

  expect_s3_class(diag, "gmov_hmmssf_diagnostics")
  expect_equal(names(diag$simulations), "markov_model")
  expect_true(all(diag$diagnostics$method == "markov_model"))
  expect_true(all(is.finite(diag$diagnostics$p_value)))
  expect_true("state_residence_geometric" %in% diag$diagnostics$diagnostic)
})
