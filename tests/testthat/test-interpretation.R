test_that("interpretation labels follow compatibility and sharpness rules", {
  df <- data.frame(
    method = "markov",
    diagnostic = "msd",
    state = NA_integer_,
    observed_statistic = 1,
    simulated_median = 1,
    simulated_q025 = 0,
    simulated_q975 = 2,
    p_value = c(0.2, 0.2, 0.01),
    sharpness_value = c(0.1, 5, 0.1),
    sharpness_label = c("sharp", "uninformative", "sharp"),
    interpretation_label = NA_character_,
    n_sims = 99,
    n_effective = 99,
    warning = "",
    stringsAsFactors = FALSE
  )
  out <- interpret_hmmssf_diagnostics(df)
  expect_equal(out$interpretation_label, c("pass_sharp", "pass_diffuse", "fail"))
})

test_that("method comparison returns transition-dynamics interpretation", {
  df <- data.frame(
    method = c("markov", "viterbi"),
    diagnostic = "msd",
    state = NA_integer_,
    interpretation_label = c("fail", "pass_sharp")
  )
  out <- compare_simulation_methods(df)
  expect_match(out$comparison_message, "latent transition")
})
