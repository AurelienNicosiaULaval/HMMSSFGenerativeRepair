test_that("missing adapter elements produce informative errors", {
  expect_error(
    as_gmov_hmmssf(list(initial = c(1, 0)), adapter = "list"),
    "initial.*transition.*kernels|transition.*kernels"
  )
})

test_that("posterior sampler returns valid paths", {
  fit <- example_hmmssf_list(n = 15, seed = 31)
  states <- simulate_posterior_states(fit, n_sims = 4, seed = 32)
  expect_equal(dim(states), c(4L, 15L))
  expect_true(all(states >= 1 & states <= 2))
})
