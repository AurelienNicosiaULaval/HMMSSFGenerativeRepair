test_that("sharpness classifies diffuse simulations as diffuse or uninformative", {
  sharp <- compute_sharpness_scalar(c(-100, 0, 100), c = 1)
  label <- classify_sharpness(sharp$value, thresholds = c(sharp = 0.1, moderate = 0.5, diffuse = 1))
  expect_true(label %in% c("diffuse", "uninformative"))
})

test_that("curve sharpness returns relative envelope area", {
  curves <- cbind(c(1, 2, 3), c(2, 3, 4), c(3, 4, 5))
  sharp <- compute_sharpness_curve(curves)
  expect_true(is.finite(sharp$value))
  expect_true(sharp$envelope_area > 0)
})
