#' Compute sharpness for scalar diagnostics
#'
#' @param simulated Simulated scalar statistics.
#' @param c Stabilizing constant.
#' @param probs Envelope probabilities.
#'
#' @return A list with envelope width, spread ratio, and label-ready value.
#' @export
compute_sharpness_scalar <- function(simulated, c = 1e-8, probs = c(0.025, 0.975)) {
  simulated <- simulated[is.finite(simulated)]
  if (length(simulated) < 2L) {
    return(list(value = NA_real_, envelope_width = NA_real_, simulated_median = NA_real_))
  }
  qs <- stats::quantile(simulated, probs = probs, names = FALSE, type = 8)
  med <- stats::median(simulated)
  width <- diff(qs)
  list(
    value = width / (abs(med) + c),
    envelope_width = width,
    simulated_median = med,
    simulated_q025 = qs[1L],
    simulated_q975 = qs[2L]
  )
}

#' Compute sharpness for curve diagnostics
#'
#' @param simulated_curves Matrix with one curve per column.
#' @param c Stabilizing constant.
#' @param probs Envelope probabilities.
#'
#' @return A list with envelope area and relative envelope area.
#' @export
compute_sharpness_curve <- function(simulated_curves, c = 1e-8, probs = c(0.025, 0.975)) {
  simulated_curves <- as.matrix(simulated_curves)
  if (ncol(simulated_curves) < 2L) {
    return(list(value = NA_real_, envelope_area = NA_real_))
  }
  lo <- apply(simulated_curves, 1L, stats::quantile, probs = probs[1L], na.rm = TRUE, type = 8)
  hi <- apply(simulated_curves, 1L, stats::quantile, probs = probs[2L], na.rm = TRUE, type = 8)
  med <- apply(simulated_curves, 1L, stats::median, na.rm = TRUE)
  area <- sum(hi - lo, na.rm = TRUE)
  list(
    value = area / (sum(abs(med), na.rm = TRUE) + c),
    envelope_area = area,
    median_curve = med,
    lo = lo,
    hi = hi
  )
}

#' Compute sharpness for distribution diagnostics
#'
#' @param simulated_distributions List of numeric vectors.
#' @param probs Probability grid.
#' @param c Stabilizing constant.
#'
#' @return A list with average quantile band width and relative width.
#' @export
compute_sharpness_distribution <- function(
  simulated_distributions,
  probs = seq(0.05, 0.95, by = 0.05),
  c = 1e-8
) {
  simulated_distributions <- lapply(simulated_distributions, function(x) x[is.finite(x)])
  simulated_distributions <- simulated_distributions[vapply(simulated_distributions, length, integer(1L)) > 0L]
  if (length(simulated_distributions) < 2L) {
    return(list(value = NA_real_, average_quantile_band_width = NA_real_))
  }
  qmat <- do.call(
    cbind,
    lapply(simulated_distributions, stats::quantile, probs = probs, names = FALSE, type = 8)
  )
  band_width <- apply(qmat, 1L, diff_range)
  med <- apply(qmat, 1L, stats::median)
  list(
    value = mean(band_width) / (mean(abs(med)) + c),
    average_quantile_band_width = mean(band_width),
    median_quantiles = med
  )
}

diff_range <- function(x) {
  diff(range(x, na.rm = TRUE))
}

#' Classify sharpness
#'
#' @param value Sharpness value.
#' @param thresholds Named vector with `sharp`, `moderate`, and `diffuse`
#'   thresholds. Larger values are less sharp.
#'
#' @return Sharpness label.
#' @export
classify_sharpness <- function(
  value,
  thresholds = c(sharp = 0.5, moderate = 1.5, diffuse = 3)
) {
  if (!is.finite(value)) {
    return("not_available")
  }
  if (value <= thresholds["sharp"]) {
    return("sharp")
  }
  if (value <= thresholds["moderate"]) {
    return("moderate")
  }
  if (value <= thresholds["diffuse"]) {
    return("diffuse")
  }
  "uninformative"
}
