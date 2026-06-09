# Utility helpers for HMM-SSF generative validation. These helpers are
# intentionally small and dependency-light.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

check_probability_vector <- function(x, name = "probability vector") {
  if (!is.numeric(x) || length(x) < 1L || any(!is.finite(x)) || any(x < 0)) {
    stop("`", name, "` must be a finite non-negative numeric vector.", call. = FALSE)
  }
  total <- sum(x)
  if (!is.finite(total) || total <= 0) {
    stop("`", name, "` must have positive total mass.", call. = FALSE)
  }
  x / total
}

check_transition_matrix <- function(x, name = "transition matrix") {
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) != ncol(x)) {
    stop("`", name, "` must be a square numeric matrix.", call. = FALSE)
  }
  if (any(!is.finite(x)) || any(x < 0)) {
    stop("`", name, "` must contain finite non-negative values.", call. = FALSE)
  }
  rs <- rowSums(x)
  if (any(rs <= 0)) {
    stop("Every row of `", name, "` must have positive total mass.", call. = FALSE)
  }
  x / rs
}

safe_log <- function(x) {
  log(pmax(x, .Machine$double.xmin))
}

log_sum_exp <- function(x) {
  if (length(x) == 0L) {
    return(-Inf)
  }
  m <- max(x)
  if (!is.finite(m)) {
    return(m)
  }
  m + log(sum(exp(x - m)))
}

normalize_log_weights <- function(log_w) {
  log_w - log_sum_exp(log_w)
}

sample_log_weights <- function(log_w) {
  probs <- exp(normalize_log_weights(log_w))
  sample(seq_along(probs), size = 1L, prob = probs)
}

as_track_df <- function(track) {
  if (!is.data.frame(track)) {
    stop("`observed_track` must be a data frame.", call. = FALSE)
  }
  if (all(c("x", "y") %in% names(track))) {
    out <- data.frame(x = as.numeric(track$x), y = as.numeric(track$y))
  } else if (all(c("x_", "y_") %in% names(track))) {
    out <- data.frame(x = as.numeric(track$x_), y = as.numeric(track$y_))
  } else {
    stop("Track data must contain either `x`, `y` or `x_`, `y_`.", call. = FALSE)
  }
  keep <- is.finite(out$x) & is.finite(out$y)
  if (!all(keep)) {
    warning("Removed rows with missing or non-finite coordinates.", call. = FALSE)
    out <- out[keep, , drop = FALSE]
  }
  if (nrow(out) < 2L) {
    stop("A track must contain at least two finite locations.", call. = FALSE)
  }
  out
}

track_step_lengths <- function(track) {
  track <- as_track_df(track)
  sqrt(diff(track$x)^2 + diff(track$y)^2)
}

track_bearings <- function(track) {
  track <- as_track_df(track)
  atan2(diff(track$y), diff(track$x))
}

wrap_angle <- function(x) {
  atan2(sin(x), cos(x))
}

track_turning_angles <- function(track) {
  bearings <- track_bearings(track)
  if (length(bearings) < 2L) {
    return(numeric())
  }
  wrap_angle(diff(bearings))
}

straightness_index <- function(track) {
  track <- as_track_df(track)
  path_length <- sum(track_step_lengths(track), na.rm = TRUE)
  if (!is.finite(path_length) || path_length <= 0) {
    return(NA_real_)
  }
  displacement <- sqrt((track$x[nrow(track)] - track$x[1])^2 + (track$y[nrow(track)] - track$y[1])^2)
  displacement / path_length
}

msd_curve <- function(track, max_lag = NULL) {
  track <- as_track_df(track)
  n <- nrow(track)
  if (is.null(max_lag)) {
    max_lag <- max(1L, floor((n - 1L) / 2L))
  }
  max_lag <- min(as.integer(max_lag), n - 1L)
  lag <- seq_len(max_lag)
  msd <- vapply(lag, function(k) {
    mean((track$x[(k + 1L):n] - track$x[1L:(n - k)])^2 +
      (track$y[(k + 1L):n] - track$y[1L:(n - k)])^2)
  }, numeric(1L))
  data.frame(lag = lag, msd = msd)
}

run_lengths_by_state <- function(states) {
  states <- as.integer(states)
  if (length(states) == 0L) {
    return(data.frame(state = integer(), run_length = integer()))
  }
  r <- rle(states)
  data.frame(state = as.integer(r$values), run_length = as.integer(r$lengths))
}

state_occupancy_vector <- function(states, n_states = NULL) {
  states <- as.integer(states)
  if (is.null(n_states)) {
    n_states <- max(states, na.rm = TRUE)
  }
  tab <- tabulate(states, nbins = n_states)
  if (sum(tab) == 0L) {
    return(rep(NA_real_, n_states))
  }
  tab / sum(tab)
}

transition_count_matrix <- function(states, n_states = NULL) {
  states <- as.integer(states)
  if (is.null(n_states)) {
    n_states <- max(states, na.rm = TRUE)
  }
  out <- matrix(0, nrow = n_states, ncol = n_states)
  if (length(states) < 2L) {
    return(out)
  }
  for (i in seq_len(length(states) - 1L)) {
    out[states[i], states[i + 1L]] <- out[states[i], states[i + 1L]] + 1L
  }
  out
}

switching_rate <- function(states) {
  states <- as.integer(states)
  if (length(states) < 2L) {
    return(NA_real_)
  }
  mean(states[-1L] != states[-length(states)])
}

wasserstein_1d <- function(x, y, probs = seq(0.05, 0.95, by = 0.05)) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) == 0L || length(y) == 0L) {
    return(NA_real_)
  }
  mean(abs(stats::quantile(x, probs = probs, names = FALSE, type = 8) -
    stats::quantile(y, probs = probs, names = FALSE, type = 8)))
}

mc_rank_p_value <- function(observed, simulated, alternative = c("greater", "less", "two.sided"), tol = 1e-12) {
  alternative <- match.arg(alternative)
  simulated <- simulated[is.finite(simulated)]
  if (!is.finite(observed) || length(simulated) == 0L) {
    return(NA_real_)
  }
  if (identical(alternative, "greater")) {
    return((sum(simulated >= observed - tol) + 1) / (length(simulated) + 1))
  }
  if (identical(alternative, "less")) {
    return((sum(simulated <= observed + tol) + 1) / (length(simulated) + 1))
  }
  center <- stats::median(simulated)
  obs_dev <- max(0, abs(observed - center) - tol)
  sim_dev <- pmax(0, abs(simulated - center) - tol)
  (sum(sim_dev >= obs_dev) + 1) / (length(simulated) + 1)
}

circular_summary <- function(angles) {
  angles <- angles[is.finite(angles)]
  if (length(angles) == 0L) {
    return(list(mean_angle = NA_real_, mean_resultant_length = NA_real_))
  }
  cbar <- mean(cos(angles))
  sbar <- mean(sin(angles))
  list(
    mean_angle = atan2(sbar, cbar),
    mean_resultant_length = sqrt(cbar^2 + sbar^2)
  )
}

make_warning <- function(...) {
  paste(..., collapse = "")
}
