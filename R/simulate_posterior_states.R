#' Sample posterior state paths with FFBS
#'
#' @param fit Adapted object or object accepted by [as_gmov_hmmssf()].
#' @param observed_track Optional observed track.
#' @param n_sims Number of posterior paths.
#' @param seed Optional random seed.
#' @param adapter Adapter choice.
#'
#' @return Integer matrix of sampled posterior state paths.
#' @export
simulate_posterior_states <- function(fit, observed_track = NULL, n_sims = 99, seed = NULL, adapter = c("auto", "hmmSSF", "list")) {
  components <- if (inherits(fit, "gmov_hmmssf_adapter")) fit else as_gmov_hmmssf(fit, observed_track, adapter)
  if (is.null(components$log_delta) || is.null(components$log_gamma) || is.null(components$log_emission)) {
    stop(
      "Posterior state sampling requires initial probabilities, transition probabilities, and state log-likelihoods.",
      call. = FALSE
    )
  }
  simulate_posterior_states_from_components(
    log_delta = components$log_delta,
    log_gamma = components$log_gamma,
    log_emission = components$log_emission,
    n_sims = n_sims,
    seed = seed
  )
}

#' FFBS from explicit HMM components
#'
#' @param log_delta Log initial probabilities.
#' @param log_gamma Log transition probability array.
#' @param log_emission Log state-dependent likelihood matrix.
#' @param n_sims Number of paths.
#' @param seed Optional seed.
#'
#' @return Integer matrix of sampled state paths.
#' @export
simulate_posterior_states_from_components <- function(log_delta, log_gamma, log_emission, n_sims = 99, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  validate_log_hmm_components(log_delta, log_gamma, log_emission)
  n_steps <- nrow(log_emission)
  n_states <- length(log_delta)
  alpha <- matrix(-Inf, nrow = n_steps, ncol = n_states)
  alpha[1L, ] <- normalize_log_weights(log_delta + log_emission[1L, ])
  if (n_steps >= 2L) {
    for (t in 2:n_steps) {
      for (j in seq_len(n_states)) {
        alpha[t, j] <- log_emission[t, j] + log_sum_exp(alpha[t - 1L, ] + log_gamma[, j, t - 1L])
      }
      alpha[t, ] <- normalize_log_weights(alpha[t, ])
    }
  }

  out <- matrix(NA_integer_, nrow = n_sims, ncol = n_steps)
  for (sim in seq_len(n_sims)) {
    out[sim, n_steps] <- sample_log_weights(alpha[n_steps, ])
    if (n_steps >= 2L) {
      for (t in seq(n_steps - 1L, 1L, by = -1L)) {
        log_w <- alpha[t, ] + log_gamma[, out[sim, t + 1L], t]
        out[sim, t] <- sample_log_weights(log_w)
      }
    }
  }
  out
}
