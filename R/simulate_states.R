#' Simulate Markov state sequences
#'
#' @param initial Initial state probabilities.
#' @param transition Transition matrix, transition array, or transition function.
#' @param n_steps Number of time steps.
#' @param n_sims Number of state paths.
#' @param seed Optional random seed.
#' @param ... Passed to a transition function.
#'
#' @return Integer matrix with one row per simulation and one column per time.
#' @export
simulate_markov_states <- function(initial, transition, n_steps, n_sims = 99, seed = NULL, ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  initial <- check_probability_vector(initial, "initial")
  n_states <- length(initial)
  n_steps <- as.integer(n_steps)
  n_sims <- as.integer(n_sims)
  if (n_steps < 1L || n_sims < 1L) {
    stop("`n_steps` and `n_sims` must be positive integers.", call. = FALSE)
  }

  out <- matrix(NA_integer_, nrow = n_sims, ncol = n_steps)
  for (sim in seq_len(n_sims)) {
    out[sim, 1L] <- sample(seq_len(n_states), 1L, prob = initial)
    if (n_steps >= 2L) {
      for (t in seq_len(n_steps - 1L)) {
        gamma_t <- transition_at(transition, t, current_state = out[sim, t], ...)
        out[sim, t + 1L] <- sample(seq_len(n_states), 1L, prob = gamma_t[out[sim, t], ])
      }
    }
  }
  out
}

#' Return a fixed Viterbi state path
#'
#' @param fit Adapted object or object accepted by [as_gmov_hmmssf()].
#' @param observed_track Optional observed track.
#' @param n_sims Number of rows to return.
#' @param adapter Adapter choice.
#'
#' @return Integer matrix with repeated Viterbi paths.
#' @export
simulate_viterbi_states <- function(fit, observed_track = NULL, n_sims = 1, adapter = c("auto", "hmmSSF", "list")) {
  components <- if (inherits(fit, "gmov_hmmssf_adapter")) fit else as_gmov_hmmssf(fit, observed_track, adapter)
  path <- components$viterbi_path
  if (is.null(path)) {
    path <- viterbi_path_from_components(components$log_delta, components$log_gamma, components$log_emission)
  }
  matrix(rep(as.integer(path), times = n_sims), nrow = n_sims, byrow = TRUE)
}

transition_at <- function(transition, t, current_state = NULL, ...) {
  if (is.function(transition)) {
    gamma <- transition(t = t, current_state = current_state, ...)
    return(check_transition_matrix(gamma, "transition function result"))
  }
  if (is.matrix(transition)) {
    return(check_transition_matrix(transition, "transition"))
  }
  if (is.array(transition) && length(dim(transition)) == 3L) {
    idx <- min(t, dim(transition)[3L])
    return(check_transition_matrix(transition[, , idx], "transition array slice"))
  }
  stop("`transition` must be a matrix, array, or function.", call. = FALSE)
}

viterbi_path_from_components <- function(log_delta, log_gamma, log_emission) {
  validate_log_hmm_components(log_delta, log_gamma, log_emission)
  n_steps <- nrow(log_emission)
  n_states <- length(log_delta)
  delta <- matrix(-Inf, nrow = n_steps, ncol = n_states)
  psi <- matrix(NA_integer_, nrow = n_steps, ncol = n_states)
  delta[1L, ] <- log_delta + log_emission[1L, ]
  if (n_steps >= 2L) {
    for (t in 2:n_steps) {
      gamma_t <- matrix(log_gamma[, , t - 1L, drop = FALSE], nrow = n_states, ncol = n_states)
      for (j in seq_len(n_states)) {
        scores <- delta[t - 1L, ] + gamma_t[, j]
        psi[t, j] <- which.max(scores)
        delta[t, j] <- max(scores) + log_emission[t, j]
      }
    }
  }
  path <- integer(n_steps)
  path[n_steps] <- which.max(delta[n_steps, ])
  if (n_steps >= 2L) {
    for (t in seq(n_steps, 2L, by = -1L)) {
      path[t - 1L] <- psi[t, path[t]]
    }
  }
  path
}

score_state_path <- function(states, log_delta, log_gamma, log_emission) {
  validate_log_hmm_components(log_delta, log_gamma, log_emission)
  states <- as.integer(states)
  if (length(states) != nrow(log_emission)) {
    stop("State path length does not match `log_emission`.", call. = FALSE)
  }
  score <- log_delta[states[1L]] + log_emission[1L, states[1L]]
  if (length(states) >= 2L) {
    for (t in seq_len(length(states) - 1L)) {
      score <- score + log_gamma[states[t], states[t + 1L], t] + log_emission[t + 1L, states[t + 1L]]
    }
  }
  score
}

validate_log_hmm_components <- function(log_delta, log_gamma, log_emission) {
  if (is.null(log_delta) || is.null(log_gamma) || is.null(log_emission)) {
    stop("`log_delta`, `log_gamma`, and `log_emission` are required.", call. = FALSE)
  }
  if (!is.matrix(log_emission)) {
    log_emission <- as.matrix(log_emission)
  }
  if (length(log_delta) != ncol(log_emission)) {
    stop("`log_delta` length must match number of emission states.", call. = FALSE)
  }
  if (!is.array(log_gamma) || length(dim(log_gamma)) != 3L) {
    stop("`log_gamma` must be a 3D log transition array.", call. = FALSE)
  }
  if (dim(log_gamma)[1L] != length(log_delta) || dim(log_gamma)[2L] != length(log_delta)) {
    stop("`log_gamma` dimensions must match the number of states.", call. = FALSE)
  }
  if (dim(log_gamma)[3L] < max(1L, nrow(log_emission) - 1L)) {
    stop("`log_gamma` has too few time slices.", call. = FALSE)
  }
  invisible(TRUE)
}
