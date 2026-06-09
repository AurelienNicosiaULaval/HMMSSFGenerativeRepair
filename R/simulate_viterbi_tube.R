#' Sample state paths inside an epsilon Viterbi tube
#'
#' The epsilon Viterbi tube is the set of paths whose log joint likelihood is
#' within `epsilon` of the Viterbi path. This implementation uses dynamic
#' programming max-suffix bounds and stochastic traceback constrained by the
#' tube threshold. It is exact for the max-score constraint used here, but it
#' samples stochastically rather than enumerating all paths.
#'
#' @param fit Adapted object or object accepted by [as_gmov_hmmssf()].
#' @param observed_track Optional observed track.
#' @param epsilon Tube tolerance. If `NULL`, a default is chosen and reported.
#' @param n_sims Number of paths to sample.
#' @param seed Optional random seed.
#' @param adapter Adapter choice.
#' @param max_attempts Maximum traceback attempts.
#'
#' @return A list containing sampled states and path diagnostics.
#' @export
simulate_viterbi_tube_states <- function(
  fit,
  observed_track = NULL,
  epsilon = NULL,
  n_sims = 99,
  seed = NULL,
  adapter = c("auto", "hmmSSF", "list"),
  max_attempts = 1000
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  components <- if (inherits(fit, "gmov_hmmssf_adapter")) fit else as_gmov_hmmssf(fit, observed_track, adapter)
  log_delta <- components$log_delta
  log_gamma <- components$log_gamma
  log_emission <- components$log_emission
  validate_log_hmm_components(log_delta, log_gamma, log_emission)

  n_steps <- nrow(log_emission)
  n_states <- length(log_delta)
  if (is.null(epsilon)) {
    epsilon <- max(1, 0.01 * n_steps)
    message("`epsilon` was NULL; using epsilon = ", signif(epsilon, 4), " for ", n_steps, " time steps.")
  }

  forward <- viterbi_forward_scores(log_delta, log_gamma, log_emission)
  backward <- viterbi_backward_scores(log_gamma, log_emission)
  v_path <- viterbi_path_from_components(log_delta, log_gamma, log_emission)
  v_score <- score_state_path(v_path, log_delta, log_gamma, log_emission)
  threshold <- v_score - epsilon

  paths <- matrix(NA_integer_, nrow = n_sims, ncol = n_steps)
  scores <- rep(NA_real_, n_sims)
  accepted <- 0L
  attempts <- 0L
  while (accepted < n_sims && attempts < max_attempts) {
    attempts <- attempts + 1L
    path <- sample_one_tube_path(log_delta, log_gamma, log_emission, forward, backward, threshold)
    score <- score_state_path(path, log_delta, log_gamma, log_emission)
    if (is.finite(score) && score >= threshold - 1e-8) {
      accepted <- accepted + 1L
      paths[accepted, ] <- path
      scores[accepted] <- score
    }
  }

  warnings <- character()
  if (accepted < n_sims) {
    warnings <- paste0("Only ", accepted, " paths were sampled inside the epsilon tube after ", attempts, " attempts.")
    paths <- paths[seq_len(max(accepted, 1L)), , drop = FALSE]
    scores <- scores[seq_len(max(accepted, 1L))]
    if (accepted == 0L) {
      paths <- matrix(v_path, nrow = 1L)
      scores <- v_score
    }
  }

  occupancy <- do.call(
    rbind,
    lapply(seq_len(nrow(paths)), function(i) {
      state_occupancy_vector(paths[i, ], n_states = n_states)
    })
  )
  colnames(occupancy) <- paste0("state_", seq_len(n_states))

  list(
    states = paths,
    method = "viterbi_tube",
    epsilon = epsilon,
    viterbi_path = v_path,
    viterbi_score = v_score,
    path_diagnostics = data.frame(
      path_id = seq_len(nrow(paths)),
      log_likelihood_path = scores,
      delta_from_viterbi = v_score - scores,
      number_of_switches = rowSums(paths[, -1L, drop = FALSE] != paths[, -ncol(paths), drop = FALSE])
    ),
    occupancy = occupancy,
    warnings = warnings
  )
}

viterbi_forward_scores <- function(log_delta, log_gamma, log_emission) {
  n_steps <- nrow(log_emission)
  n_states <- length(log_delta)
  f <- matrix(-Inf, nrow = n_steps, ncol = n_states)
  f[1L, ] <- log_delta + log_emission[1L, ]
  if (n_steps >= 2L) {
    for (t in 2:n_steps) {
      for (j in seq_len(n_states)) {
        f[t, j] <- max(f[t - 1L, ] + log_gamma[, j, t - 1L]) + log_emission[t, j]
      }
    }
  }
  f
}

viterbi_backward_scores <- function(log_gamma, log_emission) {
  n_steps <- nrow(log_emission)
  n_states <- ncol(log_emission)
  b <- matrix(0, nrow = n_steps, ncol = n_states)
  if (n_steps >= 2L) {
    for (t in seq(n_steps - 1L, 1L, by = -1L)) {
      for (i in seq_len(n_states)) {
        b[t, i] <- max(log_gamma[i, , t] + log_emission[t + 1L, ] + b[t + 1L, ])
      }
    }
  }
  b
}

sample_one_tube_path <- function(log_delta, log_gamma, log_emission, forward, backward, threshold) {
  n_steps <- nrow(log_emission)
  n_states <- length(log_delta)
  path <- integer(n_steps)
  start_scores <- forward[1L, ] + backward[1L, ]
  admissible <- which(start_scores >= threshold - 1e-10)
  if (length(admissible) == 0L) {
    admissible <- which.max(start_scores)
  }
  path[1L] <- admissible[sample_log_weights(start_scores[admissible])]
  prefix_score <- log_delta[path[1L]] + log_emission[1L, path[1L]]

  if (n_steps >= 2L) {
    for (t in seq_len(n_steps - 1L)) {
      candidate_scores <- prefix_score + log_gamma[path[t], , t] + log_emission[t + 1L, ] + backward[t + 1L, ]
      admissible <- which(candidate_scores >= threshold - 1e-10)
      if (length(admissible) == 0L) {
        admissible <- which.max(candidate_scores)
      }
      next_state <- admissible[sample_log_weights(candidate_scores[admissible])]
      prefix_score <- prefix_score + log_gamma[path[t], next_state, t] + log_emission[t + 1L, next_state]
      path[t + 1L] <- next_state
    }
  }
  path
}
