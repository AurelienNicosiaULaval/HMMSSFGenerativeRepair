#' Compute state occupancy diagnostics
#'
#' @param observed_states Observed or decoded state path.
#' @param simulated_states Matrix of simulated state paths.
#'
#' @return Diagnostic rows by state and total variation.
#' @export
diagnostic_state_occupancy <- function(observed_states, simulated_states) {
  simulated_states <- as.matrix(simulated_states)
  n_states <- max(c(observed_states, simulated_states), na.rm = TRUE)
  obs <- state_occupancy_vector(observed_states, n_states)
  sim <- do.call(
    rbind,
    lapply(seq_len(nrow(simulated_states)), function(i) {
      state_occupancy_vector(simulated_states[i, ], n_states = n_states)
    })
  )
  rows <- vector("list", n_states + 1L)
  for (s in seq_len(n_states)) {
    obs_stat <- abs(obs[s] - mean(sim[, s], na.rm = TRUE))
    sim_stat <- leave_one_out_scalar_discrepancy(sim[, s])
    sharp <- compute_sharpness_scalar(sim[, s])
    rows[[s]] <- diagnostic_row(
      diagnostic = "state_occupancy",
      state = s,
      observed = obs_stat,
      simulated = sim_stat,
      sharpness = sharp,
      alternative = "greater"
    )
  }
  tv_obs <- mean(abs(obs - colMeans(sim, na.rm = TRUE)))
  tv_sim <- vapply(seq_len(nrow(sim)), function(i) {
    mean(abs(sim[i, ] - colMeans(sim[-i, , drop = FALSE], na.rm = TRUE)))
  }, numeric(1L))
  rows[[n_states + 1L]] <- diagnostic_row("state_occupancy_total_variation", NA, tv_obs, tv_sim, compute_sharpness_scalar(tv_sim), "greater")
  out <- do.call(rbind, rows)
  out$state_reference <- state_reference_note()
  out
}

#' Compute residence-time diagnostics
#'
#' @param observed_states Observed or decoded state path.
#' @param simulated_states Matrix of simulated state paths.
#'
#' @return Diagnostic rows by state.
#' @export
diagnostic_state_residence_time <- function(observed_states, simulated_states) {
  simulated_states <- as.matrix(simulated_states)
  n_states <- max(c(observed_states, simulated_states), na.rm = TRUE)
  obs_runs <- run_lengths_by_state(observed_states)
  rows <- vector("list", n_states)
  for (s in seq_len(n_states)) {
    obs_s <- obs_runs$run_length[obs_runs$state == s]
    sim_s <- lapply(seq_len(nrow(simulated_states)), function(i) {
      r <- run_lengths_by_state(simulated_states[i, ])
      r$run_length[r$state == s]
    })
    comparison <- distribution_leave_one_out(obs_s, sim_s)
    warn <- if (length(obs_s) < 3L) "Too few observed residence bouts; low power." else ""
    rows[[s]] <- diagnostic_row("state_residence_time", s, comparison$observed, comparison$simulated, compute_sharpness_scalar(comparison$simulated), "greater", warning = warn)
  }
  out <- do.call(rbind, rows)
  out$state_reference <- state_reference_note()
  out
}

#' Compute geometric residence-time diagnostics
#'
#' This diagnostic compares decoded or supplied residence times with the
#' geometric dwell-time distribution implied by a homogeneous HMM transition
#' matrix. It is intended as a targeted check of the Markov residence-time
#' assumption, not as proof that decoded states are true states.
#'
#' @param observed_states Observed or decoded state path.
#' @param transition Homogeneous transition matrix.
#' @param n_sims Number of parametric-bootstrap samples.
#' @param seed Optional random seed.
#'
#' @return Diagnostic rows by state.
#' @export
diagnostic_state_residence_geometric <- function(observed_states, transition, n_sims = 999, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  transition <- check_transition_matrix(transition, "transition")
  n_states <- nrow(transition)
  observed_states <- as.integer(observed_states)
  obs_runs <- run_lengths_by_state(observed_states)
  rows <- vector("list", n_states)

  for (s in seq_len(n_states)) {
    obs_s <- obs_runs$run_length[obs_runs$state == s]
    warning <- character()
    if (length(obs_s) == 0L) {
      rows[[s]] <- diagnostic_row(
        "state_residence_geometric",
        s,
        NA_real_,
        numeric(),
        list(value = NA_real_),
        "greater",
        warning = "No observed residence bouts for this state."
      )
      next
    }
    if (length(obs_s) < 3L) {
      warning <- c(warning, "Too few observed residence bouts; low power.")
    }

    leave_prob <- 1 - transition[s, s]
    if (!is.finite(leave_prob) || leave_prob <= .Machine$double.eps) {
      rows[[s]] <- diagnostic_row(
        "state_residence_geometric",
        s,
        NA_real_,
        numeric(),
        list(value = NA_real_),
        "greater",
        warning = paste(
          c(warning, "Self-transition probability is too close to one; geometric bootstrap is undefined."),
          collapse = " "
        )
      )
      next
    }

    sim_s <- lapply(seq_len(as.integer(n_sims)), function(i) {
      stats::rgeom(length(obs_s), prob = leave_prob) + 1L
    })
    comparison <- distribution_leave_one_out(obs_s, sim_s)
    rows[[s]] <- diagnostic_row(
      "state_residence_geometric",
      s,
      comparison$observed,
      comparison$simulated,
      compute_sharpness_scalar(comparison$simulated),
      "greater",
      warning = paste(warning, collapse = " ")
    )
  }

  out <- do.call(rbind, rows)
  out$state_reference <- state_reference_note()
  out
}

#' Compute switching-rate diagnostics
#'
#' @param observed_states Observed or decoded state path.
#' @param simulated_states Matrix of simulated state paths.
#'
#' @return Diagnostic row.
#' @export
diagnostic_switching_rate <- function(observed_states, simulated_states) {
  simulated_states <- as.matrix(simulated_states)
  obs <- switching_rate(observed_states)
  sim <- apply(simulated_states, 1L, switching_rate)
  obs_stat <- abs(obs - mean(sim, na.rm = TRUE))
  sim_stat <- leave_one_out_scalar_discrepancy(sim)
  diagnostic_row("switching_rate", NA, obs_stat, sim_stat, compute_sharpness_scalar(sim), "greater")
}

#' Compute transition-count diagnostics
#'
#' @param observed_states Observed or decoded state path.
#' @param simulated_states Matrix of simulated state paths.
#'
#' @return Diagnostic row.
#' @export
diagnostic_transition_counts <- function(observed_states, simulated_states) {
  simulated_states <- as.matrix(simulated_states)
  n_states <- max(c(observed_states, simulated_states), na.rm = TRUE)
  obs <- transition_count_matrix(observed_states, n_states)
  obs_norm <- obs / max(1, sum(obs))
  sim_vecs <- do.call(
    rbind,
    lapply(seq_len(nrow(simulated_states)), function(i) {
      mat <- transition_count_matrix(simulated_states[i, ], n_states)
      mat_norm <- mat / max(1, sum(mat))
      as.numeric(mat_norm)
    })
  )
  mean_sim <- colMeans(sim_vecs, na.rm = TRUE)
  obs_stat <- sqrt(sum((as.numeric(obs_norm) - mean_sim)^2))
  sim_stats <- vapply(seq_len(nrow(sim_vecs)), function(i) {
    loo <- colMeans(sim_vecs[-i, , drop = FALSE], na.rm = TRUE)
    sqrt(sum((sim_vecs[i, ] - loo)^2))
  }, numeric(1L))
  diagnostic_row("transition_counts", NA, obs_stat, sim_stats, compute_sharpness_scalar(sim_stats), "greater")
}

state_conditioned_step_length <- function(observed_track, observed_states, simulated_tracks, simulated_states) {
  n_states <- max(c(observed_states, simulated_states), na.rm = TRUE)
  obs_steps <- track_step_lengths(observed_track)
  rows <- vector("list", n_states)
  for (s in seq_len(n_states)) {
    obs_idx <- which(observed_states[-1L] == s)
    obs_s <- obs_steps[obs_idx]
    sim_stat <- vapply(seq_along(simulated_tracks), function(i) {
      sim_steps <- track_step_lengths(simulated_tracks[[i]])
      sim_idx <- which(simulated_states[i, -1L] == s)
      wasserstein_1d(obs_s, sim_steps[sim_idx])
    }, numeric(1L))
    sim_dists <- lapply(seq_along(simulated_tracks), function(i) {
      sim_steps <- track_step_lengths(simulated_tracks[[i]])
      sim_idx <- which(simulated_states[i, -1L] == s)
      sim_steps[sim_idx]
    })
    comparison <- distribution_leave_one_out(obs_s, sim_dists)
    warn <- if (length(obs_s) < 5L) "Too few state-specific steps; low power." else ""
    rows[[s]] <- diagnostic_row("state_conditioned_step_length", s, comparison$observed, comparison$simulated, compute_sharpness_scalar(comparison$simulated), "greater", warning = warn)
  }
  do.call(rbind, rows)
}

state_conditioned_turning_angle <- function(observed_track, observed_states, simulated_tracks, simulated_states) {
  n_states <- max(c(observed_states, simulated_states), na.rm = TRUE)
  obs_angles <- track_turning_angles(observed_track)
  rows <- vector("list", n_states)
  for (s in seq_len(n_states)) {
    obs_idx <- which(observed_states[-c(1L, 2L)] == s)
    obs_summary <- circular_summary(obs_angles[obs_idx])
    sim_stat <- vapply(seq_along(simulated_tracks), function(i) {
      sim_angles <- track_turning_angles(simulated_tracks[[i]])
      sim_idx <- which(simulated_states[i, -c(1L, 2L)] == s)
      sim_summary <- circular_summary(sim_angles[sim_idx])
      sim_summary$mean_resultant_length
    }, numeric(1L))
    obs_stat <- abs(obs_summary$mean_resultant_length - mean(sim_stat, na.rm = TRUE))
    sim_ref <- leave_one_out_scalar_discrepancy(sim_stat)
    warn <- if (length(obs_idx) < 5L) "Too few state-specific turning angles; low power." else ""
    rows[[s]] <- diagnostic_row("state_conditioned_turning_angle", s, obs_stat, sim_ref, compute_sharpness_scalar(sim_stat), "greater", warning = warn)
  }
  do.call(rbind, rows)
}

state_conditioned_msd <- function(observed_track, observed_states, simulated_tracks, simulated_states, max_lag = 10) {
  n_states <- max(c(observed_states, simulated_states), na.rm = TRUE)
  rows <- vector("list", n_states)
  for (s in seq_len(n_states)) {
    obs_track_s <- observed_track[observed_states == s, , drop = FALSE]
    if (nrow(obs_track_s) <= max_lag + 1L) {
      rows[[s]] <- diagnostic_row("state_conditioned_msd", s, NA_real_, numeric(), list(value = NA_real_), "greater", warning = "Too few state-specific points; low power.")
      next
    }
    obs_curve <- msd_curve(obs_track_s, max_lag = max_lag)$msd
    sim_curves <- lapply(seq_along(simulated_tracks), function(i) {
      tr <- simulated_tracks[[i]][simulated_states[i, ] == s, , drop = FALSE]
      if (nrow(tr) <= max_lag + 1L) {
        return(rep(NA_real_, max_lag))
      }
      msd_curve(tr, max_lag = max_lag)$msd
    })
    sim_matrix <- do.call(cbind, sim_curves)
    mean_sim <- rowMeans(sim_matrix, na.rm = TRUE)
    obs_stat <- sum((obs_curve - mean_sim)^2)
    sim_stat <- vapply(seq_len(ncol(sim_matrix)), function(i) {
      loo <- rowMeans(sim_matrix[, -i, drop = FALSE], na.rm = TRUE)
      sum((sim_matrix[, i] - loo)^2)
    }, numeric(1L))
    rows[[s]] <- diagnostic_row("state_conditioned_msd", s, obs_stat, sim_stat, compute_sharpness_curve(sim_matrix), "greater")
  }
  do.call(rbind, rows)
}

leave_one_out_scalar_discrepancy <- function(values) {
  values <- as.numeric(values)
  vapply(seq_along(values), function(i) {
    abs(values[i] - mean(values[-i], na.rm = TRUE))
  }, numeric(1L))
}

distribution_leave_one_out <- function(observed, simulated_distributions) {
  observed <- observed[is.finite(observed)]
  simulated_distributions <- lapply(simulated_distributions, function(x) x[is.finite(x)])
  observed_distances <- vapply(simulated_distributions, function(x) {
    wasserstein_1d(observed, x)
  }, numeric(1L))
  observed_stat <- mean(observed_distances, na.rm = TRUE)
  simulated_stats <- vapply(seq_along(simulated_distributions), function(i) {
    others <- simulated_distributions[-i]
    mean(vapply(others, function(x) wasserstein_1d(simulated_distributions[[i]], x), numeric(1L)), na.rm = TRUE)
  }, numeric(1L))
  list(observed = observed_stat, simulated = simulated_stats)
}

diagnostic_row <- function(diagnostic, state, observed, simulated, sharpness, alternative = "greater", warning = "") {
  simulated <- simulated[is.finite(simulated)]
  med <- if (length(simulated) == 0L) NA_real_ else stats::median(simulated)
  qs <- if (length(simulated) == 0L) c(NA_real_, NA_real_) else stats::quantile(simulated, c(0.025, 0.975), names = FALSE, type = 8)
  sharp_value <- sharpness$value %||% NA_real_
  p_value <- mc_rank_p_value(observed, simulated, alternative = alternative)
  p_resolution <- if (length(simulated) == 0L) NA_real_ else 1 / (length(simulated) + 1)
  data.frame(
    method = NA_character_,
    diagnostic = diagnostic,
    state = if (is.na(state)) NA_integer_ else as.integer(state),
    observed_statistic = observed,
    simulated_median = med,
    simulated_q025 = qs[1L],
    simulated_q975 = qs[2L],
    observed_discrepancy = observed,
    simulated_discrepancy_median = med,
    simulated_discrepancy_q025 = qs[1L],
    simulated_discrepancy_q975 = qs[2L],
    p_value = p_value,
    mc_p_value = p_value,
    mc_p_value_resolution = p_resolution,
    statistic_type = "discrepancy",
    state_reference = NA_character_,
    sharpness_value = sharp_value,
    sharpness_label = classify_sharpness(sharp_value),
    interpretation_label = NA_character_,
    n_sims = length(simulated),
    n_effective = length(simulated),
    warning = warning,
    stringsAsFactors = FALSE
  )
}

state_reference_note <- function() {
  paste(
    "Observed states are decoded or supplied;",
    "simulated states are generated latent paths or supplied state paths."
  )
}
