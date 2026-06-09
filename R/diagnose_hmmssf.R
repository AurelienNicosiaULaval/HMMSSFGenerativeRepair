#' Diagnose HMM-SSF generative simulations
#'
#' Runs simulation methods, global movement diagnostics, HMM-specific state
#' diagnostics, sharpness diagnostics, and compatibility-sharpness
#' interpretations. The function returns a dashboard object rather than a
#' single global goodness-of-fit score.
#'
#' @param fit Fitted object or generic list object.
#' @param observed_track Observed track.
#' @param n_sims Number of simulations per method.
#' @param methods Simulation methods.
#' @param diagnostics Diagnostic names.
#' @param sharpness Should sharpness be computed?
#' @param parameter_uncertainty Should parameter uncertainty be requested?
#' @param epsilon Viterbi tube tolerance.
#' @param barriers Optional barriers for future barrier diagnostics.
#' @param raster_template Optional raster template for future UD diagnostics.
#' @param seed Optional seed.
#' @param ... Passed to simulation.
#'
#' @return Object of class `"gmov_hmmssf_diagnostics"`.
#' @export
diagnose_hmmssf <- function(
  fit,
  observed_track,
  n_sims = 99,
  methods = c("markov", "viterbi", "viterbi_tube", "posterior"),
  diagnostics = c(
    "ud_wasserstein",
    "msd",
    "sinuosity",
    "barrier_crossing",
    "state_occupancy",
    "state_residence_time",
    "switching_rate",
    "transition_counts",
    "state_conditioned_msd",
    "state_conditioned_step_length",
    "state_conditioned_turning_angle",
    "state_conditioned_ud",
    "habitat_use_by_state"
  ),
  sharpness = TRUE,
  parameter_uncertainty = TRUE,
  epsilon = NULL,
  barriers = NULL,
  raster_template = NULL,
  seed = NULL,
  ...
) {
  observed_track <- as_track_df(observed_track)
  methods <- unique(methods)
  diagnostics <- unique(diagnostics)
  simulations <- list()
  diagnostic_tables <- list()
  warnings <- character()

  components <- as_gmov_hmmssf(fit, observed_track = observed_track)
  observed_states <- components$viterbi_path
  if (is.null(observed_states)) {
    observed_states <- viterbi_path_from_components(components$log_delta, components$log_gamma, components$log_emission)
  }

  for (method_i in methods) {
    sim_i <- simulate_generative_hmmssf(
      fit = components,
      observed_track = observed_track,
      n_sims = n_sims,
      method = method_i,
      epsilon = epsilon,
      parameter_uncertainty = parameter_uncertainty,
      seed = seed,
      ...
    )
    simulations[[method_i]] <- sim_i
    warnings <- c(warnings, sim_i$warnings)
    tab_i <- compute_diagnostics_for_simulation(
      sim_i = sim_i,
      observed_states = observed_states,
      diagnostics = diagnostics,
      barriers = barriers
    )
    tab_i$method <- method_i
    diagnostic_tables[[method_i]] <- tab_i
  }

  table <- do.call(rbind, diagnostic_tables)
  rownames(table) <- NULL
  table <- interpret_hmmssf_diagnostics(table)
  comparison <- compare_simulation_methods(table)
  out <- list(
    diagnostics = table,
    simulations = simulations,
    method_comparison = comparison,
    warnings = unique(warnings),
    metadata = list(n_sims = n_sims, methods = methods, requested_diagnostics = diagnostics, sharpness = sharpness),
    plotting_data = list()
  )
  class(out) <- c("gmov_hmmssf_diagnostics", "list")
  out
}

#' Diagnose precomputed HMM-SSF simulations
#'
#' This helper is useful when trajectories are simulated outside the package,
#' for example with a project-specific HMM-SSF simulator. It computes the same
#' diagnostic table and dashboard structure as [diagnose_hmmssf()] without
#' refitting or resimulating the model.
#'
#' @param observed_track Observed track.
#' @param simulated_tracks List of simulated tracks.
#' @param observed_states Optional observed or decoded state path.
#' @param simulated_states Optional matrix of simulated state paths.
#' @param label Label used for the simulation source in plots and tables.
#' @param diagnostics Diagnostic names.
#' @param barriers Optional barriers for future barrier diagnostics.
#' @param metadata Optional metadata added to the output object.
#' @param warnings Optional warnings added to the output object.
#'
#' @return Object of class `"gmov_hmmssf_diagnostics"`.
#' @export
diagnose_hmmssf_simulations <- function(
  observed_track,
  simulated_tracks,
  observed_states = NULL,
  simulated_states = NULL,
  label = "simulation",
  diagnostics = c(
    "ud_wasserstein",
    "msd",
    "sinuosity",
    "barrier_crossing",
    "state_occupancy",
    "state_residence_time",
    "switching_rate",
    "transition_counts",
    "state_conditioned_msd",
    "state_conditioned_step_length",
    "state_conditioned_turning_angle",
    "state_conditioned_ud",
    "habitat_use_by_state"
  ),
  barriers = NULL,
  metadata = list(),
  warnings = character()
) {
  observed_track <- as_track_df(observed_track)
  if (!is.list(simulated_tracks) || length(simulated_tracks) == 0L) {
    stop("`simulated_tracks` must be a non-empty list of track data frames.", call. = FALSE)
  }
  simulated_tracks <- lapply(simulated_tracks, as_track_df)
  names(simulated_tracks) <- names(simulated_tracks) %||% paste0("sim_", seq_along(simulated_tracks))

  state_diagnostics <- c(
    "state_occupancy",
    "state_residence_time",
    "switching_rate",
    "transition_counts",
    "state_conditioned_msd",
    "state_conditioned_step_length",
    "state_conditioned_turning_angle",
    "state_conditioned_ud",
    "habitat_use_by_state"
  )
  diagnostics <- unique(diagnostics)
  requested_state_diagnostics <- intersect(diagnostics, state_diagnostics)

  if (is.null(observed_states) || is.null(simulated_states)) {
    if (length(requested_state_diagnostics) > 0L) {
      warnings <- unique(c(
        warnings,
        "State diagnostics were skipped because `observed_states` or `simulated_states` was not supplied."
      ))
    }
    diagnostics <- setdiff(diagnostics, state_diagnostics)
  } else {
    observed_states <- as.integer(observed_states)
    simulated_states <- as.matrix(simulated_states)
    if (nrow(simulated_states) != length(simulated_tracks)) {
      stop("`simulated_states` must have one row per simulated track.", call. = FALSE)
    }
    if (ncol(simulated_states) != nrow(observed_track)) {
      stop("`simulated_states` must have one column per observed location.", call. = FALSE)
    }
    if (length(observed_states) != nrow(observed_track)) {
      stop("`observed_states` must have one value per observed location.", call. = FALSE)
    }
  }

  sim_object <- list(
    observed_track = observed_track,
    simulated_tracks = simulated_tracks,
    simulated_states = simulated_states,
    method = label,
    model_parameters = NULL,
    parameter_draws = NULL,
    state_paths = NULL,
    metadata = c(
      list(
        n_sims = length(simulated_tracks),
        n_steps = nrow(observed_track),
        simulation_label = label
      ),
      metadata
    ),
    warnings = warnings,
    seed = NA_integer_
  )
  class(sim_object) <- c("gmov_hmmssf_simulations", "list")

  table <- compute_diagnostics_for_simulation(
    sim_i = sim_object,
    observed_states = observed_states,
    diagnostics = diagnostics,
    barriers = barriers
  )
  table$method <- label
  table <- interpret_hmmssf_diagnostics(table)

  out <- list(
    diagnostics = table,
    simulations = stats::setNames(list(sim_object), label),
    method_comparison = compare_simulation_methods(table),
    warnings = unique(warnings),
    metadata = c(
      list(
        n_sims = length(simulated_tracks),
        methods = label,
        requested_diagnostics = diagnostics,
        source = "precomputed_simulations"
      ),
      metadata
    ),
    plotting_data = list()
  )
  class(out) <- c("gmov_hmmssf_diagnostics", "list")
  out
}

compute_diagnostics_for_simulation <- function(sim_i, observed_states, diagnostics, barriers = NULL) {
  rows <- list()
  observed_track <- sim_i$observed_track
  simulated_tracks <- sim_i$simulated_tracks
  simulated_states <- sim_i$simulated_states

  if ("ud_wasserstein" %in% diagnostics) {
    rows[[length(rows) + 1L]] <- global_ud_diagnostic(observed_track, simulated_tracks)
  }
  if ("msd" %in% diagnostics) {
    rows[[length(rows) + 1L]] <- global_msd_diagnostic(observed_track, simulated_tracks)
  }
  if ("sinuosity" %in% diagnostics) {
    rows[[length(rows) + 1L]] <- global_sinuosity_diagnostic(observed_track, simulated_tracks)
  }
  if ("barrier_crossing" %in% diagnostics && !is.null(barriers)) {
    rows[[length(rows) + 1L]] <- diagnostic_row("barrier_crossing", NA, NA_real_, numeric(), list(value = NA_real_), warning = "Barrier diagnostics require a geometry-specific implementation.")
  } else if ("barrier_crossing" %in% diagnostics) {
    rows[[length(rows) + 1L]] <- diagnostic_row("barrier_crossing", NA, NA_real_, numeric(), list(value = NA_real_), warning = "No barrier geometry was supplied.")
  }
  if ("state_occupancy" %in% diagnostics) {
    rows[[length(rows) + 1L]] <- diagnostic_state_occupancy(observed_states, simulated_states)
  }
  if ("state_residence_time" %in% diagnostics) {
    rows[[length(rows) + 1L]] <- diagnostic_state_residence_time(observed_states, simulated_states)
  }
  if ("switching_rate" %in% diagnostics) {
    rows[[length(rows) + 1L]] <- diagnostic_switching_rate(observed_states, simulated_states)
  }
  if ("transition_counts" %in% diagnostics) {
    rows[[length(rows) + 1L]] <- diagnostic_transition_counts(observed_states, simulated_states)
  }
  if ("state_conditioned_step_length" %in% diagnostics) {
    rows[[length(rows) + 1L]] <- state_conditioned_step_length(observed_track, observed_states, simulated_tracks, simulated_states)
  }
  if ("state_conditioned_turning_angle" %in% diagnostics) {
    rows[[length(rows) + 1L]] <- state_conditioned_turning_angle(observed_track, observed_states, simulated_tracks, simulated_states)
  }
  if ("state_conditioned_msd" %in% diagnostics) {
    rows[[length(rows) + 1L]] <- state_conditioned_msd(observed_track, observed_states, simulated_tracks, simulated_states)
  }
  if ("state_conditioned_ud" %in% diagnostics) {
    rows[[length(rows) + 1L]] <- diagnostic_row("state_conditioned_ud", NA, NA_real_, numeric(), list(value = NA_real_), warning = "State-conditioned UD requires a spatial reference grid; not available in this run.")
  }
  if ("habitat_use_by_state" %in% diagnostics) {
    rows[[length(rows) + 1L]] <- diagnostic_row("habitat_use_by_state", NA, NA_real_, numeric(), list(value = NA_real_), warning = "Habitat-use diagnostics require state-specific covariates; not available in this run.")
  }
  do.call(rbind, rows)
}

global_ud_diagnostic <- function(observed_track, simulated_tracks) {
  obs_stat <- mean(vapply(simulated_tracks, function(x) track_distance_stat(observed_track, x), numeric(1L)), na.rm = TRUE)
  sim_stat <- leave_one_out_track_stats(simulated_tracks, track_distance_stat)
  diagnostic_row("ud_wasserstein", NA, obs_stat, sim_stat, compute_sharpness_scalar(sim_stat), "greater")
}

global_msd_diagnostic <- function(observed_track, simulated_tracks, max_lag = 50) {
  max_lag <- min(max_lag, nrow(observed_track) - 1L)
  obs_curve <- msd_curve(observed_track, max_lag = max_lag)$msd
  sim_curves <- do.call(cbind, lapply(simulated_tracks, function(x) msd_curve(x, max_lag = max_lag)$msd))
  sim_mean <- rowMeans(sim_curves, na.rm = TRUE)
  obs_stat <- sum((obs_curve - sim_mean)^2)
  sim_stat <- vapply(seq_len(ncol(sim_curves)), function(i) {
    loo <- rowMeans(sim_curves[, -i, drop = FALSE], na.rm = TRUE)
    sum((sim_curves[, i] - loo)^2)
  }, numeric(1L))
  diagnostic_row("msd", NA, obs_stat, sim_stat, compute_sharpness_curve(sim_curves), "greater")
}

global_sinuosity_diagnostic <- function(observed_track, simulated_tracks) {
  obs_value <- straightness_index(observed_track)
  sim_values <- vapply(simulated_tracks, straightness_index, numeric(1L))
  obs_stat <- abs(obs_value - mean(sim_values, na.rm = TRUE))
  sim_stat <- vapply(seq_along(sim_values), function(i) {
    abs(sim_values[i] - mean(sim_values[-i], na.rm = TRUE))
  }, numeric(1L))
  diagnostic_row("sinuosity", NA, obs_stat, sim_stat, compute_sharpness_scalar(sim_values), "greater")
}

track_distance_stat <- function(a, b) {
  a <- as_track_df(a)
  b <- as_track_df(b)
  wasserstein_1d(a$x, b$x) + wasserstein_1d(a$y, b$y)
}

leave_one_out_track_stats <- function(tracks, stat_fun) {
  vapply(seq_along(tracks), function(i) {
    others <- tracks[-i]
    mean(vapply(others, function(x) stat_fun(tracks[[i]], x), numeric(1L)), na.rm = TRUE)
  }, numeric(1L))
}

#' @export
print.gmov_hmmssf_diagnostics <- function(x, ...) {
  cat("gmov HMM-SSF diagnostics\n")
  cat("Methods:", paste(x$metadata$methods, collapse = ", "), "\n")
  cat("Simulations per method:", x$metadata$n_sims, "\n")
  print(x$diagnostics, row.names = FALSE)
  if (length(x$warnings) > 0L) {
    cat("\nWarnings\n")
    cat(paste0("- ", x$warnings, collapse = "\n"), "\n")
  }
  invisible(x)
}

#' @export
summary.gmov_hmmssf_diagnostics <- function(object, ...) {
  object$diagnostics
}
