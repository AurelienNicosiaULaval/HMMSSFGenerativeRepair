#' Simulate one candidate step from a state-specific SSF kernel
#'
#' @param current_location Numeric length-2 vector.
#' @param previous_location Numeric length-2 vector.
#' @param state Integer state.
#' @param kernel State-specific kernel definition.
#' @param habitat_function Optional function returning log selection weights or
#'   covariates for candidate endpoints.
#' @param n_available Number of available candidate endpoints.
#' @param domain Optional function returning `TRUE` for valid endpoints.
#' @param max_attempts Maximum attempts to find valid candidates.
#' @param ... Reserved for future extensions.
#'
#' @return Numeric length-2 next location.
#' @export
simulate_step_from_state <- function(
  current_location,
  previous_location,
  state,
  kernel,
  habitat_function = NULL,
  n_available = 100,
  domain = NULL,
  max_attempts = 1000,
  ...
) {
  current_location <- as.numeric(current_location)
  previous_location <- as.numeric(previous_location)
  if (length(current_location) != 2L || length(previous_location) != 2L) {
    stop("Locations must be numeric length-2 vectors.", call. = FALSE)
  }
  previous_bearing <- atan2(current_location[2] - previous_location[2], current_location[1] - previous_location[1])
  if (!is.finite(previous_bearing)) {
    previous_bearing <- 0
  }

  for (attempt in seq_len(max_attempts)) {
    steps <- draw_step_lengths(kernel, n_available)
    turns <- draw_turning_angles(kernel, n_available)
    bearings <- wrap_angle(previous_bearing + turns)
    candidates <- data.frame(
      x = current_location[1] + steps * cos(bearings),
      y = current_location[2] + steps * sin(bearings),
      step = steps,
      angle = turns
    )

    valid <- is.finite(candidates$x) & is.finite(candidates$y)
    if (!is.null(domain)) {
      valid <- valid & as.logical(domain(candidates[, c("x", "y")]))
    }
    candidates <- candidates[valid, , drop = FALSE]
    if (nrow(candidates) == 0L) {
      next
    }

    log_weights <- rep(0, nrow(candidates))
    if (!is.null(habitat_function)) {
      habitat_out <- habitat_function(candidates[, c("x", "y")], state = state, kernel = kernel)
      if (is.numeric(habitat_out) && length(habitat_out) == nrow(candidates)) {
        log_weights <- log_weights + habitat_out
      } else if (is.data.frame(habitat_out) && !is.null(kernel$selection)) {
        mm <- stats::model.matrix(kernel$selection$formula, habitat_out)
        beta <- kernel$selection$coef
        log_weights <- log_weights + as.numeric(mm[, names(beta), drop = FALSE] %*% beta)
      } else {
        stop("`habitat_function` must return log weights or covariates compatible with `kernel$selection`.", call. = FALSE)
      }
    }

    idx <- sample_log_weights(log_weights)
    return(as.numeric(candidates[idx, c("x", "y")]))
  }

  stop("Could not simulate a valid step after `max_attempts` attempts.", call. = FALSE)
}

#' Simulate complete trajectories from a fitted HMM-SSF-like model
#'
#' @param fit Fitted object or generic list object.
#' @param observed_track Observed track used for time horizon and starting point.
#' @param n_sims Number of simulated trajectories.
#' @param method Simulation method.
#' @param epsilon Viterbi tube tolerance.
#' @param n_state_paths Optional number of state paths.
#' @param parameter_uncertainty Should parameter uncertainty be propagated?
#' @param state_uncertainty Should latent-state uncertainty be represented when
#'   the chosen method supports it?
#' @param start Start rule.
#' @param max_attempts Maximum simulation attempts.
#' @param seed Optional seed.
#' @param adapter Adapter choice.
#' @param ... Reserved for future extensions.
#'
#' @return Object of class `"gmov_hmmssf_simulations"`.
#' @export
simulate_generative_hmmssf <- function(
  fit,
  observed_track,
  n_sims = 99,
  method = c("markov", "viterbi", "viterbi_tube", "posterior"),
  epsilon = NULL,
  n_state_paths = NULL,
  parameter_uncertainty = TRUE,
  state_uncertainty = TRUE,
  start = c("observed", "stationary"),
  max_attempts = 1000,
  seed = NULL,
  adapter = c("auto", "hmmSSF", "list"),
  ...
) {
  method <- match.arg(method)
  start <- match.arg(start)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  observed_track <- as_track_df(observed_track)
  components <- as_gmov_hmmssf(fit, observed_track = observed_track, adapter = adapter)
  n_steps <- nrow(observed_track)
  n_sims <- as.integer(n_sims)
  n_state_paths <- as.integer(n_state_paths %||% n_sims)
  warnings <- character()
  parameter_draws <- NULL

  if (parameter_uncertainty) {
    if (is.null(components$covariance)) {
      warnings <- c(warnings, "Parameter uncertainty was requested but no covariance matrix was available; fixed fitted parameters were used.")
    } else {
      warnings <- c(warnings, "Covariance-based parameter draws are recorded as planned, but kernel-specific parameter remapping is not implemented for this adapter; fixed fitted parameters were used.")
    }
  }

  state_result <- simulate_state_paths_for_method(
    components = components,
    method = method,
    n_steps = n_steps,
    n_sims = n_state_paths,
    epsilon = epsilon,
    start = start,
    seed = seed,
    max_attempts = max_attempts
  )
  state_paths <- state_result$states
  warnings <- c(warnings, state_result$warnings)
  if (nrow(state_paths) < n_sims) {
    state_paths <- state_paths[rep(seq_len(nrow(state_paths)), length.out = n_sims), , drop = FALSE]
  } else {
    state_paths <- state_paths[seq_len(n_sims), , drop = FALSE]
  }

  simulated_tracks <- vector("list", n_sims)
  for (i in seq_len(n_sims)) {
    simulated_tracks[[i]] <- simulate_track_from_states(
      states = state_paths[i, ],
      kernels = components$kernels,
      observed_track = observed_track,
      habitat_function = components$habitat_function,
      domain = components$domain,
      start = start,
      max_attempts = max_attempts
    )
  }
  names(simulated_tracks) <- paste0("sim_", seq_len(n_sims))

  if (identical(method, "viterbi")) {
    warnings <- c(
      warnings,
      "Viterbi-conditioned simulations are diagnostic only and do not treat Viterbi states as true states."
    )
  }

  out <- list(
    observed_track = observed_track,
    simulated_tracks = simulated_tracks,
    simulated_states = state_paths,
    method = method,
    model_parameters = components$parameters,
    parameter_draws = parameter_draws,
    state_paths = state_result,
    metadata = list(n_sims = n_sims, n_steps = n_steps, start = start, state_uncertainty = state_uncertainty),
    warnings = unique(warnings),
    seed = seed
  )
  class(out) <- c("gmov_hmmssf_simulations", "list")
  out
}

#' @export
print.gmov_hmmssf_simulations <- function(x, ...) {
  cat("gmov HMM-SSF simulations\n")
  cat("Method:", x$method, "\n")
  cat("Observed locations:", nrow(x$observed_track), "\n")
  cat("Simulated tracks:", length(x$simulated_tracks), "\n")
  if (length(x$warnings) > 0L) {
    cat("Warnings:\n")
    cat(paste0("- ", x$warnings, collapse = "\n"), "\n")
  }
  invisible(x)
}

simulate_state_paths_for_method <- function(components, method, n_steps, n_sims, epsilon, start, seed, max_attempts) {
  warnings <- character()
  if (method == "markov") {
    initial <- if (start == "stationary" && is.matrix(components$transition)) {
      stationary_distribution(components$transition)
    } else {
      components$initial
    }
    states <- simulate_markov_states(initial, components$transition, n_steps = n_steps, n_sims = n_sims, seed = seed)
    return(list(states = states, method = method, warnings = warnings))
  }
  if (method == "viterbi") {
    states <- simulate_viterbi_states(components, n_sims = n_sims)
    return(list(states = states, method = method, warnings = warnings))
  }
  if (method == "viterbi_tube") {
    tube <- simulate_viterbi_tube_states(components, epsilon = epsilon, n_sims = n_sims, seed = seed, max_attempts = max_attempts)
    return(list(states = tube$states, method = method, tube = tube, warnings = tube$warnings))
  }
  if (method == "posterior") {
    states <- simulate_posterior_states(components, n_sims = n_sims, seed = seed)
    return(list(states = states, method = method, warnings = warnings))
  }
  stop("Unknown simulation method.", call. = FALSE)
}

simulate_track_from_states <- function(
  states,
  kernels,
  observed_track = NULL,
  start_location = NULL,
  habitat_function = NULL,
  domain = NULL,
  start = "observed",
  max_attempts = 1000,
  seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  states <- as.integer(states)
  n_steps <- length(states)
  if (is.null(start_location)) {
    if (!is.null(observed_track) && start == "observed") {
      observed_track <- as_track_df(observed_track)
      start_location <- as.numeric(observed_track[1L, c("x", "y")])
      previous_location <- if (nrow(observed_track) >= 2L) {
        as.numeric(observed_track[1L, c("x", "y")]) - (as.numeric(observed_track[2L, c("x", "y")]) - as.numeric(observed_track[1L, c("x", "y")]))
      } else {
        start_location - c(1, 0)
      }
    } else {
      start_location <- c(0, 0)
      previous_location <- c(-1, 0)
    }
  } else {
    previous_location <- start_location - c(1, 0)
  }

  track <- matrix(NA_real_, nrow = n_steps, ncol = 2L)
  track[1L, ] <- start_location
  if (n_steps >= 2L) {
    for (t in seq_len(n_steps - 1L)) {
      kernel <- kernels[[states[t + 1L]]]
      next_location <- simulate_step_from_state(
        current_location = track[t, ],
        previous_location = previous_location,
        state = states[t + 1L],
        kernel = kernel,
        habitat_function = habitat_function,
        domain = domain,
        max_attempts = max_attempts
      )
      previous_location <- track[t, ]
      track[t + 1L, ] <- next_location
    }
  }
  data.frame(x = track[, 1L], y = track[, 2L])
}

draw_step_lengths <- function(kernel, n) {
  step <- kernel$step %||% kernel$step_length
  if (is.function(step)) {
    return(step(n))
  }
  dist <- step$dist %||% "gamma"
  if (dist == "gamma") {
    return(stats::rgamma(n, shape = step$shape %||% 2, scale = step$scale %||% 1))
  }
  if (dist == "lognormal") {
    return(stats::rlnorm(n, meanlog = step$meanlog %||% 0, sdlog = step$sdlog %||% 1))
  }
  if (dist == "exponential") {
    return(stats::rexp(n, rate = step$rate %||% 1))
  }
  if (dist == "fixed") {
    return(rep(step$value %||% 1, n))
  }
  stop("Unsupported step-length distribution: ", dist, call. = FALSE)
}

draw_turning_angles <- function(kernel, n) {
  angle <- kernel$angle %||% kernel$turning_angle
  if (is.function(angle)) {
    return(wrap_angle(angle(n)))
  }
  dist <- angle$dist %||% "uniform"
  if (dist == "uniform") {
    return(stats::runif(n, -pi, pi))
  }
  if (dist == "wrapped_normal") {
    return(wrap_angle(stats::rnorm(n, mean = angle$mean %||% 0, sd = angle$sd %||% 1)))
  }
  if (dist == "fixed") {
    return(rep(angle$value %||% 0, n))
  }
  stop("Unsupported turning-angle distribution: ", dist, call. = FALSE)
}

stationary_distribution <- function(transition) {
  transition <- check_transition_matrix(transition)
  stat <- Re(eigen(t(transition))$vectors[, 1L])
  stat <- abs(stat)
  check_probability_vector(stat)
}
