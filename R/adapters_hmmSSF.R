#' Adapt fitted HMM-SSF objects to the gmov HMM-SSF interface
#'
#' `as_gmov_hmmssf()` converts either a generic list object or, when possible, a
#' fitted `hmmSSF` object into a defensive internal representation. The adapter
#' does not rewrite `hmmSSF` internals. If required quantities are unavailable,
#' it returns clear errors rather than guessing.
#'
#' @param fit A fitted object or generic list object.
#' @param observed_track Optional observed track.
#' @param adapter Adapter choice.
#' @param ... Reserved for future adapter-specific arguments.
#'
#' @return A list with standardized HMM-SSF components.
#' @export
as_gmov_hmmssf <- function(fit, observed_track = NULL, adapter = c("auto", "hmmSSF", "list"), ...) {
  adapter <- match.arg(adapter)
  if (adapter == "auto") {
    adapter <- if (inherits(fit, "hmmSSF")) "hmmSSF" else "list"
  }

  if (adapter == "list") {
    out <- as_gmov_hmmssf_list(fit, observed_track = observed_track)
  } else {
    out <- as_gmov_hmmssf_hmmSSF(fit, observed_track = observed_track, ...)
  }

  validate_hmmssf_adapter(out)
  out
}

as_gmov_hmmssf_list <- function(fit, observed_track = NULL) {
  if (!is.list(fit)) {
    stop("The generic list adapter requires `fit` to be a list.", call. = FALSE)
  }

  initial <- fit$initial %||% fit$delta
  transition <- fit$transition %||% fit$gamma
  kernels <- fit$kernels
  if (is.null(initial) || is.null(transition) || is.null(kernels)) {
    stop(
      "A generic HMM-SSF list must provide `initial`, `transition`, and `kernels`.",
      call. = FALSE
    )
  }

  initial <- check_probability_vector(initial, "initial")
  if (is.matrix(transition)) {
    transition <- check_transition_matrix(transition, "transition")
  } else if (is.array(transition) && length(dim(transition)) == 3L) {
    for (i in seq_len(dim(transition)[3L])) {
      transition[, , i] <- check_transition_matrix(transition[, , i], "transition array slice")
    }
  } else if (!is.function(transition)) {
    stop("`transition` must be a matrix, 3D array, or function.", call. = FALSE)
  }

  observed_track <- observed_track %||% fit$observed_track
  if (!is.null(observed_track)) {
    observed_track <- as_track_df(observed_track)
  }

  n_states <- length(initial)
  if (length(kernels) != n_states) {
    stop("`kernels` must contain one state-specific kernel per state.", call. = FALSE)
  }

  log_emission <- fit$log_emission %||% fit$emission_loglik
  if (!is.null(log_emission)) {
    log_emission <- as.matrix(log_emission)
  }

  log_gamma <- fit$log_gamma
  if (is.null(log_gamma) && !is.function(transition)) {
    n_steps <- if (!is.null(log_emission)) nrow(log_emission) else if (!is.null(observed_track)) nrow(observed_track) else NULL
    if (!is.null(n_steps)) {
      log_gamma <- transition_to_log_array(transition, n_steps)
    }
  }

  structure(
    list(
      model = "list",
      n_states = n_states,
      initial = initial,
      transition = transition,
      log_delta = safe_log(initial),
      log_gamma = log_gamma,
      log_emission = log_emission,
      kernels = kernels,
      habitat_function = fit$habitat_function,
      domain = fit$domain,
      viterbi_path = fit$viterbi_path,
      posterior = fit$posterior,
      covariance = fit$covariance,
      parameters = fit$parameters,
      observed_track = observed_track
    ),
    class = c("gmov_hmmssf_adapter", "list")
  )
}

as_gmov_hmmssf_hmmSSF <- function(fit, observed_track = NULL, ...) {
  if (!inherits(fit, "hmmSSF")) {
    stop("The `hmmSSF` adapter requires an object inheriting from `hmmSSF`.", call. = FALSE)
  }
  if (!requireNamespace("hmmSSF", quietly = TRUE)) {
    stop("Package `hmmSSF` is required to adapt `hmmSSF` objects.", call. = FALSE)
  }

  data <- fit$args$data
  observed_track <- observed_track %||% hmmssf_observed_track(data)
  initial <- extract_initial_distribution(fit)
  transition <- extract_transition_probabilities(fit)
  log_emission <- extract_state_likelihoods(fit, required = FALSE)
  log_gamma <- if (!is.null(transition)) transition_to_log_array(transition, nrow(log_emission %||% observed_track)) else NULL

  structure(
    list(
      model = "hmmSSF",
      n_states = fit$args$n_states,
      initial = initial,
      transition = transition,
      log_delta = safe_log(initial),
      log_gamma = log_gamma,
      log_emission = log_emission,
      kernels = extract_state_specific_kernels(fit, required = FALSE),
      habitat_function = NULL,
      domain = NULL,
      viterbi_path = extract_viterbi_path(fit, required = FALSE),
      posterior = try_extract_local_decoding(fit),
      covariance = try_extract_covariance(fit),
      parameters = extract_hmmssf_parameters(fit),
      observed_track = observed_track
    ),
    class = c("gmov_hmmssf_adapter", "list")
  )
}

#' Validate an adapted HMM-SSF object
#'
#' @param x Adapted object returned by [as_gmov_hmmssf()].
#'
#' @return Invisibly returns `x`.
#' @export
validate_hmmssf_adapter <- function(x) {
  required <- c("n_states", "initial", "transition", "kernels")
  missing <- required[!vapply(required, function(name) !is.null(x[[name]]), logical(1L))]
  if (length(missing) > 0L) {
    stop("Adapted object is missing required element(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  check_probability_vector(x$initial, "initial")
  if (!is.function(x$transition)) {
    if (is.matrix(x$transition)) {
      check_transition_matrix(x$transition, "transition")
    } else if (!is.array(x$transition) || length(dim(x$transition)) != 3L) {
      stop("`transition` must be a matrix, 3D array, or function.", call. = FALSE)
    }
  }
  if (length(x$kernels) != x$n_states) {
    stop("`kernels` must contain one kernel per state.", call. = FALSE)
  }
  invisible(x)
}

#' Extract HMM-SSF parameters
#'
#' @param fit A fitted object.
#'
#' @return A list of available parameters.
#' @export
extract_hmmssf_parameters <- function(fit) {
  if (inherits(fit, "hmmSSF")) {
    return(fit$par)
  }
  if (is.list(fit)) {
    return(fit$parameters %||% list(initial = fit$initial, transition = fit$transition))
  }
  stop("No parameter extractor is available for this object.", call. = FALSE)
}

#' Extract transition probabilities
#'
#' @param fit A fitted object.
#'
#' @return A transition matrix, transition array, or transition function.
#' @export
extract_transition_probabilities <- function(fit) {
  if (inherits(fit, "hmmSSF")) {
    if (!requireNamespace("hmmSSF", quietly = TRUE)) {
      stop("Package `hmmSSF` is required to extract transition probabilities.", call. = FALSE)
    }
    gamma <- hmmSSF::predict_tpm(fit)
    if (is.list(gamma) && !is.null(gamma$est)) {
      gamma <- gamma$est
    }
    return(gamma)
  }
  if (is.list(fit)) {
    return(fit$transition %||% fit$gamma)
  }
  stop("No transition-probability extractor is available for this object.", call. = FALSE)
}

#' Extract the initial distribution
#'
#' @param fit A fitted object.
#'
#' @return A probability vector.
#' @export
extract_initial_distribution <- function(fit) {
  if (inherits(fit, "hmmSSF")) {
    gamma <- extract_transition_probabilities(fit)
    gamma1 <- if (is.array(gamma)) gamma[, , 1L] else gamma
    return(check_probability_vector(solve(t(diag(nrow(gamma1)) - gamma1 + 1), rep(1, nrow(gamma1)))))
  }
  if (is.list(fit)) {
    return(check_probability_vector(fit$initial %||% fit$delta, "initial"))
  }
  stop("No initial-distribution extractor is available for this object.", call. = FALSE)
}

#' Extract Viterbi path
#'
#' @param fit A fitted object.
#' @param required Should missing paths trigger an error?
#'
#' @return Integer state path or `NULL`.
#' @export
extract_viterbi_path <- function(fit, required = TRUE) {
  out <- NULL
  if (inherits(fit, "hmmSSF") && requireNamespace("hmmSSF", quietly = TRUE)) {
    out <- tryCatch(hmmSSF::viterbi_decoding(fit), error = function(e) NULL)
  } else if (is.list(fit)) {
    out <- fit$viterbi_path
  }
  if (is.null(out) && required) {
    stop("A Viterbi path is unavailable. Provide `viterbi_path` or an object exposing Viterbi decoding.", call. = FALSE)
  }
  if (is.null(out)) NULL else as.integer(out)
}

#' Extract state likelihoods
#'
#' @param fit A fitted object.
#' @param required Should missing likelihoods trigger an error?
#'
#' @return Matrix of log likelihoods or `NULL`.
#' @export
extract_state_likelihoods <- function(fit, required = TRUE) {
  out <- NULL
  if (is.list(fit) && !inherits(fit, "hmmSSF")) {
    out <- fit$log_emission %||% fit$emission_loglik
  }
  if (inherits(fit, "hmmSSF")) {
    out <- tryCatch(extract_hmmSSF_log_emission(fit), error = function(e) NULL)
  }
  if (is.null(out) && required) {
    stop(
      "State log-likelihoods are unavailable. Provide `log_emission` or an object whose state likelihoods can be computed.",
      call. = FALSE
    )
  }
  if (is.null(out)) NULL else as.matrix(out)
}

#' Extract state-specific movement kernels
#'
#' @param fit A fitted object.
#' @param required Should missing kernels trigger an error?
#'
#' @return List of kernels or `NULL`.
#' @export
extract_state_specific_kernels <- function(fit, required = TRUE) {
  out <- if (is.list(fit) && !inherits(fit, "hmmSSF")) fit$kernels else NULL
  if (is.null(out) && required) {
    stop(
      "State-specific kernels are unavailable. Provide a `kernels` list for generic simulation.",
      call. = FALSE
    )
  }
  out
}

transition_to_log_array <- function(transition, n_steps) {
  if (is.null(n_steps) || !is.finite(n_steps) || n_steps < 1L) {
    stop("`n_steps` must be available to build log transition arrays.", call. = FALSE)
  }
  if (is.matrix(transition)) {
    transition <- check_transition_matrix(transition)
    arr <- array(NA_real_, dim = c(nrow(transition), ncol(transition), max(1L, n_steps - 1L)))
    for (i in seq_len(dim(arr)[3L])) {
      arr[, , i] <- safe_log(transition)
    }
    return(arr)
  }
  if (is.array(transition) && length(dim(transition)) == 3L) {
    return(safe_log(transition))
  }
  NULL
}

hmmssf_observed_track <- function(data) {
  if (is.null(data) || !is.data.frame(data) || !"obs" %in% names(data)) {
    stop("Cannot extract observed track from this `hmmSSF` object.", call. = FALSE)
  }
  obs <- data[data$obs == 1, , drop = FALSE]
  if (all(c("x", "y") %in% names(obs))) {
    return(as_track_df(obs[, c("x", "y")]))
  }
  if (all(c("x_", "y_") %in% names(obs))) {
    return(as_track_df(obs[, c("x_", "y_")]))
  }
  stop("Observed hmmSSF data must contain coordinates.", call. = FALSE)
}

try_extract_local_decoding <- function(fit) {
  if (!inherits(fit, "hmmSSF") || !requireNamespace("hmmSSF", quietly = TRUE)) {
    return(NULL)
  }
  tryCatch(hmmSSF::local_decoding(fit), error = function(e) NULL)
}

try_extract_covariance <- function(fit) {
  hessian <- fit$fit$hessian
  if (is.null(hessian)) {
    return(NULL)
  }
  tryCatch(solve(hessian), error = function(e) NULL)
}

extract_hmmSSF_log_emission <- function(fit) {
  if (!requireNamespace("hmmSSF", quietly = TRUE)) {
    stop("Package `hmmSSF` is required.", call. = FALSE)
  }
  data <- fit$args$data
  obs <- data[data$obs == 1, , drop = FALSE]
  if (is.null(obs) || nrow(obs) == 0L) {
    stop("No observed rows are available in the hmmSSF model data.", call. = FALSE)
  }
  ssf_mm <- stats::model.matrix(fit$args$ssf_formula, data)
  ssf_mm <- ssf_mm[, colnames(ssf_mm) != "(Intercept)", drop = FALSE]
  linear_pred <- ssf_mm %*% as.matrix(fit$par$ssf)
  state_dens <- utils::getFromNamespace("state_dens_rcpp", "hmmSSF")
  dens <- state_dens(
    linear_pred = linear_pred,
    stratum = data$stratum,
    n_states = fit$args$n_states,
    sampling_densities = data$w,
    n_obs = nrow(obs)
  )
  safe_log(dens)
}

#' Create a small simulated generic HMM-SSF-like object
#'
#' This object is for tests, examples, and vignettes only. It is not empirical
#' data.
#'
#' @param n Number of observed locations.
#' @param seed Optional random seed.
#'
#' @return A generic list object accepted by [as_gmov_hmmssf()].
#' @export
example_hmmssf_list <- function(n = 80, seed = 1) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  transition <- matrix(c(0.92, 0.08, 0.12, 0.88), nrow = 2, byrow = TRUE)
  initial <- c(0.6, 0.4)
  kernels <- list(
    list(step = list(dist = "gamma", shape = 2.5, scale = 0.35), angle = list(dist = "wrapped_normal", mean = 0, sd = 1.8)),
    list(step = list(dist = "gamma", shape = 3.0, scale = 1.1), angle = list(dist = "wrapped_normal", mean = 0, sd = 0.35))
  )

  states <- simulate_markov_states(initial, transition, n_steps = n, n_sims = 1, seed = seed)[1, ]
  track <- simulate_track_from_states(states, kernels, start_location = c(0, 0), seed = seed)
  log_emission <- matrix(0, nrow = n, ncol = 2)
  log_emission[cbind(seq_len(n), states)] <- 1
  posterior <- exp(log_emission)
  posterior <- posterior / rowSums(posterior)

  list(
    initial = initial,
    transition = transition,
    kernels = kernels,
    observed_track = track,
    viterbi_path = states,
    posterior = posterior,
    log_emission = log_emission,
    parameters = list(note = "simulated toy object")
  )
}
