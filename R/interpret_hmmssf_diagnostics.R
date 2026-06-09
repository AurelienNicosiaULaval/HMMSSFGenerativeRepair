#' Interpret HMM-SSF diagnostic rows
#'
#' @param diagnostics Diagnostic table.
#' @param alpha Compatibility threshold.
#'
#' @return Diagnostic table with interpretation labels.
#' @export
interpret_hmmssf_diagnostics <- function(diagnostics, alpha = 0.05) {
  if (!is.data.frame(diagnostics)) {
    stop("`diagnostics` must be a data frame.", call. = FALSE)
  }
  if (!"interpretation_label" %in% names(diagnostics)) {
    diagnostics$interpretation_label <- NA_character_
  }
  for (i in seq_len(nrow(diagnostics))) {
    diagnostics$interpretation_label[i] <- interpret_one_row(
      p_value = diagnostics$p_value[i],
      sharpness_label = diagnostics$sharpness_label[i],
      warning = diagnostics$warning[i],
      alpha = alpha
    )
  }
  diagnostics
}

interpret_one_row <- function(p_value, sharpness_label, warning = "", alpha = 0.05) {
  if (is.na(p_value) || identical(sharpness_label, "not_available")) {
    return("not_available")
  }
  if (!is.na(warning) && nzchar(warning) && grepl("low power|too few|few", warning, ignore.case = TRUE)) {
    return("low_power")
  }
  if (p_value < alpha) {
    return("fail")
  }
  if (sharpness_label == "sharp") {
    return("pass_sharp")
  }
  if (sharpness_label == "moderate") {
    return("pass_moderate")
  }
  if (sharpness_label %in% c("diffuse", "uninformative")) {
    return("pass_diffuse")
  }
  "ambiguous"
}

#' Compare simulation methods
#'
#' @param diagnostics Diagnostic table returned by [diagnose_hmmssf()].
#'
#' @return A data frame with method-comparison interpretation rules.
#' @export
compare_simulation_methods <- function(diagnostics) {
  if (!is.data.frame(diagnostics) || !all(c("method", "diagnostic", "interpretation_label") %in% names(diagnostics))) {
    stop("`diagnostics` must contain method, diagnostic, and interpretation_label columns.", call. = FALSE)
  }
  methods <- unique(diagnostics$method)
  rows <- list()
  diagnostics_names <- unique(diagnostics$diagnostic)
  for (diag_name in diagnostics_names) {
    subset_i <- diagnostics[diagnostics$diagnostic == diag_name & is.na(diagnostics$state), , drop = FALSE]
    label <- stats::setNames(subset_i$interpretation_label, subset_i$method)
    markov <- label["markov"]
    viterbi <- label["viterbi"]
    tube <- label["viterbi_tube"]
    posterior <- label["posterior"]
    message <- NA_character_
    if (!is.na(markov) && !is.na(viterbi) && startsWith(markov, "fail") && startsWith(viterbi, "pass")) {
      message <- "Likely failure source: latent transition dynamics, residence times, state occupancy, or switching structure."
    } else if (!is.na(markov) && !is.na(viterbi) && startsWith(markov, "fail") && startsWith(viterbi, "fail")) {
      message <- "Likely failure source: state-specific SSF kernels, missing covariates, movement distributions, state definition, or structural misspecification."
    } else if (!is.na(markov) && !is.na(viterbi) && startsWith(markov, "pass") && startsWith(viterbi, "fail")) {
      message <- "Possible compensation across states; global patterns may be reproduced for the wrong reason."
    } else if (!is.na(markov) && startsWith(markov, "pass_diffuse")) {
      message <- "Model is compatible but too variable; do not claim strong generative realism."
    }
    if (!is.na(tube) && !is.na(markov) && startsWith(markov, "pass") && tube %in% c("fail", "ambiguous")) {
      message <- "Latent-state uncertainty affects validation conclusions."
    }
    if (!is.na(tube) && !is.na(posterior) && tube != posterior) {
      message <- paste(c(message, "Posterior and Viterbi-tube diagnostics disagree; flag latent-path sensitivity."), collapse = " ")
    }
    rows[[length(rows) + 1L]] <- data.frame(diagnostic = diag_name, methods = paste(methods, collapse = ","), comparison_message = message)
  }
  do.call(rbind, rows)
}
