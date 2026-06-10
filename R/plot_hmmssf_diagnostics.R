#' Plot HMM-SSF diagnostics
#'
#' @param x Object returned by [diagnose_hmmssf()].
#' @param type Plot type.
#' @param method Simulation method to display for method-specific dashboards.
#' @param ... Currently unused.
#'
#' @return A ggplot object or list of ggplot objects.
#' @export
plot.gmov_hmmssf_diagnostics <- function(
  x,
  type = c("dashboard", "p_values", "msd", "state_occupancy", "residence_time", "method_comparison", "tracks"),
  method = NULL,
  ...
) {
  type <- match.arg(type)
  if (type == "dashboard") {
    return(plot_full_hmmssf_dashboard(x, method = method))
  }
  if (type == "p_values") {
    return(plot_diagnostic_p_values(x))
  }
  if (type == "msd") {
    return(plot_msd_dashboard(x))
  }
  if (type == "state_occupancy") {
    return(plot_state_occupancy_dashboard(x))
  }
  if (type == "residence_time") {
    return(plot_residence_time_dashboard(x))
  }
  if (type == "tracks") {
    return(plot_tracks_dashboard(x, method = method))
  }
  plot_method_comparison(x)
}

plot_diagnostic_p_values <- function(x) {
  df <- x$diagnostics
  df <- df[is.na(df$state), , drop = FALSE]
  ggplot2::ggplot(df, ggplot2::aes(x = diagnostic, y = p_value, fill = interpretation_label)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(~method) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Diagnostic", y = "Monte Carlo p-value", fill = "Interpretation")
}

plot_full_hmmssf_dashboard <- function(x, method = NULL) {
  method <- choose_simulation_method(x, method)
  sim <- x$simulations[[method]]
  plots <- list(
    trajectories = plot_tracks_panel(sim, method = method),
    ud = plot_ud_panel(sim, x$diagnostics, method),
    msd = plot_msd_panel(sim, x$diagnostics, method),
    sinuosity = plot_sinuosity_panel(sim, x$diagnostics, method),
    state_occupancy = plot_state_occupancy_panel(x$diagnostics, method),
    state_diagnostics = plot_state_diagnostic_panel(x$diagnostics, method),
    sharpness = plot_sharpness_panel(x$diagnostics, method),
    interpretation = plot_interpretation_panel(x$diagnostics, method)
  )
  draw_hmmssf_dashboard(plots)
  invisible(plots)
}

plot_tracks_dashboard <- function(x, method = NULL) {
  method <- choose_simulation_method(x, method)
  plot_tracks_panel(x$simulations[[method]], method = method)
}

plot_tracks_panel <- function(sim, method = NULL) {
  simulated_tracks <- sim$simulated_tracks
  track_df <- do.call(
    rbind,
    Map(
      function(track, id) data.frame(sim_id = id, x = track$x, y = track$y),
      simulated_tracks,
      names(simulated_tracks)
    )
  )
  obs <- sim$observed_track
  text <- simulation_panel_text(sim, method = method)
  ggplot2::ggplot() +
    ggplot2::geom_path(
      data = track_df,
      ggplot2::aes(x = x, y = y, group = sim_id),
      color = "#4C78A8",
      alpha = 0.18,
      linewidth = 0.25
    ) +
    ggplot2::geom_path(
      data = obs,
      ggplot2::aes(x = x, y = y),
      color = "#D55E00",
      linewidth = 0.9
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(
      x = "x-coordinate",
      y = "y-coordinate",
      title = paste0("Observed and simulated trajectories: ", text$label),
      subtitle = text$subtitle,
      caption = text$caption
    )
}

plot_ud_panel <- function(sim, diagnostics, method) {
  sim_stats <- leave_one_out_track_stats(sim$simulated_tracks, track_distance_stat)
  observed <- mean(vapply(sim$simulated_tracks, function(track) {
    track_distance_stat(sim$observed_track, track)
  }, numeric(1L)), na.rm = TRUE)
  row <- get_diagnostic_row(diagnostics, method, "ud_wasserstein")
  ggplot2::ggplot(data.frame(statistic = sim_stats), ggplot2::aes(x = statistic)) +
    ggplot2::geom_histogram(bins = 18, fill = "#4C78A8", color = "white", alpha = 0.75) +
    ggplot2::geom_vline(xintercept = observed, color = "#D55E00", linewidth = 1) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(
      x = "UD discrepancy",
      y = "Number of simulations",
      title = "Emergent utilization distribution",
      subtitle = diagnostic_subtitle(row)
    )
}

plot_msd_panel <- function(sim, diagnostics, method, max_lag = 50) {
  max_lag <- min(max_lag, nrow(sim$observed_track) - 1L)
  obs <- msd_curve(sim$observed_track, max_lag = max_lag)
  sim_curves <- do.call(cbind, lapply(sim$simulated_tracks, function(track) {
    msd_curve(track, max_lag = max_lag)$msd
  }))
  df <- data.frame(
    lag = obs$lag,
    observed = obs$msd,
    median = apply(sim_curves, 1L, stats::median, na.rm = TRUE),
    lo = apply(sim_curves, 1L, stats::quantile, probs = 0.025, na.rm = TRUE, type = 8),
    hi = apply(sim_curves, 1L, stats::quantile, probs = 0.975, na.rm = TRUE, type = 8)
  )
  row <- get_diagnostic_row(diagnostics, method, "msd")
  ggplot2::ggplot(df, ggplot2::aes(x = lag)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "#4C78A8", alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(y = median), color = "#4C78A8", linetype = "dashed") +
    ggplot2::geom_line(ggplot2::aes(y = observed), color = "#D55E00", linewidth = 0.9) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(
      x = "Lag",
      y = "Mean squared displacement",
      title = "Mean squared displacement",
      subtitle = diagnostic_subtitle(row)
    )
}

plot_sinuosity_panel <- function(sim, diagnostics, method) {
  sim_values <- vapply(sim$simulated_tracks, straightness_index, numeric(1L))
  observed <- straightness_index(sim$observed_track)
  row <- get_diagnostic_row(diagnostics, method, "sinuosity")
  ggplot2::ggplot(data.frame(straightness = sim_values), ggplot2::aes(x = straightness)) +
    ggplot2::geom_histogram(bins = 18, fill = "#4C78A8", color = "white", alpha = 0.75) +
    ggplot2::geom_vline(xintercept = observed, color = "#D55E00", linewidth = 1) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(
      x = "Straightness index",
      y = "Number of simulations",
      title = "Path sinuosity",
      subtitle = diagnostic_subtitle(row)
    )
}

plot_state_occupancy_panel <- function(diagnostics, method) {
  df <- diagnostics[
    diagnostics$method == method &
      diagnostics$diagnostic == "state_occupancy" &
      !is.na(diagnostics$state),
    ,
    drop = FALSE
  ]
  if (nrow(df) == 0L) {
    return(empty_panel("State occupancy", "Not available"))
  }
  plot_df <- rbind(
    data.frame(state = df$state, source = "observed discrepancy", value = df$observed_statistic),
    data.frame(state = df$state, source = "simulated median", value = df$simulated_median)
  )
  ggplot2::ggplot(plot_df, ggplot2::aes(x = factor(state), y = value, fill = source)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(
      x = "State",
      y = "Occupancy discrepancy",
      title = "State occupancy",
      fill = NULL
    )
}

plot_state_diagnostic_panel <- function(diagnostics, method) {
  keep <- diagnostics$method == method &
    diagnostics$diagnostic %in% c(
      "state_residence_time",
      "switching_rate",
      "transition_counts",
      "state_conditioned_step_length",
      "state_conditioned_turning_angle",
      "state_conditioned_msd"
    )
  df <- diagnostics[keep, , drop = FALSE]
  if (nrow(df) == 0L) {
    return(empty_panel("State diagnostics", "Not available"))
  }
  df$label <- ifelse(is.na(df$state), df$diagnostic, paste0(df$diagnostic, " S", df$state))
  df$label <- factor(df$label, levels = rev(df$label))
  ggplot2::ggplot(df, ggplot2::aes(x = p_value, y = label, color = interpretation_label)) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_vline(xintercept = 0.05, linetype = "dashed", color = "grey50") +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::labs(
      x = "Monte Carlo p-value",
      y = NULL,
      title = "State-specific diagnostics",
      color = "Interpretation"
    )
}

plot_sharpness_panel <- function(diagnostics, method) {
  df <- diagnostics[diagnostics$method == method & is.na(diagnostics$state), , drop = FALSE]
  if (nrow(df) == 0L) {
    return(empty_panel("Sharpness", "Not available"))
  }
  df$diagnostic <- factor(df$diagnostic, levels = rev(df$diagnostic))
  ggplot2::ggplot(df, ggplot2::aes(x = sharpness_value, y = diagnostic, fill = sharpness_label)) +
    ggplot2::geom_col() +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::labs(
      x = "Sharpness value",
      y = NULL,
      title = "Sharpness of global diagnostics",
      fill = "Sharpness"
    )
}

plot_interpretation_panel <- function(diagnostics, method) {
  df <- diagnostics[diagnostics$method == method & is.na(diagnostics$state), , drop = FALSE]
  if (nrow(df) == 0L) {
    return(empty_panel("Interpretation", "Not available"))
  }
  df$diagnostic <- factor(df$diagnostic, levels = rev(df$diagnostic))
  ggplot2::ggplot(df, ggplot2::aes(x = method, y = diagnostic, fill = interpretation_label)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = "Compatibility and sharpness labels",
      fill = "Interpretation"
    )
}

draw_hmmssf_dashboard <- function(plots) {
  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = 3,
    ncol = 3,
    widths = grid::unit(c(1.45, 1, 1), "null"),
    heights = grid::unit(c(1, 1, 1), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  print(plots$trajectories, vp = grid::viewport(layout.pos.row = 1:3, layout.pos.col = 1))
  print(plots$ud, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(plots$msd, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 3))
  print(plots$sinuosity, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
  print(plots$state_occupancy, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 3))
  print(plots$state_diagnostics, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 2))
  print(plots$sharpness, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 3))
  grid::popViewport()
  invisible(plots)
}

choose_simulation_method <- function(x, method = NULL) {
  available <- names(x$simulations)
  if (length(available) == 0L) {
    stop("No simulations are available in this diagnostics object.", call. = FALSE)
  }
  if (is.null(method)) {
    return(available[1L])
  }
  if (!method %in% available) {
    stop(
      "`method` must be one of: ",
      paste(available, collapse = ", "),
      call. = FALSE
    )
  }
  method
}

simulation_panel_text <- function(sim, method = NULL) {
  metadata <- sim$metadata %||% list()
  label <- metadata$simulation_label %||% metadata$method_label %||% sim$method %||% method %||% "simulation"
  description <- metadata$simulation_description %||% ""
  detail <- metadata$simulation_detail %||% ""
  n_sims <- length(sim$simulated_tracks)
  subtitle <- paste0(n_sims, " simulated trajectories")
  if (nzchar(description)) {
    subtitle <- paste0(subtitle, "\n", paste(strwrap(description, width = 74), collapse = "\n"))
  }
  caption <- ""
  if (nzchar(detail)) {
    caption <- paste(strwrap(detail, width = 96), collapse = "\n")
  }
  list(label = label, subtitle = subtitle, caption = caption)
}

get_diagnostic_row <- function(diagnostics, method, diagnostic_name) {
  row <- diagnostics[
    diagnostics$method == method &
      diagnostics$diagnostic == diagnostic_name &
      is.na(diagnostics$state),
    ,
    drop = FALSE
  ]
  if (nrow(row) == 0L) {
    return(NULL)
  }
  row[1L, , drop = FALSE]
}

diagnostic_subtitle <- function(row) {
  if (is.null(row) || nrow(row) == 0L) {
    return("")
  }
  paste0(
    "p = ",
    format(round(row$p_value, 3), nsmall = 3),
    "; ",
    row$interpretation_label,
    "; ",
    row$sharpness_label
  )
}

empty_panel <- function(title, subtitle = "") {
  ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_blank() +
    ggplot2::theme_void() +
    ggplot2::labs(title = title, subtitle = subtitle)
}

plot_msd_dashboard <- function(x) {
  rows <- list()
  for (method in names(x$simulations)) {
    sim <- x$simulations[[method]]
    obs <- msd_curve(sim$observed_track)
    rows[[length(rows) + 1L]] <- data.frame(method = method, source = "observed", lag = obs$lag, msd = obs$msd)
    for (i in seq_along(sim$simulated_tracks)) {
      curve <- msd_curve(sim$simulated_tracks[[i]], max_lag = max(obs$lag))
      rows[[length(rows) + 1L]] <- data.frame(method = method, source = "simulated", sim_id = i, lag = curve$lag, msd = curve$msd)
    }
  }
  df <- do.call(rbind, rows)
  ggplot2::ggplot() +
    ggplot2::geom_line(
      data = df[df$source == "simulated", ],
      ggplot2::aes(x = lag, y = msd, group = interaction(method, sim_id)),
      color = "#4C78A8",
      alpha = 0.2
    ) +
    ggplot2::geom_line(
      data = df[df$source == "observed", ],
      ggplot2::aes(x = lag, y = msd),
      color = "#D55E00",
      linewidth = 0.9
    ) +
    ggplot2::facet_wrap(~method) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Lag", y = "Mean squared displacement")
}

plot_state_occupancy_dashboard <- function(x) {
  df <- x$diagnostics[x$diagnostics$diagnostic == "state_occupancy" & !is.na(x$diagnostics$state), , drop = FALSE]
  ggplot2::ggplot(df, ggplot2::aes(x = factor(state), y = observed_statistic, fill = method)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::geom_point(ggplot2::aes(y = simulated_median), position = ggplot2::position_dodge(width = 0.9), color = "black") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "State", y = "Occupancy discrepancy", fill = "Method")
}

plot_residence_time_dashboard <- function(x) {
  df <- x$diagnostics[x$diagnostics$diagnostic == "state_residence_time", , drop = FALSE]
  ggplot2::ggplot(df, ggplot2::aes(x = factor(state), y = simulated_median, fill = interpretation_label)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~method) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "State", y = "Residence-time discrepancy", fill = "Interpretation")
}

plot_method_comparison <- function(x) {
  df <- x$diagnostics
  ggplot2::ggplot(df, ggplot2::aes(x = method, y = diagnostic, fill = interpretation_label)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Simulation method", y = "Diagnostic", fill = "Interpretation")
}
