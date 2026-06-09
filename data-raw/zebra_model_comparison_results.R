# Generate embedded zebra model-comparison results for the pkgdown vignette.
#
# This script is intentionally kept outside the built package. It requires the
# local hmmSSF repository and can be slow because it fits and simulates several
# zebra models.

library(hmmSSF)
library(raster)
library(survival)

script_file <- tryCatch(
  normalizePath(sys.frame(1)$ofile),
  error = function(error) NA_character_
)
if (is.na(script_file) || !nzchar(script_file)) {
  script_file <- normalizePath(
    file.path(getwd(), "data-raw", "zebra_model_comparison_results.R"),
    mustWork = FALSE
  )
}

package_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
root_dir <- normalizePath(
  Sys.getenv("VALIDHMMSSF_ROOT", file.path(package_dir, "..")),
  mustWork = FALSE
)
hmmssf_repo <- normalizePath(
  Sys.getenv("HMMSSF_REPO", file.path(root_dir, "hmmSSF")),
  mustWork = FALSE
)
out_dir <- Sys.getenv(
  "ZEBRA_MODEL_COMPARISON_OUT_DIR",
  file.path(root_dir, "zebra_validation", "model_comparison")
)
extdata_dir <- file.path(package_dir, "inst", "extdata")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(extdata_dir, showWarnings = FALSE, recursive = TRUE)

devtools::load_all(package_dir, quiet = TRUE)

source(file.path(hmmssf_repo, "inst/functions/cov_df.R"))
source(file.path(hmmssf_repo, "inst/functions/simHMMSSF.R"))

n_validate <- as.integer(Sys.getenv("ZEBRA_MODEL_COMPARISON_N_VALIDATE", "250"))
n_sims <- as.integer(Sys.getenv("ZEBRA_MODEL_COMPARISON_N_SIMS", "99"))
n_zeros <- as.integer(Sys.getenv("ZEBRA_MODEL_COMPARISON_N_ZEROS", "500"))
rmax <- as.numeric(Sys.getenv("ZEBRA_MODEL_COMPARISON_RMAX", "5"))
optim_maxit <- as.integer(Sys.getenv("ZEBRA_MODEL_COMPARISON_MAXIT", "1500"))
sim_max_attempts <- as.integer(Sys.getenv("ZEBRA_MODEL_COMPARISON_SIM_ATTEMPTS", "50"))
sim_cores <- as.integer(Sys.getenv(
  "ZEBRA_MODEL_COMPARISON_CORES",
  as.character(min(4L, max(1L, parallel::detectCores() - 1L)))
))
reuse_fits <- tolower(Sys.getenv("ZEBRA_MODEL_COMPARISON_REUSE_FITS", "true")) %in%
  c("true", "1", "yes")
seed <- as.integer(Sys.getenv("ZEBRA_MODEL_COMPARISON_SEED", "250"))

set.seed(seed)

model_data_file <- Sys.getenv(
  "ZEBRA_MODEL_DATA_FILE",
  file.path(root_dir, "zebra_validation", "zebra_model_data.rds")
)
if (!file.exists(model_data_file)) {
  stop("Missing model data file: ", model_data_file, call. = FALSE)
}

model_data <- readRDS(model_data_file)
habitat <- raster::raster(file.path(hmmssf_repo, "inst/zebra/vegetation2.grd"))
observed_full <- subset(model_data, obs == 1)
complete_observed_index <- which(is.finite(observed_full$x) & is.finite(observed_full$y))
n_validate <- min(n_validate, length(complete_observed_index))
validation_index <- complete_observed_index[seq_len(n_validate)]
observed_validate <- observed_full[validation_index, c("x", "y", "time", "tod")]
observed_track <- data.frame(x = observed_validate$x, y = observed_validate$y)

ssf_formula <- ~ step + log(step) + cos(angle) + veg
tpm_fixed <- ~ 1
tpm_time <- ~ cos(2 * pi * tod / 24) + sin(2 * pi * tod / 24)

initial_par_list <- readRDS(file.path(hmmssf_repo, "inst/zebra/initial_par.RData"))
initial_par_raw <- initial_par_list[[3]]

make_three_state_betas <- function(two_state_betas) {
  beta_jitter <- seq(-0.05, 0.05, length.out = nrow(two_state_betas))
  cbind(
    two_state_betas,
    rowMeans(two_state_betas) + beta_jitter
  )
}

make_tpm_start <- function(n_states, tpm_formula, observed_data, raw_alphas) {
  tpm_mm <- stats::model.matrix(tpm_formula, observed_data)
  n_tpm_cov <- ncol(tpm_mm)
  n_cols <- n_states * (n_states - 1L)
  out <- matrix(0, nrow = n_tpm_cov, ncol = n_cols)

  if (n_states == 2L && n_tpm_cov == nrow(raw_alphas)) {
    out[, ] <- raw_alphas
  } else if (n_states == 2L && n_tpm_cov == 1L) {
    out[1L, ] <- raw_alphas[1L, ]
  } else if (n_states == 3L && n_tpm_cov == 1L) {
    out[1L, 1L] <- raw_alphas[1L, 1L]
    out[1L, 3L] <- raw_alphas[1L, 2L]
  }

  out
}

fit_one_state_ssf <- function(data, ssf_formula) {
  formula <- stats::as.formula(
    paste0("obs ~ ", deparse(ssf_formula[[2L]]), " + strata(stratum)")
  )
  fit <- survival::clogit(formula, data = data)
  betas <- matrix(stats::coef(fit), ncol = 1)
  rownames(betas) <- names(stats::coef(fit))
  colnames(betas) <- "S1"
  list(
    par = list(
      ssf = betas,
      tpm = matrix(numeric(0), nrow = 1, ncol = 0)
    ),
    fit = fit,
    args = list(
      ssf_formula = ssf_formula,
      tpm_formula = ~ 1,
      data = data,
      n_states = 1L
    )
  )
}

fit_hmmssf_model <- function(n_states, tpm_formula, ssf_par0, tpm_par0, data) {
  hmmSSF::hmmSSF(
    ssf_formula = ssf_formula,
    tpm_formula = tpm_formula,
    data = data,
    ssf_par0 = ssf_par0,
    tpm_par0 = tpm_par0,
    n_states = n_states,
    optim_opts = list(trace = 0, maxit = optim_maxit)
  )
}

fit_hmmssf_multistart <- function(model_name, n_states, tpm_formula, initial_pars, data) {
  message("Fitting ", model_name, " with ", length(initial_pars), " starts")

  fit_attempts <- lapply(seq_along(initial_pars), function(start_id) {
    message("  start ", start_id, "/", length(initial_pars))
    fit <- try(
      fit_hmmssf_model(
        n_states = n_states,
        tpm_formula = tpm_formula,
        ssf_par0 = initial_pars[[start_id]]$betas,
        tpm_par0 = initial_pars[[start_id]]$alphas,
        data = data
      ),
      silent = TRUE
    )

    if (inherits(fit, "try-error")) {
      return(list(
        start = start_id,
        fit = NULL,
        convergence = NA_integer_,
        objective = Inf,
        error = conditionMessage(attr(fit, "condition"))
      ))
    }

    list(
      start = start_id,
      fit = fit,
      convergence = fit$fit$convergence,
      objective = fit$fit$value,
      error = NA_character_
    )
  })

  start_results <- do.call(
    rbind,
    lapply(fit_attempts, function(attempt) {
      data.frame(
        start = attempt$start,
        convergence = attempt$convergence,
        objective = attempt$objective,
        error = attempt$error,
        stringsAsFactors = FALSE
      )
    })
  )

  successful <- is.finite(start_results$objective)
  if (!any(successful)) {
    stop(
      "All starting values failed for ", model_name, ": ",
      paste(start_results$error, collapse = "; "),
      call. = FALSE
    )
  }

  best_row <- which.min(start_results$objective)
  best_fit <- fit_attempts[[best_row]]$fit
  best_fit$generative_repair <- list(
    fit_strategy = "multi_start_lowest_objective",
    n_starts = length(initial_pars),
    selected_start = start_results$start[best_row],
    start_results = start_results
  )

  message(
    "Selected start ", start_results$start[best_row],
    " for ", model_name,
    " with objective ", signif(start_results$objective[best_row], 8)
  )

  best_fit
}

model_specs <- list(
  M1_one_state = list(
    label = "M1 one-state SSF",
    n_states = 1L,
    tpm_formula = ~ 1,
    fit_type = "one_state_ssf"
  ),
  M2_two_state_fixed = list(
    label = "M2 two-state fixed TPM",
    n_states = 2L,
    tpm_formula = tpm_fixed,
    fit_type = "hmmssf",
    ssf_par0 = initial_par_raw$betas,
    tpm_par0 = make_tpm_start(2L, tpm_fixed, observed_validate, initial_par_raw$alphas)
  ),
  M3_two_state_time = list(
    label = "M3 two-state time TPM",
    n_states = 2L,
    tpm_formula = tpm_time,
    fit_type = "hmmssf",
    fit_strategy = "multi_start_lowest_objective",
    initial_pars = initial_par_list
  ),
  M4_three_state_fixed = list(
    label = "M4 three-state fixed TPM",
    n_states = 3L,
    tpm_formula = tpm_fixed,
    fit_type = "hmmssf",
    ssf_par0 = make_three_state_betas(initial_par_raw$betas),
    tpm_par0 = make_tpm_start(3L, tpm_fixed, observed_validate, initial_par_raw$alphas)
  )
)

fit_model <- function(model_name) {
  spec <- model_specs[[model_name]]
  fit_file <- file.path(out_dir, paste0(model_name, "_fit.rds"))
  if (reuse_fits && file.exists(fit_file)) {
    cached_fit <- readRDS(fit_file)
    if (!identical(spec$fit_strategy, "multi_start_lowest_objective") ||
        identical(
          cached_fit$generative_repair$fit_strategy,
          "multi_start_lowest_objective"
        )) {
      return(cached_fit)
    }
    message(
      "Ignoring cached ", model_name,
      " fit because it was not created with the article-aligned multi-start strategy."
    )
  }
  message("Fitting ", model_name)
  fit <- if (identical(spec$fit_type, "one_state_ssf")) {
    fit_one_state_ssf(model_data, ssf_formula)
  } else if (identical(spec$fit_strategy, "multi_start_lowest_objective")) {
    fit_hmmssf_multistart(
      model_name = model_name,
      n_states = spec$n_states,
      tpm_formula = spec$tpm_formula,
      initial_pars = spec$initial_pars,
      data = model_data
    )
  } else {
    fit_hmmssf_model(
      n_states = spec$n_states,
      tpm_formula = spec$tpm_formula,
      ssf_par0 = spec$ssf_par0,
      tpm_par0 = spec$tpm_par0,
      data = model_data
    )
  }
  saveRDS(fit, fit_file)
  fit
}

ssf_cov <- list(veg = habitat)
levels_cov <- list(veg = levels(model_data$veg))
y1 <- as.numeric(observed_validate[1L, c("x", "y")])

simulate_one_markov <- function(model_name, fit, sim_id) {
  spec <- model_specs[[model_name]]
  fit_par <- list(betas = fit$par$ssf, alphas = fit$par$tpm)
  last_error <- NULL
  for (attempt in seq_len(sim_max_attempts)) {
    sim <- try(
      simHMMSSF(
        ssf_formula = ssf_formula,
        tpm_formula = spec$tpm_formula,
        ssf_cov = ssf_cov,
        tpm_cov = observed_validate,
        par = fit_par,
        n_states = spec$n_states,
        n_obs = n_validate,
        n_zeros = n_zeros,
        y1 = y1,
        rmax = rmax,
        levels_cov = levels_cov,
        print = FALSE
      ),
      silent = TRUE
    )
    if (inherits(sim, "try-error")) {
      last_error <- conditionMessage(attr(sim, "condition"))
      next
    }
    if (is.data.frame(sim) &&
        nrow(sim) == n_validate &&
        all(is.finite(sim$x)) &&
        all(is.finite(sim$y))) {
      sim$x[1L] <- observed_track$x[1L]
      sim$y[1L] <- observed_track$y[1L]
      return(sim)
    }
  }
  stop(
    "Could not simulate ", model_name, " trajectory ", sim_id,
    " after ", sim_max_attempts, " attempts. Last error: ", last_error,
    call. = FALSE
  )
}

simulate_model <- function(model_name, fit) {
  sim_file <- file.path(out_dir, paste0(model_name, "_markov_simulations.rds"))
  fit_file <- file.path(out_dir, paste0(model_name, "_fit.rds"))
  fit_mtime <- if (file.exists(fit_file)) {
    file.info(fit_file)$mtime
  } else {
    as.POSIXct(0, origin = "1970-01-01")
  }

  if (reuse_fits &&
      file.exists(sim_file) &&
      file.info(sim_file)$mtime >= fit_mtime) {
    return(readRDS(sim_file))
  }

  message("Simulating ", model_name)
  sim_part_dir <- file.path(out_dir, paste0(model_name, "_markov_simulation_parts"))
  dir.create(sim_part_dir, showWarnings = FALSE, recursive = TRUE)

  simulate_and_save <- function(i) {
    part_file <- file.path(sim_part_dir, sprintf("sim_%03d.rds", i))
    if (reuse_fits &&
        file.exists(part_file) &&
        file.info(part_file)$mtime >= fit_mtime) {
      return(readRDS(part_file))
    }
    set.seed(seed + 10000L * which(names(model_specs) == model_name) + i)
    sim <- simulate_one_markov(model_name, fit, i)
    saveRDS(sim, part_file)
    sim
  }

  sim_ids <- seq_len(n_sims)
  sims <- vector("list", length(sim_ids))
  batch_size <- max(sim_cores, 1L)
  batches <- split(sim_ids, ceiling(seq_along(sim_ids) / batch_size))
  for (batch_index in seq_along(batches)) {
    batch <- batches[[batch_index]]
    message(
      "  batch ", batch_index, "/", length(batches),
      " for ", model_name,
      " (simulations ", min(batch), "-", max(batch), ")"
    )
    batch_sims <- if (sim_cores > 1L) {
      parallel::mclapply(
        batch,
        simulate_and_save,
        mc.cores = min(sim_cores, length(batch))
      )
    } else {
      lapply(batch, simulate_and_save)
    }
    sims[batch] <- batch_sims
  }
  simulated_tracks <- lapply(sims, function(x) data.frame(x = x$x, y = x$y))
  names(simulated_tracks) <- paste0("sim_", seq_len(n_sims))
  simulated_states <- do.call(rbind, lapply(sims, function(x) as.integer(x$state)))
  out <- list(simulated_tracks = simulated_tracks, simulated_states = simulated_states)
  saveRDS(out, sim_file)
  out
}

get_observed_states <- function(model_name, fit) {
  spec <- model_specs[[model_name]]
  if (spec$n_states == 1L) {
    return(rep(1L, nrow(observed_track)))
  }
  hmmSSF::viterbi_decoding(fit)[validation_index]
}

fit_summary_row <- function(model_name, fit) {
  spec <- model_specs[[model_name]]
  fit_metadata <- fit$generative_repair
  convergence <- NA_integer_
  objective <- NA_real_
  fit_strategy <- "single_start"
  n_starts <- 1L
  selected_start <- NA_integer_

  if (identical(spec$fit_type, "one_state_ssf")) {
    convergence <- NA_integer_
    objective <- -as.numeric(stats::logLik(fit$fit))
    fit_strategy <- "conditional_logistic"
  } else {
    convergence <- fit$fit$convergence
    objective <- fit$fit$value
    if (is.list(fit_metadata)) {
      fit_strategy <- fit_metadata$fit_strategy
      n_starts <- fit_metadata$n_starts
      selected_start <- fit_metadata$selected_start
    }
  }

  data.frame(
    model = model_name,
    label = spec$label,
    n_states = spec$n_states,
    tpm = deparse(spec$tpm_formula),
    generation = "markov",
    n_validate = n_validate,
    n_sims = n_sims,
    n_zeros = n_zeros,
    rmax = rmax,
    convergence = convergence,
    objective = objective,
    fit_strategy = fit_strategy,
    n_starts = n_starts,
    selected_start = selected_start,
    stringsAsFactors = FALSE
  )
}

fits <- lapply(names(model_specs), fit_model)
names(fits) <- names(model_specs)
simulations <- Map(simulate_model, names(fits), fits)

diagnoses <- lapply(names(fits), function(model_name) {
  sim <- simulations[[model_name]]
  diagnose_hmmssf_simulations(
    observed_track = observed_track,
    simulated_tracks = sim$simulated_tracks,
    observed_states = get_observed_states(model_name, fits[[model_name]]),
    simulated_states = sim$simulated_states,
    label = model_name,
    diagnostics = c(
      "ud_wasserstein",
      "msd",
      "sinuosity",
      "state_occupancy",
      "state_residence_time",
      "switching_rate",
      "transition_counts",
      "state_conditioned_step_length",
      "state_conditioned_turning_angle",
      "state_conditioned_msd"
    ),
    metadata = list(
      simulation_label = model_specs[[model_name]]$label,
      simulation_description = "Markov-generated latent states and SSF-generated spatial steps.",
      generation = "markov",
      n_states = model_specs[[model_name]]$n_states,
      tpm_formula = deparse(model_specs[[model_name]]$tpm_formula)
    )
  )
})
names(diagnoses) <- names(fits)

diagnostics <- do.call(
  rbind,
  lapply(names(diagnoses), function(model_name) {
    out <- diagnoses[[model_name]]$diagnostics
    out$model <- model_name
    out$model_label <- model_specs[[model_name]]$label
    out
  })
)
rownames(diagnostics) <- NULL

model_summary <- do.call(
  rbind,
  Map(fit_summary_row, names(fits), fits)
)

global_diagnostics <- diagnostics[
  is.na(diagnostics$state) &
    diagnostics$diagnostic %in% c("ud_wasserstein", "msd", "sinuosity"),
  ,
  drop = FALSE
]

utils::write.csv(
  diagnostics,
  file.path(extdata_dir, "zebra_model_comparison_diagnostics.csv"),
  row.names = FALSE
)
utils::write.csv(
  model_summary,
  file.path(extdata_dir, "zebra_model_comparison_model_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  global_diagnostics,
  file.path(extdata_dir, "zebra_model_comparison_global_diagnostics.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    model_summary = model_summary,
    diagnostics = diagnostics,
    diagnoses = diagnoses
  ),
  file.path(out_dir, "zebra_model_comparison_results.rds")
)

print(model_summary)
print(global_diagnostics[, c("model", "diagnostic", "p_value", "observed_statistic", "simulated_median")])
