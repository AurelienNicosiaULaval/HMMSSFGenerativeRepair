# HMMSSFGenerativeRepair

Minimal R package for generative validation of hidden Markov
step-selection functions (HMM-SSF). The package simulates trajectories from
fitted HMM-SSF-like models, compares observed and simulated movement patterns,
and reports Monte Carlo diagnostic p-values.

The package is designed as a private GitHub research package. It keeps the
simulation strategy explicit because different strategies answer different
diagnostic questions.

## Installation

From a local checkout:

```r
devtools::install("HMMSSF_GenerativeRepair")
```

After creating a private GitHub repository under the expected name:

```r
remotes::install_github("AurelienNicosiaULaval/HMMSSFGenerativeRepair")
```

If the private repository uses a different name, replace the owner and
repository in the command above.

## Quick start

```r
library(HMMSSFGenerativeRepair)

fit <- example_hmmssf_list(n = 60, seed = 123)
observed_track <- fit$observed_track

diagnosis <- diagnose_hmmssf(
  fit = fit,
  observed_track = observed_track,
  n_sims = 99,
  methods = c("markov", "viterbi", "viterbi_tube", "posterior"),
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
  parameter_uncertainty = FALSE,
  epsilon = 10,
  seed = 123
)

summary(diagnosis)
plot(diagnosis, method = "markov")
```

The main outputs are:

- `observed_statistic`: diagnostic statistic for the observed trajectory.
- `simulated_median`: median of the simulated diagnostic statistics.
- `simulated_q025` and `simulated_q975`: simulation envelope.
- `p_value`: Monte Carlo p-value.
- `interpretation_label`: compact diagnostic label.

## Four simulation strategies

`diagnose_hmmssf()` can use one or several simulation strategies:

```r
methods = c("markov", "viterbi", "viterbi_tube", "posterior")
```

The strategies differ only in how the latent state path is generated. Once a
state path has been chosen, the package simulates spatial steps using the
state-specific movement and selection kernel.

| Method | Latent state path | Main diagnostic question |
| --- | --- | --- |
| `markov` | Simulated from the fitted transition model | Does the full HMM-SSF generative model reproduce the trajectory? |
| `viterbi` | Fixed to the Viterbi path | If the decoded states are fixed, do the state-specific SSF kernels reproduce movement? |
| `viterbi_tube` | Sampled among paths close to Viterbi | Are conclusions sensitive to near-optimal latent paths? |
| `posterior` | Sampled from posterior state information | Are conclusions sensitive to latent-state uncertainty? |

For a model with one state, all four latent-state strategies collapse to the
same single-state path. The relevant validation then becomes the original
trajectory-level generative validation: observed trajectory versus simulated
trajectories, mainly through UD, MSD, sinuosity, and other global movement
diagnostics.

## pkgdown site

The package includes a minimal pkgdown configuration. To build the site locally:

```r
pkgdown::build_site("HMMSSF_GenerativeRepair")
```

The main article is:

```r
vignette("hmmssf-generative-validation")
```

The zebra model-comparison article is:

```r
vignette("zebra-model-comparison")
```

## Private GitHub setup

Example SSH workflow:

```sh
cd HMMSSF_GenerativeRepair
git init
git add .
git commit -m "Initial HMM-SSF generative validation package"
git remote add origin git@github.com:AurelienNicosiaULaval/HMMSSFGenerativeRepair.git
git push -u origin main
```

Create the GitHub repository as private before pushing.
