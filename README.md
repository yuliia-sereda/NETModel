# NETModel: treatment sequencing in small bowel neuroendocrine tumours

A continuous-time microsimulation model comparing sequences of peptide receptor
radionuclide therapy (PRRT; **P**), everolimus (**E**), and cabozantinib (**C**)
for patients with small bowel neuroendocrine tumours (SB-NETs) who progress on
somatostatin analogues, and a decision analysis under deep uncertainty built on it.

Delaying PRRT may reduce the risk of post-PRRT myelodysplastic syndrome and acute
myeloid leukaemia (MDS/AML) but lessen its effectiveness; the effectiveness of
everolimus and cabozantinib may also decline with later use. The model propagates
these trade-offs, which are deeply uncertain, to mean overall survival (OS),
progression-free survival (PFS), and MDS/AML incidence.

## Model

### Decision problem

- **Population:** a hypothetical population of patients with SB-NETs who progressed
  on somatostatin analogues and are eligible for everolimus, PRRT, or cabozantinib;
  50% female, age at entry normally distributed (mean 63, SD 9 years, truncated to
  18–84 years; based on NETTER-1), birth cohort 1984 (`make_population()`).
- **Strategies:** the 6 orderings of all three drugs (PEC, PCE, EPC, ECP, CPE, CEP)
  and the 2 orderings without cabozantinib (PE, EP), which may be unavailable or not
  standard of care in some settings (`treatment_sequences()`). Patients start on the
  first treatment and move to the next upon disease progression (RECIST 1.1).
- **Treatment:** everolimus and cabozantinib are given continuously until progression
  or death; PRRT is given as 4 doses, one every 2 months. Each PRRT dose permanently
  increases the MDS/AML hazard rate by a constant amount; doses stop accruing when
  PRRT stops (e.g. at progression).
- **Time horizon:** lifetime, with a maximum age of 100 years.
- **Outcomes:** mean OS (time to death), mean PFS (time to the earliest of salvage
  treatment initiation, MDS/AML, and death), and the 5-year cumulative incidence of
  MDS/AML (death as a competing event).

### States and transitions

All patients start with stable disease on the first treatment of a sequence
(state 1). States are mutually exclusive and exhaustive:

| State | Description |
|---|---|
| 1, 2, 3 | Stable disease on the 1st, 2nd, 3rd treatment |
| 4 | Salvage treatment (after failure of the last treatment) |
| 5 | MDS/AML |
| 6 | Death from disease progression (only after starting salvage treatment) |
| 7 | Death from MDS/AML (only after developing MDS/AML) |
| 8 | Death from other causes (possible at any time) |

Patients who develop MDS/AML are treated for the complication: all SB-NET
treatments stop, and no further treatments in the sequence are started. Two-drug strategies move from state 2 directly to
salvage.

The model simulates potential times for each transition and the earliest happens
first. Transition times are drawn from separate non-homogeneous Poisson point
processes (NHPPPs), generated from intensity (hazard rate) functions with the
[`nhppp`](https://github.com/bladder-ca/nhppp) package. Ages of death from other
causes are drawn from sex-, age-, and birth-cohort-specific Human Mortality Database
rates with the [`mortality`](https://github.com/bladder-ca/mortality) package.
Each person is simulated under every strategy, and simulations use common random
numbers to reduce simulation error when comparing strategies and design points.

### Input parameters

Parameters are read from `Data/input_params_<set>.csv`. Two parameter sets differ
in the sources for first-line PFS:

| Parameter | Transition | `NETTER&RADIANT` | `COMPETE` |
|---|---|---|---|
| `lambda_everolimus_12` | Progression on everolimus | RADIANT-4: median PFS 11.0 (9.2–13.3) months | COMPETE: 12.0 (8.8–23.7) months |
| `lambda_prrt_12` | Progression on PRRT | NETTER-1: median PFS 28.4 (28.4–NE) months | COMPETE: 23.9 (14.0–NE) months |
| `lambda_cabozantinib_23` | Progression on cabozantinib | CABINET: median PFS 8.4 (7.6–12.7) months | same |
| `lambda_57` | Death from MDS/AML | Median OS for very-high IPSS-M: 12 (6.0–21.6) months | same |
| `lambda_46` | Death from disease progression | Calibrated (see below) | Calibrated |
| `lambda_dot5` | MDS/AML, per PRRT dose | Calibrated (see below) | Calibrated |

Hazard rates are constant, except the MDS/AML hazard (stepwise, one step per PRRT
dose) and other-cause mortality. Rates fitted to median survival times assume the
median is lognormally distributed with matched mean and variance, so the rate
λ = log(2) / median is also lognormal; rates given directly as incidence rates are
lognormal on the log-rate scale (`sample_lambda()`, Supplement A.2.2). Rates are
sampled once per person and reused across strategies and design points.

`lambda_46` and `lambda_dot5` are calibrated jointly to the median OS in NETTER-1
(48.0, 95% CI 37.4–55.2 months) and the 5-year MDS/AML cumulative incidence in the
Washington University in St. Louis data (5.1%, 95% CI 1.8–8.4%), under the PE
sequence, using a Latin hypercube design and GAM emulators solved for the targets.

### Uncertainty

- **Probabilistic uncertainty:** the lognormal distributions above
  (`analysis/calculate_lambda_dist_moments.R` gives their parameters for Table 1).
- **Deep uncertainty (uncertainty sets):** treatment-effect attenuation and MDS/AML
  incidence, which the trials do not identify. A treatment can lose up to 50% of its
  effectiveness with each subsequent position, i.e. its progression rate is up to
  1.5 times that in the previous position; the same attenuation applies to all three
  drugs. Progression rates in the 2nd and 3rd positions are the 1st-position rates
  multiplied by `discount1` ∈ [1, 1.5] and `discount2` ∈ [`discount1`,
  1.5 × `discount1`] (so at most 2.25); cabozantinib, observed in position 2, is
  adjusted relative to that position. The MDS/AML hazard is multiplied by
  `aml_mds_multiplier` (γ) ∈ [1, 3].
- **Propagation:** the model is run over a 3,000-point space-filling (maxPro)
  design of the per-position attenuation factors and γ, padded by 20% of each range
  beyond the set (`design/01_create_maxPro_design.R`). GAM emulators then predict
  each strategy's outcomes over a regular grid restricted to the uncertainty set.
- **Bounding scenarios:** favourable (no attenuation, γ = 1) and unfavourable
  (`discount1` = 1.5, `discount2` = 2.25, γ = 3), each simulated 1000 times.

### Choosing strategies

Under deep uncertainty each strategy's mean survival is known only up to an
interval. Strategies are compared with choice functions (`R/choice_functions.R`):
*not-uniformly-worse* strategies (maximality, interval dominance) remain after
excluding any strategy with shorter survival than a competitor in every state of
the world; the *best-worst-case* strategy (pessimism) has the longest survival in
its most adversarial state of the world. Optimism and Hurwicz criteria represent
other attitudes to ambiguity. Results are shown as Hasse diagrams, with and without
a 2-month clinical indifference margin.

## Repository structure

```
R/                   Package code
  Population.R         R6 class: simulates a population under each strategy
  intensity_functions.R  Hazard-rate (intensity) functions
  simulate_utils.R     Parameter sampling, sim_data(), make_population(), treatment_sequences()
  choice_functions.R   Bounds, choice functions, Hasse diagrams
analysis/            Analysis scripts (run from the repository root)
design/              maxPro design of the uncertainty set
Data/                Input parameters, calibration results, simulation outputs
Fig/                 Figures and tables
tests/               testthat tests
app/                 Shiny clinician app (see app/README.md)
run_oscar_pipeline.sh  Submits the simulation pipeline on OSCAR (SLURM)
```

## Workflow

| Step | Script | Where | Output |
|---|---|---|---|
| 0 | `analysis/00_calibrate_lambda46&dot5.R` | Local | Calibrated `lambda_46`, `lambda_dot5` per set (`Data/calibrate_lambda_46&dot5_*.csv`); enter them in `Data/input_params_<set>.csv` |
| 1a | `analysis/01_oscar_simulate_uncertainty_sets.R` | OSCAR | `NEC_data_uncertainty_sets.rds` |
| 1b | `analysis/01_oscar_simulate_uncertainty_bounds.R` | OSCAR | `NEC_data_uncertainty_bounds.rds` |
| 2 | `analysis/02_learn_GAM_emulator.R` (after 1a) | OSCAR | `NEC_GAM_predictions.rds` |
| 3 | `analysis/03_calculate_tx_differences.R` (after 2) | OSCAR | `NEC_GAM_differences_part_XX.rds` |
| 4 | `analysis/Pareto_front.R`, `bounds_figures.R`, `Hasse_graphs.R` | Local | `Fig/<set>/` |
| – | `analysis/calculate_lambda_dist_moments.R` | Local | `Fig/Table1_lambda_moments.csv` |

Simulation outputs are written to `Data/simulation_outputs/<set>/`.

### Running on OSCAR

On an OSCAR login node, from an up-to-date checkout (`git pull`):

```bash
bash run_oscar_pipeline.sh COMPETE
bash run_oscar_pipeline.sh 'NETTER&RADIANT'
```

The script installs `NECModel` from the checkout, submits steps 1a and 1b in
parallel, step 2 after 1a succeeds, and step 3 after 2 succeeds. Jobs use the
`priority+` QOS; resources are set at the top of the script. Outputs go to
`Data/simulation_outputs/<set>/` (with `run_info.txt`: commit, time, job IDs) and
logs to `logs/`. To re-run only some steps: `STEPS="02 03" bash run_oscar_pipeline.sh COMPETE`.
Installing the GitHub dependencies requires `GITHUB_PAT` in `~/.Renviron`.

Copy `Data/simulation_outputs/<set>/` to the same location locally to run step 4.

### Clinician-app simulations

A separate OSCAR pipeline simulates the outcomes shown in the clinician app, for
the COMPETE parameter set only, over a 180-point grid
(`design/02_create_app_design_grid.R`): attenuation of 0-50% at position 2 and,
on top of it, 0-50% at position 3 (in 10% steps), and MDS/AML multipliers
1-3 (in steps of 0.5). Parameter uncertainty is propagated with 1000 parameter
sets of a 100,000-person population each, the same sets at every design point.

```bash
bash run_oscar_app_pipeline.sh
```

The script submits a 25-task array (`analysis/app/01_oscar_app_simulate.R`),
then combines and summarises the results (`02_oscar_app_combine.R`) and
computes pairwise differences between strategies (`03_oscar_app_differences.R`).
Outputs go to `Data/app_outputs/COMPETE/`:

| File | Content |
|---|---|
| `app_outcomes_summary.rds/csv` | Mean and 95% interval across parameter sets of mean and median OS and PFS, OS and PFS at 1, 3 and 5 years, and MDS/AML cumulative incidence at 1, 3 and 5 years, for each strategy and design point |
| `app_differences_summary.rds/csv` | Differences between each pair of strategies (months or percentage points), with 95% intervals and the share of parameter sets in which the first strategy is better |
| `app_sets_part_XX.rds` | Set-level outcomes |
| `app_parameter_draws.rds` | Sampled parameter values of each set |


To update the clinician app (`app/`), copy `Data/app_outputs/COMPETE/` from
OSCAR and run, from the repository root:

```bash
Rscript analysis/app/04_prepare_app_data.R
```

This writes app-ready summaries to `app/data/`. The app only displays these
precomputed results and does not depend on NECModel, so only `app/` is deployed:

```r
rsconnect::deployApp("app")
```

## Installation

```r
# install.packages("remotes")
remotes::install_github("bladder-ca/nhppp")
remotes::install_github("bladder-ca/mortality")
devtools::install()   # or devtools::load_all(), from the repository root
```

Analysis scripts additionally use `mgcv`, `lhs`, `ggplot2`, `ggrepel`, `MASS`,
and `igraph`.

The clinician app is available at <https://ys-shiny.shinyapps.io/NEC-app/>.

## Citation

See `CITATION.cff`. Release notes are in `NEWS.md`.

## Authors

Thomas A. Trikalinos & Yuliia Sereda.
