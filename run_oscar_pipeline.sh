#!/bin/bash
# ==============================================================================
# NEC model pipeline on OSCAR
#
#   01a  analysis/01_oscar_simulate_uncertainty_sets.R    ─┐
#   01b  analysis/01_oscar_simulate_uncertainty_bounds.R   │ (independent)
#   02   analysis/02_learn_GAM_emulator.R      after 01a succeeds
#   03   analysis/03_calculate_tx_differences.R after 02 succeeds
#
# Usage (run on a login node from anywhere; quote set names containing '&'):
#   bash run_oscar_pipeline.sh COMPETE
#   bash run_oscar_pipeline.sh 'NETTER&RADIANT'
#
# The parameter set selects Data/input_params_<SET>.csv; all outputs go to
# Data/simulation_outputs/<SET>/ (overwritten on each run).
# To run only some steps, set STEPS, e.g. STEPS="02 03" bash run_oscar_pipeline.sh COMPETE
# ==============================================================================

set -euo pipefail

# ==============================================================================
# Settings
# ==============================================================================

PACKAGE_DIR="${PACKAGE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"   # this checkout
R_LIBRARY="${R_LIBRARY:-${HOME}/R/x86_64-pc-linux-gnu-library/4.6}"
R_MODULE="r/4.6.0"
PARTITION="batch"
QOS="priority+"
STEPS="${STEPS:-01a 01b 02 03}"

# Resources per step. Simulations use ~3.6 GB per core at N = 100,000.
SETS_CPUS=16;   SETS_MEM="96G";   SETS_TIME="06:00:00"
BOUNDS_CPUS=16; BOUNDS_MEM="96G"; BOUNDS_TIME="04:00:00"
GAM_CPUS=8;     GAM_MEM="32G";    GAM_TIME="02:00:00"
DIFF_CPUS=4;    DIFF_MEM="32G";   DIFF_TIME="01:00:00"

# ==============================================================================
# Parameter set and paths
# ==============================================================================

if [[ $# -ne 1 ]]; then
  echo "Usage: bash run_oscar_pipeline.sh <PARAM_SET>   (e.g. COMPETE or 'NETTER&RADIANT')" >&2
  exit 1
fi

PARAM_SET="$1"
# Safe version of the set name for SLURM job and log names
SET_TAG="${PARAM_SET//[^A-Za-z0-9]/_}"

export NEC_PACKAGE_DIR="${PACKAGE_DIR}"
export NEC_INPUT_FILE="${PACKAGE_DIR}/Data/input_params_${PARAM_SET}.csv"
export NEC_DESIGN_FILE="${PACKAGE_DIR}/design/maxPro_design_3000.rds"
export NEC_OUTPUT_DIR="${PACKAGE_DIR}/Data/simulation_outputs/${PARAM_SET}"
export R_LIBS_USER="${R_LIBRARY}"
LOG_DIR="${PACKAGE_DIR}/logs"

for f in "${NEC_INPUT_FILE}" "${NEC_DESIGN_FILE}"; do
  if [[ ! -f "${f}" ]]; then
    echo "File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "${R_LIBRARY}" "${LOG_DIR}" "${NEC_OUTPUT_DIR}"

# ==============================================================================
# Load R and install packages
# ==============================================================================

module load "${R_MODULE}"
echo "Using R library: ${R_LIBS_USER}"

Rscript - <<'RSCRIPT'
lib <- Sys.getenv("R_LIBS_USER")
.libPaths(c(lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

install_if_missing <- function(pkg, repo = NULL) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    return(invisible())
  }
  message("Installing ", pkg, "...")
  if (is.null(repo)) {
    install.packages(pkg, lib = lib)
  } else {
    github_pat <- Sys.getenv("GITHUB_PAT")
    if (!nzchar(github_pat)) {
      stop("GITHUB_PAT is not available. Add it to ~/.Renviron before running this script.")
    }
    remotes::install_github(repo, auth_token = github_pat, lib = lib,
                            upgrade = "never", dependencies = TRUE)
  }
}

# CRAN packages
for (pkg in c("remotes", "cmprsk", "data.table", "mgcv", "R6", "survival", "truncnorm")) {
  install_if_missing(pkg)
}

# GitHub packages
install_if_missing("nhppp", "bladder-ca/nhppp")
install_if_missing("mortality", "bladder-ca/mortality")
RSCRIPT

# Always reinstall NECModel from this checkout, so the jobs run the code in
# PACKAGE_DIR rather than an older installed copy.
echo "Installing NECModel from ${PACKAGE_DIR}"
R CMD INSTALL --no-docs --no-multiarch -l "${R_LIBRARY}" "${PACKAGE_DIR}" > "${LOG_DIR}/${SET_TAG}_install.log" 2>&1 || {
  echo "NECModel installation failed; see ${LOG_DIR}/${SET_TAG}_install.log" >&2
  exit 1
}

# Record what is being run
GIT_COMMIT="$(git -C "${PACKAGE_DIR}" rev-parse --short HEAD 2>/dev/null || echo unknown)"
GIT_DIRTY="$( (git -C "${PACKAGE_DIR}" status --porcelain -- R DESCRIPTION NAMESPACE analysis 2>/dev/null || true) | head -1)"
if [[ -n "${GIT_DIRTY}" ]]; then
  echo "Warning: uncommitted changes in R/, analysis/ or package files; outputs will not match commit ${GIT_COMMIT}." >&2
fi

# ==============================================================================
# Submit jobs
# ==============================================================================

# Submit one step as a SLURM job and print its job ID.
# Arguments: step name, R script, CPUs, memory, time, dependency job IDs (optional)
submit() {
  local name="$1" script="$2" cpus="$3" mem="$4" time="$5" deps="${6:-}"
  local dependency=()
  if [[ -n "${deps}" ]]; then
    dependency=(--dependency="afterok:${deps}" --kill-on-invalid-dep=yes)
  fi
  sbatch --parsable \
    --job-name="nec_${SET_TAG}_${name}" \
    --partition="${PARTITION}" \
    --qos="${QOS}" \
    --nodes=1 --ntasks=1 \
    --cpus-per-task="${cpus}" \
    --mem="${mem}" \
    --time="${time}" \
    --output="${LOG_DIR}/${SET_TAG}_${name}_%j.out" \
    --error="${LOG_DIR}/${SET_TAG}_${name}_%j.err" \
    --export=ALL \
    ${dependency[@]+"${dependency[@]}"} \
    --wrap="module load ${R_MODULE}; Rscript --vanilla ${PACKAGE_DIR}/analysis/${script}"
}

has_step() { [[ " ${STEPS} " == *" $1 "* ]]; }

JOB_SETS=""; JOB_BOUNDS=""; JOB_GAM=""; JOB_DIFF=""

if has_step 01a; then
  JOB_SETS=$(submit sets 01_oscar_simulate_uncertainty_sets.R "${SETS_CPUS}" "${SETS_MEM}" "${SETS_TIME}")
fi
if has_step 01b; then
  JOB_BOUNDS=$(submit bounds 01_oscar_simulate_uncertainty_bounds.R "${BOUNDS_CPUS}" "${BOUNDS_MEM}" "${BOUNDS_TIME}")
fi
if has_step 02; then
  JOB_GAM=$(submit gam 02_learn_GAM_emulator.R "${GAM_CPUS}" "${GAM_MEM}" "${GAM_TIME}" "${JOB_SETS}")
fi
if has_step 03; then
  JOB_DIFF=$(submit diff 03_calculate_tx_differences.R "${DIFF_CPUS}" "${DIFF_MEM}" "${DIFF_TIME}" "${JOB_GAM}")
fi

cat > "${NEC_OUTPUT_DIR}/run_info.txt" <<EOF
parameter_set: ${PARAM_SET}
submitted: $(date '+%Y-%m-%d %H:%M:%S')
git_commit: ${GIT_COMMIT}${GIT_DIRTY:+ (with uncommitted changes)}
input_file: ${NEC_INPUT_FILE}
steps: ${STEPS}
jobs: sets=${JOB_SETS:-none} bounds=${JOB_BOUNDS:-none} gam=${JOB_GAM:-none} diff=${JOB_DIFF:-none}
EOF

echo "Parameter set: ${PARAM_SET} (commit ${GIT_COMMIT})"
echo "  01a sets:        ${JOB_SETS:-not submitted}"
echo "  01b bounds:      ${JOB_BOUNDS:-not submitted}"
echo "  02  GAM:         ${JOB_GAM:-not submitted}${JOB_SETS:+ (after ${JOB_SETS})}"
echo "  03  differences: ${JOB_DIFF:-not submitted}${JOB_GAM:+ (after ${JOB_GAM})}"
echo "Outputs: ${NEC_OUTPUT_DIR}"
echo "Monitor with: squeue -u \$USER"
