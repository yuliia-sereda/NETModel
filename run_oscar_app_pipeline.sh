#!/bin/bash
# ==============================================================================
# NEC clinician-app pipeline on OSCAR (COMPETE parameter set only)
#
#   sim      analysis/app/01_oscar_app_simulate.R    SLURM array, one share of the
#                                                    parameter sets per task
#   combine  analysis/app/02_oscar_app_combine.R     after all sim tasks succeed
#   diff     analysis/app/03_oscar_app_differences.R after combine succeeds
#
# Simulates the 180-point app design grid (design/app_design_grid.rds) under
# parameter uncertainty: 1000 parameter sets, each with a 100,000-person
# population, at every design point. Outputs go to Data/app_outputs/COMPETE/
# (overwritten on each run), separate from the main pipeline's outputs.
#
# Usage (run on a login node from anywhere):
#   bash run_oscar_app_pipeline.sh
# To run only some steps, set STEPS, e.g. STEPS="combine diff" bash run_oscar_app_pipeline.sh
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
PARAM_SET="COMPETE"
STEPS="${STEPS:-sim combine diff}"

# Parameter uncertainty: number of parameter sets, split across array tasks
N_SETS=1000
N_TASKS=25
SEED=5678

# Resources. Each sim task runs N_SETS / N_TASKS sets on SIM_CPUS cores; each set
# runs all 180 design points (about 48-90 min per set at N = 100,000), and each
# core needs ~3.6 GB. 25 tasks x 24 cores = 600 cores, within the priority+
# limit of 624 CPUs and 4500 GB per user.
SIM_CPUS=24;    SIM_MEM="120G";    SIM_TIME="08:00:00"
COMBINE_CPUS=2; COMBINE_MEM="64G"; COMBINE_TIME="01:00:00"
DIFF_CPUS=4;    DIFF_MEM="64G";    DIFF_TIME="02:00:00"

# ==============================================================================
# Paths
# ==============================================================================

export NEC_PACKAGE_DIR="${PACKAGE_DIR}"
export NEC_INPUT_FILE="${PACKAGE_DIR}/Data/input_params_${PARAM_SET}.csv"
export NEC_DESIGN_FILE="${PACKAGE_DIR}/design/app_design_grid.rds"
export NEC_OUTPUT_DIR="${PACKAGE_DIR}/Data/app_outputs/${PARAM_SET}"
export NEC_N_SETS="${N_SETS}"
export NEC_N_TASKS="${N_TASKS}"
export NEC_SEED="${SEED}"
export R_LIBS_USER="${R_LIBRARY}"
LOG_DIR="${PACKAGE_DIR}/logs"

for f in "${NEC_INPUT_FILE}" "${NEC_DESIGN_FILE}"; do
  if [[ ! -f "${f}" ]]; then
    echo "File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "${R_LIBRARY}" "${LOG_DIR}" "${NEC_OUTPUT_DIR}/parts"

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
for (pkg in c("remotes", "cmprsk", "data.table", "R6", "survival", "truncnorm")) {
  install_if_missing(pkg)
}

# GitHub packages
install_if_missing("nhppp", "bladder-ca/nhppp")
install_if_missing("mortality", "bladder-ca/mortality")
RSCRIPT

# Always reinstall NECModel from this checkout, so the jobs run the code in
# PACKAGE_DIR rather than an older installed copy.
echo "Installing NECModel from ${PACKAGE_DIR}"
R CMD INSTALL --no-docs --no-multiarch -l "${R_LIBRARY}" "${PACKAGE_DIR}" > "${LOG_DIR}/app_install.log" 2>&1 || {
  echo "NECModel installation failed; see ${LOG_DIR}/app_install.log" >&2
  exit 1
}

# Record what is being run
GIT_COMMIT="$(git -C "${PACKAGE_DIR}" rev-parse --short HEAD 2>/dev/null || echo unknown)"
GIT_DIRTY="$( (git -C "${PACKAGE_DIR}" status --porcelain -- R DESCRIPTION NAMESPACE analysis design 2>/dev/null || true) | head -1)"
if [[ -n "${GIT_DIRTY}" ]]; then
  echo "Warning: uncommitted changes in R/, analysis/, design/ or package files; outputs will not match commit ${GIT_COMMIT}." >&2
fi

# ==============================================================================
# Submit jobs
# ==============================================================================

# Submit one step as a SLURM job and print its job ID.
# Arguments: step name, R script, CPUs, memory, time, dependency job ID (or ""),
#            array specification (optional, e.g. "1-25")
submit() {
  local name="$1" script="$2" cpus="$3" mem="$4" time="$5" deps="${6:-}" array="${7:-}"
  local extra=()
  if [[ -n "${deps}" ]]; then
    # For an array job, afterok waits for every task to succeed
    extra+=(--dependency="afterok:${deps}" --kill-on-invalid-dep=yes)
  fi
  local log_id="%j"
  if [[ -n "${array}" ]]; then
    extra+=(--array="${array}")
    log_id="%A_%a"
  fi
  sbatch --parsable \
    --job-name="nec_app_${name}" \
    --partition="${PARTITION}" \
    --qos="${QOS}" \
    --nodes=1 --ntasks=1 \
    --cpus-per-task="${cpus}" \
    --mem="${mem}" \
    --time="${time}" \
    --output="${LOG_DIR}/app_${name}_${log_id}.out" \
    --error="${LOG_DIR}/app_${name}_${log_id}.err" \
    --export=ALL \
    ${extra[@]+"${extra[@]}"} \
    --wrap="module load ${R_MODULE}; Rscript --vanilla ${PACKAGE_DIR}/analysis/app/${script}"
}

has_step() { [[ " ${STEPS} " == *" $1 "* ]]; }

JOB_SIM=""; JOB_COMBINE=""; JOB_DIFF=""

if has_step sim; then
  # Remove task files from earlier runs so combine only sees this run's output
  rm -f "${NEC_OUTPUT_DIR}"/parts/app_sets_task_*.rds
  JOB_SIM=$(submit sim 01_oscar_app_simulate.R "${SIM_CPUS}" "${SIM_MEM}" "${SIM_TIME}" "" "1-${N_TASKS}")
fi
if has_step combine; then
  JOB_COMBINE=$(submit combine 02_oscar_app_combine.R "${COMBINE_CPUS}" "${COMBINE_MEM}" "${COMBINE_TIME}" "${JOB_SIM}")
fi
if has_step diff; then
  JOB_DIFF=$(submit diff 03_oscar_app_differences.R "${DIFF_CPUS}" "${DIFF_MEM}" "${DIFF_TIME}" "${JOB_COMBINE}")
fi

cat > "${NEC_OUTPUT_DIR}/run_info.txt" <<EOF
pipeline: clinician app
parameter_set: ${PARAM_SET}
submitted: $(date '+%Y-%m-%d %H:%M:%S')
git_commit: ${GIT_COMMIT}${GIT_DIRTY:+ (with uncommitted changes)}
input_file: ${NEC_INPUT_FILE}
design_file: ${NEC_DESIGN_FILE}
parameter_sets: ${N_SETS} (seed ${SEED}) in ${N_TASKS} array tasks
steps: ${STEPS}
jobs: sim=${JOB_SIM:-none} combine=${JOB_COMBINE:-none} diff=${JOB_DIFF:-none}
EOF

echo "App pipeline, parameter set ${PARAM_SET} (commit ${GIT_COMMIT})"
echo "  sim (array 1-${N_TASKS}): ${JOB_SIM:-not submitted}"
echo "  combine:            ${JOB_COMBINE:-not submitted}${JOB_SIM:+ (after all of ${JOB_SIM})}"
echo "  diff:               ${JOB_DIFF:-not submitted}${JOB_COMBINE:+ (after ${JOB_COMBINE})}"
echo "Outputs: ${NEC_OUTPUT_DIR}"
echo "Monitor with: squeue -u \$USER"
