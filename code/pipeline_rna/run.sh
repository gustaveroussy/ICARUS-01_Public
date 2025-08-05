#!/usr/bin/env bash
set -e

# This adds several functions and variable in the environment
PIPELINE_PREFIX="."
# shellcheck source=/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/bash/messages.sh
source "${PIPELINE_PREFIX}/common/bash/messages.sh"
# shellcheck source=/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/bash/environment.sh
source "${PIPELINE_PREFIX}/common/bash/environment.sh"

# Define pipeline related variables
declare -x SNAKEMAKE_PROFILE_PATH="${PIPELINE_PREFIX}/common/profiles"
declare -x PIPELINE_PATH="${PIPELINE_PREFIX}"
export SNAKEMAKE_PROFILE_PATH PIPELINE_PATH

CONDA_ENV_YAML="${PIPELINE_PATH}/envs/snakemake.yaml"
CONDA_ENV_NAME="bigr_snakemake"

message INFO "Checking for conda environment: ${CONDA_ENV_NAME}"
if conda env list | grep -q -w "^${CONDA_ENV_NAME}\s"; then
    message INFO "Conda environment '${CONDA_ENV_NAME}' already exists. Skipping creation."
else
    message INFO "Conda environment '${CONDA_ENV_NAME}' not found."
    message INFO "Creating environment from file: ${CONDA_ENV_YAML}..."

    # Check if the YAML file exists before trying to create the environment
    if [ ! -f "${CONDA_ENV_YAML}" ]; then
        message INFO "Error: Conda environment definition file not found at '${CONDA_ENV_YAML}'"
        exit 1
    fi

    # Create the environment from the YAML file
    conda env create --name "${CONDA_ENV_NAME}" --file "${CONDA_ENV_YAML}" && message INFO "Conda environment '${CONDA_ENV_NAME}' created" || error_handling "${LINENO}" 1 "Failed to create Conda environment"

fi
message INFO "Conda environment check complete"

SNAKEFILE_PATH="${PIPELINE_PATH}/Snakefile"
CONFIG_PATH="${PIPELINE_PATH}/config/config.hg38.yaml"
CONDA_PATH="${PIPELINE_PATH}/conda"
WRAPPERS_PATH="https://github.com/tdayris/snakemake-wrappers/raw/Unofficial"
SNAKE_ARGS=("--wrapper-prefix" "${WRAPPERS_PATH}/" "--conda-prefix" "${CONDA_PATH}")
PROFILE="slurm"
SUMMARY=""
GRAPH=""

while [ "$#" -gt 0 ]; do
  case "${1}" in
    -p|--profile) PROFILE="${2}"; shift 2;;
    --summary) SUMMARY="${2}"; shift 2;;
    --rulegraph|--dag) GRAPH="${2}"; shift 2;;
    hg38|HG38|GRCh38) CONFIG_PATH="${PIPELINE_PATH}/config.hg38.yaml"; shift;;
    DESeq2|deseq2|DGE|dge) SNAKE_ARGS+=("--until deseq2_results"); shift;;
    salmon|Salmon|quant) SNAKE_ARGS+=("--until salmon_quant_results"); shift;;
    fusions|fusion) SNAKE_ARGS+=("--until star_fusion_results"); shift;;
    qc|QC) SNAKE_ARGS+=("--until quality_control_results"); shift;;
    immu|deconv) SNAKE_ARGS+=("--until immunedeconv_results"); shift;;
    *) SNAKE_ARGS+=("${1}"); shift;;
  esac
done
message INFO "Environment loaded"

if [ ! -f "config.yaml" ]; then
  rsync -cv "${CONFIG_PATH}" "config.yaml"
else
  message INFO "Config file already provided"
fi

# Run pipeline
message CMD "conda_activate ${CONDA_ENV_NAME}"
conda_activate "${CONDA_ENV_NAME}" && message INFO "Conda loaded" || error_handling "${LINENO}" 1 "Could not activate conda environment"

BASE_CMD="snakemake -s ${SNAKEFILE_PATH} --profile ${SNAKEMAKE_PROFILE_PATH}/${PROFILE} ${SNAKE_ARGS[*]}"
if [ "${SUMMARY}" != "" ]; then
  SUMMARY_CMD="${BASE_CMD} --summary > ${SUMMARY}"
  message CMD "${SUMMARY_CMD}"
  eval SUMMARY_CMD
elif [ "${GRAPH}" != "" ]; then
  RULEGRAPH_CMD="${BASE_CMD} --rulegraph | dot -Tpng > ${GRAPH}"
  message CMD "${RULEGRAPH_CMD}"
  eval ${SUMMARY_CMD}
else
  message CMD "${BASE_CMD}"
  eval ${BASE_CMD}
fi
