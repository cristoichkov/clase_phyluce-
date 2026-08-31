#!/bin/bash

#SBATCH --job-name=clean_uce
#SBATCH --partition=ripley
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"

SAMPLE_LIST="${1:-metadata/samples.txt}"

CONTAINER="${FASTP_CONTAINER:-containers/fastp_1.3.6.sif}"
RAW_DIR="${RAW_DIR:-data/raw-fastq}"
CLEAN_DIR="${CLEAN_DIR:-data/clean-fastq}"
REPORT_DIR="${REPORT_DIR:-out/fastp-report}"

if [[ ! -f "${SAMPLE_LIST}" ]]; then
    echo "ERROR: no existe la lista de muestras: ${SAMPLE_LIST}" >&2
    exit 1
fi

if [[ ! -f "${CONTAINER}" ]]; then
    echo "ERROR: no existe el contenedor de fastp: ${CONTAINER}" >&2
    exit 1
fi

mkdir -p "${CLEAN_DIR}" "${REPORT_DIR}" logs

mapfile -t SAMPLES < <(
    grep -vE '^[[:space:]]*(#|$)' "${SAMPLE_LIST}" |
    sed 's/^[[:space:]]*//;s/[[:space:]]*$//'
)

if [[ "${#SAMPLES[@]}" -eq 0 ]]; then
    echo "ERROR: la lista no contiene muestras." >&2
    exit 1
fi

echo "Muestras a limpiar: ${#SAMPLES[@]}"
printf '  %s\n' "${SAMPLES[@]}"

for ID in "${SAMPLES[@]}"
do
    R1="${RAW_DIR}/${ID}_1.fastq.gz"
    R2="${RAW_DIR}/${ID}_2.fastq.gz"

    OUTDIR="${CLEAN_DIR}/${ID}"
    OUT1="${OUTDIR}/${ID}-READ1.fastq.gz"
    OUT2="${OUTDIR}/${ID}-READ2.fastq.gz"

    HTML="${REPORT_DIR}/${ID}-fastp.html"
    JSON="${REPORT_DIR}/${ID}-fastp.json"

    echo
    echo "============================================================"
    echo "Procesando: ${ID}"
    echo "============================================================"

    if [[ ! -s "${R1}" ]]; then
        echo "ERROR: no existe R1: ${R1}" >&2
        exit 1
    fi

    if [[ ! -s "${R2}" ]]; then
        echo "ERROR: no existe R2: ${R2}" >&2
        exit 1
    fi

    mkdir -p "${OUTDIR}"

    if [[ -s "${OUT1}" && -s "${OUT2}" ]]; then
        echo "Los reads limpios ya existen. Se omite ${ID}."
        continue
    fi

    apptainer exec "${CONTAINER}" \
        fastp \
        -i "${R1}" \
        -I "${R2}" \
        -o "${OUT1}" \
        -O "${OUT2}" \
        --html="${HTML}" \
        --json="${JSON}" \
        --qualified_quality_phred=20 \
        --detect_adapter_for_pe \
        --thread="${SLURM_CPUS_PER_TASK}"

    test -s "${OUT1}"
    test -s "${OUT2}"

    echo "Limpieza terminada: ${ID}"
done

echo
echo "============================================================"
echo "Limpieza terminada para todas las muestras."
echo "============================================================"
