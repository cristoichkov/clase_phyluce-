#!/bin/bash

#SBATCH --job-name=clean_uce
#SBATCH --partition=ripley
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

# ============================================================
# 01_clean_reads.slurm
#
# Limpieza de todas las muestras de metadata/samples.txt
# utilizando fastp.
#
# Uso:
#
#   sbatch bin/01_clean_reads.slurm
#
# Opcionalmente puede indicarse otro archivo de configuración:
#
#   sbatch bin/01_clean_reads.slurm metadata/otro_config.sh
# ============================================================

cd "${SLURM_SUBMIT_DIR}"

PROJECT="$(pwd)"
CONFIG="${1:-metadata/config.sh}"


# ------------------------------------------------------------
# Leer configuración
# ------------------------------------------------------------

if [[ ! -f "${CONFIG}" ]]; then
    echo "ERROR: no existe el archivo de configuración:"
    echo "       ${CONFIG}"
    exit 1
fi

source "${CONFIG}"


# ------------------------------------------------------------
# Crear directorios locales
# ------------------------------------------------------------

mkdir -p \
    "${CLEAN_DIR}" \
    "${FASTP_REPORT_DIR}" \
    "${LOG_DIR}"


# ------------------------------------------------------------
# Verificar rutas
# ------------------------------------------------------------

if [[ ! -f "${FASTP_CONTAINER}" ]]; then
    echo "ERROR: no existe el contenedor de fastp:"
    echo "       ${FASTP_CONTAINER}"
    exit 1
fi

if [[ ! -d "${RAW_SOURCE_DIR}" ]]; then
    echo "ERROR: no existe el directorio de FASTQ:"
    echo "       ${RAW_SOURCE_DIR}"
    exit 1
fi

if [[ ! -f "${SAMPLE_LIST}" ]]; then
    echo "ERROR: no existe la lista de muestras:"
    echo "       ${SAMPLE_LIST}"
    exit 1
fi


# ------------------------------------------------------------
# Leer muestras
# ------------------------------------------------------------

mapfile -t SAMPLES < <(
    grep -vE '^[[:space:]]*(#|$)' "${SAMPLE_LIST}" |
    sed 's/^[[:space:]]*//;s/[[:space:]]*$//'
)

if [[ "${#SAMPLES[@]}" -eq 0 ]]; then
    echo "ERROR: la lista de muestras está vacía."
    exit 1
fi


# ------------------------------------------------------------
# Información del job
# ------------------------------------------------------------

echo "============================================================"
echo "LIMPIEZA DE READS CON fastp"
echo "============================================================"
echo
echo "Proyecto:           ${PROJECT}"
echo "Configuración:      ${CONFIG}"
echo "Lista de muestras:  ${SAMPLE_LIST}"
echo "FASTQ de entrada:   ${RAW_SOURCE_DIR}"
echo "Contenedor fastp:   ${FASTP_CONTAINER}"
echo "Reads limpios:      ${CLEAN_DIR}"
echo "Número de muestras: ${#SAMPLES[@]}"
echo

printf '  %s\n' "${SAMPLES[@]}"


# ------------------------------------------------------------
# Procesar cada muestra
# ------------------------------------------------------------

for ID in "${SAMPLES[@]}"
do
    R1="${RAW_SOURCE_DIR}/${ID}_1.fastq.gz"
    R2="${RAW_SOURCE_DIR}/${ID}_2.fastq.gz"

    OUTDIR="${CLEAN_DIR}/${ID}"

    OUT1="${OUTDIR}/${ID}-READ1.fastq.gz"
    OUT2="${OUTDIR}/${ID}-READ2.fastq.gz"

    HTML="${FASTP_REPORT_DIR}/${ID}-fastp.html"
    JSON="${FASTP_REPORT_DIR}/${ID}-fastp.json"

    echo
    echo "============================================================"
    echo "Muestra: ${ID}"
    echo "============================================================"
    echo "R1: ${R1}"
    echo "R2: ${R2}"

    if [[ ! -s "${R1}" ]]; then
        echo
        echo "ERROR: no existe o está vacío R1:"
        echo "       ${R1}"
        exit 1
    fi

    if [[ ! -s "${R2}" ]]; then
        echo
        echo "ERROR: no existe o está vacío R2:"
        echo "       ${R2}"
        exit 1
    fi

    mkdir -p "${OUTDIR}"

    # Si ambos archivos de salida ya existen, no repetir.
    if [[ -s "${OUT1}" && -s "${OUT2}" ]]; then
        echo
        echo "Los reads limpios ya existen."
        echo "Se omite ${ID}."
        continue
    fi

    # --------------------------------------------------------
    # IMPORTANTE
    #
    # RAW_SOURCE_DIR está fuera del proyecto.
    # Se monta explícitamente dentro de Apptainer como
    # solo lectura.
    #
    # El proyecto también se monta explícitamente para
    # garantizar acceso a las rutas de salida.
    # --------------------------------------------------------

    apptainer exec \
        --bind "${RAW_SOURCE_DIR}:${RAW_SOURCE_DIR}:ro" \
        --bind "${PROJECT}:${PROJECT}" \
        "${FASTP_CONTAINER}" \
        fastp \
        -i "${R1}" \
        -I "${R2}" \
        -o "${OUT1}" \
        -O "${OUT2}" \
        --html="${HTML}" \
        --json="${JSON}" \
        --qualified_quality_phred="${FASTP_Q}" \
        --detect_adapter_for_pe \
        --thread="${SLURM_CPUS_PER_TASK}"

    # Verificar resultados
    if [[ ! -s "${OUT1}" ]]; then
        echo "ERROR: fastp no produjo:"
        echo "       ${OUT1}"
        exit 1
    fi

    if [[ ! -s "${OUT2}" ]]; then
        echo "ERROR: fastp no produjo:"
        echo "       ${OUT2}"
        exit 1
    fi

    echo
    echo "Limpieza terminada: ${ID}"
done


# ------------------------------------------------------------
# Resumen
# ------------------------------------------------------------

echo
echo "============================================================"
echo "LIMPIEZA TERMINADA"
echo "============================================================"
echo

find "${CLEAN_DIR}" \
    -type f \
    -name "*.fastq.gz" \
    | sort
