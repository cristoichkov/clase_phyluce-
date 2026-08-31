#!/bin/bash

#SBATCH --job-name=phyluce_all
#SBATCH --partition=ripley
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

# ============================================================
# 02_phyluce_all.slurm
#
# Pipeline global PHYLUCE.
#
# Parte de reads YA LIMPIOS y ejecuta:
#
#   1. assembly-all.conf
#   2. SPAdes por muestra
#   3. verificación de contigs
#   4. copia de probes
#   5. contigs vs probes con LASTZ
#   6. probe.matches.sqlite
#   7. taxon-set-all.conf
#   8. matriz muestra x locus
#   9. extracción de secuencias UCE
#  10. alineamiento por locus con MAFFT
#
# Uso:
#
#   sbatch bin/02_phyluce_all.slurm
#
# Opcionalmente:
#
#   sbatch bin/02_phyluce_all.slurm metadata/otro_config.sh
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
# Definir rutas
# ------------------------------------------------------------

ASSEMBLY_DIR="${GLOBAL_DIR}/spades-assemblies"
UCE_DIR="${GLOBAL_DIR}/uce-search"
LOCI_DIR="${GLOBAL_DIR}/loci"
ALIGN_DIR="${GLOBAL_DIR}/alignments"

ASSEMBLY_CONF="metadata/assembly-all.conf"
TAXON_CONF="metadata/taxon-set-all.conf"


# ------------------------------------------------------------
# Crear directorios
# ------------------------------------------------------------

mkdir -p \
    "${PROBES_DIR}" \
    "${GLOBAL_DIR}" \
    "${LOG_DIR}" \
    metadata


# ------------------------------------------------------------
# Verificar archivos generales
# ------------------------------------------------------------

if [[ ! -f "${PHYLUCE_CONTAINER}" ]]; then
    echo "ERROR: no existe el contenedor PHYLUCE:"
    echo "       ${PHYLUCE_CONTAINER}"
    exit 1
fi

if [[ ! -f "${PROBES_SOURCE}" ]]; then
    echo "ERROR: no existe el archivo de probes:"
    echo "       ${PROBES_SOURCE}"
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

N_TAXA="${#SAMPLES[@]}"

if [[ "${N_TAXA}" -lt 3 ]]; then
    echo "ERROR: se requieren al menos 3 muestras."
    exit 1
fi


# ------------------------------------------------------------
# Mostrar configuración
# ------------------------------------------------------------

echo "============================================================"
echo "PIPELINE GLOBAL PHYLUCE"
echo "============================================================"
echo
echo "Proyecto:           ${PROJECT}"
echo "Configuración:      ${CONFIG}"
echo "Lista de muestras:  ${SAMPLE_LIST}"
echo "PHYLUCE container:  ${PHYLUCE_CONTAINER}"
echo "Probes:             ${PROBES_SOURCE}"
echo "Reads limpios:      ${CLEAN_DIR}"
echo "Salida global:      ${GLOBAL_DIR}"
echo "Número de muestras: ${N_TAXA}"
echo

printf '  %s\n' "${SAMPLES[@]}"


# ------------------------------------------------------------
# Configurar SPAdes
# ------------------------------------------------------------
#
# PHYLUCE 1.6.8 obtiene max_memory de ~/.phyluce.conf.
#
# Solo creamos el archivo si no existe.
# ------------------------------------------------------------

if [[ ! -f "${HOME}/.phyluce.conf" ]]; then

    cat > "${HOME}/.phyluce.conf" <<EOF
[spades]
max_memory:${SPADES_MAX_MEMORY}
EOF

    echo
    echo "Se creó ${HOME}/.phyluce.conf:"
    cat "${HOME}/.phyluce.conf"

else

    echo
    echo "Se utilizará ${HOME}/.phyluce.conf existente:"
    cat "${HOME}/.phyluce.conf"

fi


# ------------------------------------------------------------
# Verificar reads limpios
# ------------------------------------------------------------

echo
echo "============================================================"
echo "0. VERIFICAR READS LIMPIOS"
echo "============================================================"

for ID in "${SAMPLES[@]}"
do
    R1="${CLEAN_DIR}/${ID}/${ID}-READ1.fastq.gz"
    R2="${CLEAN_DIR}/${ID}/${ID}-READ2.fastq.gz"

    if [[ ! -s "${R1}" ]]; then
        echo "ERROR: falta:"
        echo "       ${R1}"
        echo
        echo "Ejecute primero:"
        echo "sbatch bin/01_clean_reads.slurm"
        exit 1
    fi

    if [[ ! -s "${R2}" ]]; then
        echo "ERROR: falta:"
        echo "       ${R2}"
        echo
        echo "Ejecute primero:"
        echo "sbatch bin/01_clean_reads.slurm"
        exit 1
    fi

    echo "${ID}: OK"
done


# ------------------------------------------------------------
# Crear assembly-all.conf
# ------------------------------------------------------------

echo
echo "============================================================"
echo "1. CREAR assembly-all.conf"
echo "============================================================"

{
    echo "[samples]"

    for ID in "${SAMPLES[@]}"
    do
        echo "${ID}:${PROJECT}/${CLEAN_DIR}/${ID}"
    done

} > "${ASSEMBLY_CONF}"

cat "${ASSEMBLY_CONF}"


# ------------------------------------------------------------
# Evitar sobrescribir ensamblajes anteriores
# ------------------------------------------------------------

if [[ -e "${ASSEMBLY_DIR}" ]]; then
    echo
    echo "ERROR: ya existe:"
    echo "       ${ASSEMBLY_DIR}"
    echo
    echo "Renómbrelo o elimínelo antes de repetir el pipeline."
    exit 1
fi


# ------------------------------------------------------------
# 2. Ensamblaje
# ------------------------------------------------------------

echo
echo "============================================================"
echo "2. ENSAMBLAR CADA MUESTRA CON SPADES"
echo "============================================================"

# El proyecto se monta explícitamente para garantizar que
# PHYLUCE y SPAdes puedan acceder a reads, configuraciones
# y resultados.

apptainer exec \
    --bind "${PROJECT}:${PROJECT}" \
    "${PHYLUCE_CONTAINER}" \
    phyluce_assembly_assemblo_spades \
    --config "${ASSEMBLY_CONF}" \
    --output "${ASSEMBLY_DIR}" \
    --cores "${SLURM_CPUS_PER_TASK}" \
    --log-path "${LOG_DIR}"


# ------------------------------------------------------------
# 3. Verificar contigs
# ------------------------------------------------------------

echo
echo "============================================================"
echo "3. VERIFICAR CONTIGS"
echo "============================================================"

for ID in "${SAMPLES[@]}"
do
    CONTIG="${ASSEMBLY_DIR}/contigs/${ID}.contigs.fasta"

    if [[ ! -s "${CONTIG}" ]]; then
        echo "ERROR: no se produjo:"
        echo "       ${CONTIG}"
        exit 1
    fi

    N_CONTIGS=$(grep -c '^>' "${CONTIG}" || true)

    echo "${ID}: ${N_CONTIGS} contigs"
done


# ------------------------------------------------------------
# 4. Copiar probes
# ------------------------------------------------------------

echo
echo "============================================================"
echo "4. PREPARAR PROBES"
echo "============================================================"

# El archivo de probes está fuera del proyecto.
# Lo copiamos físicamente para evitar problemas con symlinks
# y binds hacia directorios externos.

cp "${PROBES_SOURCE}" "${PROBES}"

if [[ ! -s "${PROBES}" ]]; then
    echo "ERROR: no se pudieron copiar los probes."
    exit 1
fi

echo "Probes copiados a:"
echo "${PROBES}"

echo
echo "Verificación desde Apptainer:"

apptainer exec \
    --bind "${PROJECT}:${PROJECT}" \
    "${PHYLUCE_CONTAINER}" \
    ls -lh "${PROBES}"

echo
echo "Primeros encabezados:"

grep '^>' "${PROBES}" | head || true


# ------------------------------------------------------------
# 5. Match contigs contra probes
# ------------------------------------------------------------

echo
echo "============================================================"
echo "5. BUSCAR UCEs EN LOS CONTIGS"
echo "============================================================"

if [[ -e "${UCE_DIR}" ]]; then
    echo "ERROR: ya existe:"
    echo "       ${UCE_DIR}"
    echo
    echo "Renómbrelo o elimínelo antes de repetir."
    exit 1
fi

MATCH_ARGS=(
    phyluce_assembly_match_contigs_to_probes
    --contigs "${ASSEMBLY_DIR}/contigs"
    --probes "${PROBES}"
    --output "${UCE_DIR}"
    --log-path "${LOG_DIR}"
)

# Utilizar regex personalizado solamente si se definió.
if [[ -n "${PROBE_REGEX}" ]]; then
    MATCH_ARGS+=(--regex "${PROBE_REGEX}")
fi

apptainer exec \
    --bind "${PROJECT}:${PROJECT}" \
    "${PHYLUCE_CONTAINER}" \
    "${MATCH_ARGS[@]}"

if [[ ! -s "${UCE_DIR}/probe.matches.sqlite" ]]; then
    echo "ERROR: no se generó probe.matches.sqlite"
    exit 1
fi

echo
echo "Base de datos generada:"
ls -lh "${UCE_DIR}/probe.matches.sqlite"


# ------------------------------------------------------------
# 6. Crear taxon-set-all.conf
# ------------------------------------------------------------

echo
echo "============================================================"
echo "6. CREAR taxon-set-all.conf"
echo "============================================================"

{
    echo "[all]"
    printf '%s\n' "${SAMPLES[@]}"
} > "${TAXON_CONF}"

cat "${TAXON_CONF}"


# ------------------------------------------------------------
# 7. Matriz muestra x locus
# ------------------------------------------------------------

echo
echo "============================================================"
echo "7. OBTENER MATRIZ MUESTRA x LOCUS"
echo "============================================================"

mkdir -p "${LOCI_DIR}"

apptainer exec \
    --bind "${PROJECT}:${PROJECT}" \
    "${PHYLUCE_CONTAINER}" \
    phyluce_assembly_get_match_counts \
    --locus-db "${UCE_DIR}/probe.matches.sqlite" \
    --taxon-list-config "${TAXON_CONF}" \
    --taxon-group all \
    --output "${LOCI_DIR}/all-uce-loci.conf" \
    --incomplete-matrix

if [[ ! -s "${LOCI_DIR}/all-uce-loci.conf" ]]; then
    echo "ERROR: no se generó all-uce-loci.conf"
    exit 1
fi


# ------------------------------------------------------------
# 8. Extraer secuencias UCE
# ------------------------------------------------------------

echo
echo "============================================================"
echo "8. EXTRAER SECUENCIAS UCE"
echo "============================================================"

apptainer exec \
    --bind "${PROJECT}:${PROJECT}" \
    "${PHYLUCE_CONTAINER}" \
    phyluce_assembly_get_fastas_from_match_counts \
    --contigs "${ASSEMBLY_DIR}/contigs" \
    --locus-db "${UCE_DIR}/probe.matches.sqlite" \
    --match-count-output "${LOCI_DIR}/all-uce-loci.conf" \
    --incomplete-matrix "${LOCI_DIR}/all-uce-loci.incomplete" \
    --output "${LOCI_DIR}/all-uce-loci.fasta" \
    --log-path "${LOG_DIR}"

if [[ ! -s "${LOCI_DIR}/all-uce-loci.fasta" ]]; then
    echo "ERROR: no se generó all-uce-loci.fasta"
    exit 1
fi

N_SEQS=$(grep -c '^>' "${LOCI_DIR}/all-uce-loci.fasta" || true)

echo
echo "Secuencias UCE extraídas: ${N_SEQS}"


# ------------------------------------------------------------
# 9. Alinear por locus
# ------------------------------------------------------------

echo
echo "============================================================"
echo "9. ALINEAR CADA LOCUS CON MAFFT"
echo "============================================================"

if [[ -e "${ALIGN_DIR}" ]]; then
    echo "ERROR: ya existe:"
    echo "       ${ALIGN_DIR}"
    exit 1
fi

apptainer exec \
    --bind "${PROJECT}:${PROJECT}" \
    "${PHYLUCE_CONTAINER}" \
    phyluce_align_seqcap_align \
    --fasta "${LOCI_DIR}/all-uce-loci.fasta" \
    --output "${ALIGN_DIR}" \
    --taxa "${N_TAXA}" \
    --aligner mafft \
    --cores "${SLURM_CPUS_PER_TASK}" \
    --incomplete-matrix \
    --output-format fasta \
    --no-trim \
    --log-path "${LOG_DIR}"


# ------------------------------------------------------------
# Resumen final
# ------------------------------------------------------------

N_ALIGN=$(find "${ALIGN_DIR}" \
    -type f \
    -name "*.fasta" \
    | wc -l)

echo
echo "============================================================"
echo "PIPELINE TERMINADO"
echo "============================================================"
echo
echo "Número de muestras:       ${N_TAXA}"
echo "Secuencias UCE extraídas: ${N_SEQS}"
echo "Alineamientos FASTA:      ${N_ALIGN}"
echo
echo "Resultados principales:"
echo
echo "Contigs:"
echo "  ${ASSEMBLY_DIR}/contigs/"
echo
echo "Base UCE:"
echo "  ${UCE_DIR}/probe.matches.sqlite"
echo
echo "Matriz muestra x locus:"
echo "  ${LOCI_DIR}/all-uce-loci.conf"
echo
echo "FASTA monolítico:"
echo "  ${LOCI_DIR}/all-uce-loci.fasta"
echo
echo "Alineamientos por locus:"
echo "  ${ALIGN_DIR}/"
