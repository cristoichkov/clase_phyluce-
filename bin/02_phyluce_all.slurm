#!/bin/bash

#SBATCH --job-name=phyluce_all
#SBATCH --partition=ripley
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"

SAMPLE_LIST="${1:-metadata/samples.txt}"
PROBES_SOURCE="${2:-/srv/bishop/phylogenomics/data/spiders/RTA-v3-probe-combine-spider-color-DUPE-SCREENED.fasta}"

CONTAINER="${PHYLUCE_CONTAINER:-containers/phyluce_1.6.8.sif}"
CLEAN_DIR="${CLEAN_DIR:-data/clean-fastq}"

GLOBAL_DIR="${GLOBAL_DIR:-out/global}"
ASSEMBLY_DIR="${GLOBAL_DIR}/spades-assemblies"
UCE_DIR="${GLOBAL_DIR}/uce-search"
LOCI_DIR="${GLOBAL_DIR}/loci"
ALIGN_DIR="${GLOBAL_DIR}/alignments"

PROBES_DIR="data/probes"
PROBES="${PROBES_DIR}/RTA-v3-probes.fasta"

ASSEMBLY_CONF="metadata/assembly-all.conf"
TAXON_CONF="metadata/taxon-set-all.conf"

mkdir -p logs metadata "${PROBES_DIR}" "${GLOBAL_DIR}"

if [[ ! -f "${SAMPLE_LIST}" ]]; then
    echo "ERROR: no existe la lista: ${SAMPLE_LIST}" >&2
    exit 1
fi

if [[ ! -f "${CONTAINER}" ]]; then
    echo "ERROR: no existe el contenedor PHYLUCE: ${CONTAINER}" >&2
    exit 1
fi

if [[ ! -f "${PROBES_SOURCE}" ]]; then
    echo "ERROR: no existe el archivo de probes: ${PROBES_SOURCE}" >&2
    exit 1
fi

mapfile -t SAMPLES < <(
    grep -vE '^[[:space:]]*(#|$)' "${SAMPLE_LIST}" |
    sed 's/^[[:space:]]*//;s/[[:space:]]*$//'
)

N_TAXA="${#SAMPLES[@]}"

if [[ "${N_TAXA}" -lt 3 ]]; then
    echo "ERROR: para el flujo filogenómico se requieren al menos 3 muestras." >&2
    exit 1
fi

echo "Número de muestras: ${N_TAXA}"
printf '  %s\n' "${SAMPLES[@]}"

# -------------------------------------------------------------------
# Verificar reads limpios
# -------------------------------------------------------------------

for ID in "${SAMPLES[@]}"
do
    R1="${CLEAN_DIR}/${ID}/${ID}-READ1.fastq.gz"
    R2="${CLEAN_DIR}/${ID}/${ID}-READ2.fastq.gz"

    if [[ ! -s "${R1}" || ! -s "${R2}" ]]; then
        echo "ERROR: faltan reads limpios para ${ID}" >&2
        echo "Esperados:" >&2
        echo "  ${R1}" >&2
        echo "  ${R2}" >&2
        exit 1
    fi
done

# -------------------------------------------------------------------
# Crear assembly-all.conf
# -------------------------------------------------------------------

PROJECT="$(pwd)"

{
    echo "[samples]"
    for ID in "${SAMPLES[@]}"
    do
        echo "${ID}:${PROJECT}/${CLEAN_DIR}/${ID}"
    done
} > "${ASSEMBLY_CONF}"

echo
echo "Archivo generado: ${ASSEMBLY_CONF}"
cat "${ASSEMBLY_CONF}"

# -------------------------------------------------------------------
# Recordatorio de configuración de SPAdes
# -------------------------------------------------------------------

echo
echo "Configuración PHYLUCE del usuario:"
if [[ -f "${HOME}/.phyluce.conf" ]]; then
    cat "${HOME}/.phyluce.conf"
else
    echo "ADVERTENCIA: no existe ${HOME}/.phyluce.conf"
    echo "Para esta práctica se recomienda:"
    echo
    echo "[spades]"
    echo "max_memory:8"
    echo
fi

# -------------------------------------------------------------------
# Ensamblaje de todas las muestras
# -------------------------------------------------------------------

if [[ -e "${ASSEMBLY_DIR}" ]]; then
    echo "ERROR: ya existe ${ASSEMBLY_DIR}" >&2
    echo "Renómbrelo o elimínelo antes de repetir el pipeline." >&2
    exit 1
fi

echo
echo "============================================================"
echo "1. ENSAMBLAJE CON SPADES"
echo "============================================================"

apptainer exec "${CONTAINER}" \
    phyluce_assembly_assemblo_spades \
    --config "${ASSEMBLY_CONF}" \
    --output "${ASSEMBLY_DIR}" \
    --cores "${SLURM_CPUS_PER_TASK}" \
    --log-path logs

# -------------------------------------------------------------------
# Verificar contigs
# -------------------------------------------------------------------

echo
echo "Verificando contigs..."

for ID in "${SAMPLES[@]}"
do
    CONTIG="${ASSEMBLY_DIR}/contigs/${ID}.contigs.fasta"

    if [[ ! -s "${CONTIG}" ]]; then
        echo "ERROR: no se produjo el contig esperado para ${ID}" >&2
        echo "${CONTIG}" >&2
        exit 1
    fi

    N_CONTIGS=$(grep -c '^>' "${CONTIG}" || true)
    echo "${ID}: ${N_CONTIGS} contigs"
done

# -------------------------------------------------------------------
# Copiar probes
# -------------------------------------------------------------------

echo
echo "============================================================"
echo "2. PREPARAR PROBES"
echo "============================================================"

cp "${PROBES_SOURCE}" "${PROBES}"

test -s "${PROBES}"

# Verificación desde dentro de Apptainer
apptainer exec "${CONTAINER}" \
    ls -lh "${PROBES}"

echo "Primeros encabezados:"
grep '^>' "${PROBES}" | head || true

# -------------------------------------------------------------------
# Match contigs vs probes
# -------------------------------------------------------------------

if [[ -e "${UCE_DIR}" ]]; then
    echo "ERROR: ya existe ${UCE_DIR}" >&2
    exit 1
fi

echo
echo "============================================================"
echo "3. MATCH CONTIGS CONTRA PROBES"
echo "============================================================"

apptainer exec "${CONTAINER}" \
    phyluce_assembly_match_contigs_to_probes \
    --contigs "${ASSEMBLY_DIR}/contigs" \
    --probes "${PROBES}" \
    --output "${UCE_DIR}" \
    --log-path logs

test -s "${UCE_DIR}/probe.matches.sqlite"

# -------------------------------------------------------------------
# Crear taxon-set
# -------------------------------------------------------------------

{
    echo "[all]"
    printf '%s\n' "${SAMPLES[@]}"
} > "${TAXON_CONF}"

echo
echo "Archivo generado: ${TAXON_CONF}"
cat "${TAXON_CONF}"

# -------------------------------------------------------------------
# Obtener matriz muestra x locus
# -------------------------------------------------------------------

mkdir -p "${LOCI_DIR}"

echo
echo "============================================================"
echo "4. OBTENER LOCI COMPARTIDOS/INCOMPLETOS"
echo "============================================================"

apptainer exec "${CONTAINER}" \
    phyluce_assembly_get_match_counts \
    --locus-db "${UCE_DIR}/probe.matches.sqlite" \
    --taxon-list-config "${TAXON_CONF}" \
    --taxon-group all \
    --output "${LOCI_DIR}/all-uce-loci.conf" \
    --incomplete-matrix

test -s "${LOCI_DIR}/all-uce-loci.conf"

# -------------------------------------------------------------------
# Extraer FASTA monolítico
# -------------------------------------------------------------------

echo
echo "============================================================"
echo "5. EXTRAER SECUENCIAS UCE"
echo "============================================================"

apptainer exec "${CONTAINER}" \
    phyluce_assembly_get_fastas_from_match_counts \
    --contigs "${ASSEMBLY_DIR}/contigs" \
    --locus-db "${UCE_DIR}/probe.matches.sqlite" \
    --match-count-output "${LOCI_DIR}/all-uce-loci.conf" \
    --incomplete-matrix "${LOCI_DIR}/all-uce-loci.incomplete" \
    --output "${LOCI_DIR}/all-uce-loci.fasta"

test -s "${LOCI_DIR}/all-uce-loci.fasta"

N_SEQS=$(grep -c '^>' "${LOCI_DIR}/all-uce-loci.fasta" || true)
echo "Secuencias UCE extraídas: ${N_SEQS}"

# -------------------------------------------------------------------
# Alinear por locus
# -------------------------------------------------------------------

if [[ -e "${ALIGN_DIR}" ]]; then
    echo "ERROR: ya existe ${ALIGN_DIR}" >&2
    exit 1
fi

echo
echo "============================================================"
echo "6. ALINEAR CADA LOCUS CON MAFFT"
echo "============================================================"

apptainer exec "${CONTAINER}" \
    phyluce_align_seqcap_align \
    --fasta "${LOCI_DIR}/all-uce-loci.fasta" \
    --output "${ALIGN_DIR}" \
    --taxa "${N_TAXA}" \
    --aligner mafft \
    --cores "${SLURM_CPUS_PER_TASK}" \
    --incomplete-matrix \
    --output-format fasta \
    --no-trim

N_ALIGN=$(find "${ALIGN_DIR}" -type f -name '*.fasta' | wc -l)

echo
echo "============================================================"
echo "PIPELINE TERMINADO"
echo "============================================================"
echo "Muestras: ${N_TAXA}"
echo "Alineamientos FASTA: ${N_ALIGN}"
echo
echo "Resultados principales:"
echo "  ${ASSEMBLY_DIR}/contigs/"
echo "  ${UCE_DIR}/probe.matches.sqlite"
echo "  ${LOCI_DIR}/all-uce-loci.conf"
echo "  ${LOCI_DIR}/all-uce-loci.fasta"
echo "  ${ALIGN_DIR}/"
