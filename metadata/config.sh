#!/bin/bash

# ============================================================
# Configuración general del pipeline PHYLUCE
# ============================================================

# Lista de muestras
SAMPLE_LIST="metadata/samples.txt"

# FASTQ de entrada.
# Deben existir:
#   ID_1.fastq.gz
#   ID_2.fastq.gz
RAW_SOURCE_DIR="/srv/bishop/phylogenomics/data/class-fastq"

# Contenedores
CONTAINER_DIR="/srv/bishop/phylogenomics/contenedores"

FASTP_CONTAINER="${CONTAINER_DIR}/fastp_1.3.6.sif"
PHYLUCE_CONTAINER="${CONTAINER_DIR}/phyluce_1.6.8.sif"

# Probe set
PROBES_SOURCE="/srv/bishop/phylogenomics/data/spiders/RTA-v3-probe-combine-spider-color-DUPE-SCREENED.fasta"

# Regex de probes.
# Vacío = usar el regex por defecto de PHYLUCE.
#
# Para RTA-v3:
#   >uce-1006_p1 ...
# es compatible con:
#   ^(uce-\d+)(?:_p\d+.*)
PROBE_REGEX=""

# Directorios locales del proyecto
CLEAN_DIR="data/clean-fastq"
PROBES_DIR="data/probes"

FASTP_REPORT_DIR="out/fastp-report"
GLOBAL_DIR="out/global"
LOG_DIR="logs"

# Archivo local de probes dentro del proyecto
PROBES="${PROBES_DIR}/probes.fasta"

# fastp
FASTP_Q=20

# PHYLUCE / SPAdes
SPADES_MAX_MEMORY=8
