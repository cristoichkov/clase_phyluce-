# Práctica: procesamiento de UCEs con PHYLUCE

## Objetivo

En esta práctica seguiremos paso a paso un flujo de trabajo para procesar datos de **elementos ultraconservados (UCEs)** con PHYLUCE.

La primera parte está diseñada para comprender qué ocurre con **una muestra individual**:

1. seleccionar una muestra;
2. copiar los FASTQ;
3. limpiar los reads con `fastp`;
4. ensamblar los reads con SPAdes mediante PHYLUCE;
5. identificar UCEs comparando los contigs contra los probes;
6. extraer las regiones UCE recuperadas.

Después utilizaremos **todas las muestras juntas** para mostrar la lógica filogenómica:

7. ensamblar cada muestra;
8. reunir sus contigs;
9. identificar UCEs de todas las muestras en una base común;
10. determinar qué loci están presentes en cada muestra;
11. extraer las secuencias de los loci;
12. alinear las secuencias **por locus UCE** con MAFFT.

Finalmente se incluyen dos scripts generales:

```text
01_clean_reads.slurm
02_phyluce_all.slurm
```

Estos scripts leen una lista de muestras desde:

```text
samples.txt
```

Por lo tanto, el mismo flujo puede reutilizarse con otras muestras siempre que los FASTQ sigan el patrón:

```text
MUESTRA_1.fastq.gz
MUESTRA_2.fastq.gz
```

---

# Datos utilizados

Las muestras pertenecen al proyecto:

**Ultraconserved element phylogenomics reveals novel insights into the historical biogeography of euophryine jumping spiders (Araneae: Salticidae)**

- BioProject: `PRJNA1506394`
- SRA Study: `SRP725708`

Las cuatro corridas utilizadas durante la práctica son:

| Muestra | SRA Experiment | Uso |
|---|---|---|
| `SRR40095305` | `SRX34736494` | práctica |
| `SRR40095304` | `SRX34736495` | práctica |
| `SRR40095292` | `SRX34736507` | práctica |
| `SRR40095270` | `SRX34736529` | práctica |

Los metadatos de las bibliotecas indican:

```text
Strategy:   Targeted-Capture
Source:     GENOMIC
Selection:  Hybrid Selection
Layout:     PAIRED
Platform:   ILLUMINA
Instrument: Illumina NovaSeq 6000
```

---

# Datos preparados para clase

Los FASTQ originales contienen varios millones de reads.

Para reducir el tiempo de cómputo durante la práctica se preparó un subconjunto de:

```text
500,000 reads R1
500,000 reads R2
```

por muestra.

Los archivos están disponibles en:

```text
/srv/bishop/phylogenomics/data/class-fastq/
```

Contenido:

```text
SRR40095270_1.fastq.gz
SRR40095270_2.fastq.gz
SRR40095292_1.fastq.gz
SRR40095292_2.fastq.gz
SRR40095304_1.fastq.gz
SRR40095304_2.fastq.gz
SRR40095305_1.fastq.gz
SRR40095305_2.fastq.gz
```

> El subsampling se utiliza con fines docentes. En un análisis real debe evaluarse si reducir la cobertura es apropiado, porque puede reducir el número de loci recuperados.

---

# Flujo general

```text
FASTQ paired-end
       │
       ▼
     fastp
       │
       ▼
Reads limpios
       │
       ▼
SPAdes mediante PHYLUCE
       │
       ▼
Contigs por muestra
       │
       ▼
Probes UCE
       │
       ▼
LASTZ
       │
       ▼
probe.matches.sqlite
       │
       ▼
Matriz muestra × locus
       │
       ▼
Extracción de secuencias UCE
       │
       ▼
FASTA con todas las muestras
       │
       ▼
MAFFT por locus
       │
       ▼
un alineamiento por UCE
```

---

# PARTE I. Procesamiento paso a paso de una muestra

La primera parte de la práctica se realiza con **una sola muestra** para observar claramente qué produce cada etapa.

---

# 1. Crear el proyecto

```bash
mkdir -p phyluce_uce/{data/raw-fastq,data/clean-fastq,data/probes,out,bin,metadata,logs,containers}

cd phyluce_uce
```

Estructura:

```text
phyluce_uce/
├── data/
│   ├── raw-fastq/
│   ├── clean-fastq/
│   └── probes/
├── out/
├── bin/
├── metadata/
├── logs/
└── containers/
```

---

# 2. Preparar los contenedores

Los contenedores compartidos se encuentran en:

```text
/srv/bishop/phylogenomics/contenedores/
```

Creamos enlaces simbólicos para los contenedores:

```bash
ln -s \
/srv/bishop/phylogenomics/contenedores/fastp_1.3.6.sif \
containers/fastp_1.3.6.sif
```

```bash
ln -s \
/srv/bishop/phylogenomics/contenedores/phyluce_1.6.8.sif \
containers/phyluce_1.6.8.sif
```

Comprobamos:

```bash
ls -lh containers/
```

---

# 3. Seleccionar una muestra

Cada estudiante trabajará inicialmente con una sola accesión.

Por ejemplo:

```bash
ID=SRR40095305
```

También puede utilizarse cualquiera de:

```text
SRR40095304
SRR40095292
SRR40095270
```

Comprobamos:

```bash
echo ${ID}
```

> La variable `ID` existe solamente en la terminal actual. Si abre una nueva terminal deberá definirla nuevamente.

---

# 4. Copiar los FASTQ

```bash
cp \
/srv/bishop/phylogenomics/data/class-fastq/${ID}_1.fastq.gz \
/srv/bishop/phylogenomics/data/class-fastq/${ID}_2.fastq.gz \
data/raw-fastq/
```

Comprobamos:

```bash
ls -lh data/raw-fastq/
```

Debemos tener:

```text
${ID}_1.fastq.gz
${ID}_2.fastq.gz
```

---

# 5. Revisar los FASTQ

Observar los primeros registros:

```bash
zcat data/raw-fastq/${ID}_1.fastq.gz | head
```

Un registro FASTQ contiene cuatro líneas:

```text
@identificador
SECUENCIA
+
CALIDAD
```

Contamos los reads:

```bash
zcat data/raw-fastq/${ID}_1.fastq.gz | \
awk 'END {print NR/4}'
```

```bash
zcat data/raw-fastq/${ID}_2.fastq.gz | \
awk 'END {print NR/4}'
```

Esperamos:

```text
500000
```

en ambos archivos.

---

# 6. Limpieza con fastp

El objetivo de `fastp` es producir reads limpios antes del ensamblaje.

```text
FASTQ
  │
  ▼
fastp
  │
  ▼
reads limpios
```

Creamos:

```bash
mkdir -p out/fastp-report
```

## Script `bin/fastp.slurm`

```bash
nano bin/fastp.slurm
```

Contenido:

```bash
#!/bin/bash

#SBATCH --job-name=fastp
#SBATCH --partition=ripley
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"

ID=$1

CONTAINER="containers/fastp_1.3.6.sif"
INPUT="data/raw-fastq"
OUTDIR="data/clean-fastq/${ID}"
REPORT="out/fastp-report"

mkdir -p "${OUTDIR}" "${REPORT}"

apptainer exec ${CONTAINER} \
    fastp \
    -i ${INPUT}/${ID}_1.fastq.gz \
    -I ${INPUT}/${ID}_2.fastq.gz \
    -o ${OUTDIR}/${ID}-READ1.fastq.gz \
    -O ${OUTDIR}/${ID}-READ2.fastq.gz \
    --html=${REPORT}/${ID}-fastp.html \
    --json=${REPORT}/${ID}-fastp.json \
    --qualified_quality_phred=20 \
    --detect_adapter_for_pe \
    --thread=${SLURM_CPUS_PER_TASK}
```

Ejecutamos:

```bash
sbatch bin/fastp.slurm ${ID}
```

Revisamos:

```bash
squeue -u $USER
```

Cuando termine:

```bash
ls -lh data/clean-fastq/${ID}/
```

Esperamos:

```text
${ID}-READ1.fastq.gz
${ID}-READ2.fastq.gz
```

---

# 7. ¿Por qué cambiamos READ1 y READ2?

Los FASTQ originales tienen:

```text
${ID}_1.fastq.gz
${ID}_2.fastq.gz
```

Después de `fastp` los nombramos:

```text
${ID}-READ1.fastq.gz
${ID}-READ2.fastq.gz
```

Esto permite que el wrapper de ensamblaje de PHYLUCE reconozca correctamente los mates.

---

# 8. Crear `assembly.conf`

PHYLUCE necesita conocer qué muestra debe ensamblar y dónde están sus reads.

```bash
printf "[samples]\n%s:%s/data/clean-fastq/%s\n" \
    "${ID}" "$(pwd)" "${ID}" \
    > metadata/assembly.conf
```

Revisamos:

```bash
cat metadata/assembly.conf
```

Ejemplo:

```ini
[samples]
SRR40095305:/srv/bishop/phylogenomics/phylogen20/phyluce_uce/data/clean-fastq/SRR40095305
```

---

# 9. Configurar la memoria de SPAdes

PHYLUCE 1.6.8 lee opciones de SPAdes desde:

```text
~/.phyluce.conf
```

Para los datos reducidos de esta práctica utilizamos:

```ini
[spades]
max_memory:8
```

Si el archivo no existe:

```bash
cat > ~/.phyluce.conf <<'EOF'
[spades]
max_memory:8
EOF
```

Comprobamos:

```bash
cat ~/.phyluce.conf
```

> Si el archivo ya existe, revise su contenido antes de modificarlo.

---

# 10. Ensamblar la muestra con SPAdes

El ensamblaje transforma millones de reads cortos en secuencias más largas llamadas **contigs**.

```text
reads
 │
 ├─────────────┐
 ▼             ▼
fragmentos solapantes
       │
       ▼
     SPAdes
       │
       ▼
     contigs
```

## Script `bin/spades.slurm`

```bash
nano bin/spades.slurm
```

Contenido:

```bash
#!/bin/bash

#SBATCH --job-name=phyluce_spades
#SBATCH --partition=ripley
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"

CONTAINER="containers/phyluce_1.6.8.sif"

apptainer exec ${CONTAINER} \
    phyluce_assembly_assemblo_spades \
    --config metadata/assembly.conf \
    --output out/spades-assemblies \
    --cores ${SLURM_CPUS_PER_TASK} \
    --log-path logs
```

Ejecutamos:

```bash
sbatch bin/spades.slurm
```

Revisamos:

```bash
squeue -u $USER
```

---

# 11. Revisar los contigs

Cuando termine:

```bash
ls -lh out/spades-assemblies/contigs/
```

PHYLUCE crea un FASTA de contigs asociado a la muestra:

```text
${ID}.contigs.fasta
```

Comprobamos que el archivo exista:

```bash
test -s out/spades-assemblies/${ID}_spades/contigs.fasta \
    && echo "Ensamblaje correcto" \
    || echo "ERROR: SPAdes no produjo contigs.fasta"
```

Definimos:

```bash
CONTIGS="out/spades-assemblies/contigs/${ID}.contigs.fasta"
```

Contamos:

```bash
grep -c "^>" ${CONTIGS}
```

Cada `>` corresponde a un contig ensamblado.

---

# 12. Preparar los probes UCE

Los probes disponibles para esta práctica están en:

```text
/srv/bishop/phylogenomics/data/spiders/
```

Utilizaremos:

```text
RTA-v3-probe-combine-spider-color-DUPE-SCREENED.fasta
```

## Importante: copiar, no enlazar

Copiaremos físicamente el archivo al proyecto:

```bash
cp \
/srv/bishop/phylogenomics/data/spiders/RTA-v3-probe-combine-spider-color-DUPE-SCREENED.fasta \
data/probes/RTA-v3-probes.fasta
```

Esto evita problemas con enlaces simbólicos que apuntan fuera del directorio visible para Apptainer.

Comprobamos desde el host:

```bash
ls -lh data/probes/RTA-v3-probes.fasta
```

Y desde el contenedor:

```bash
apptainer exec containers/phyluce_1.6.8.sif \
    ls -lh data/probes/RTA-v3-probes.fasta
```

---

# 13. Entender los nombres de los probes

Revisamos:

```bash
grep "^>" data/probes/RTA-v3-probes.fasta | head
```

Ejemplo:

```text
>uce-1006_p1 |design:RTA-v1,...
>uce-1006_p2 |design:RTA-v1,...
>uce-1006_p3 |design:RTA-v1,...
```

PHYLUCE utiliza por defecto:

```text
^(uce-\d+)(?:_p\d+.*)
```

Por ejemplo:

```text
uce-1006_p1
uce-1006_p2
uce-1006_p3
      │
      ▼
   uce-1006
```

Varios probes pueden pertenecer al mismo locus UCE.

---

# 14. Comparar contigs contra probes

PHYLUCE utiliza:

```text
phyluce_assembly_match_contigs_to_probes
```

y realiza las comparaciones mediante **LASTZ**.

```text
contigs de la muestra
          +
      probes UCE
          │
          ▼
        LASTZ
          │
          ▼
UCEs identificados
```

## Script `bin/uce_search.slurm`

```bash
nano bin/uce_search.slurm
```

Contenido:

```bash
#!/bin/bash

#SBATCH --job-name=uce_search
#SBATCH --partition=ripley
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"

CONTAINER="containers/phyluce_1.6.8.sif"
PROBES="data/probes/RTA-v3-probes.fasta"

apptainer exec ${CONTAINER} \
    phyluce_assembly_match_contigs_to_probes \
    --contigs out/spades-assemblies/contigs \
    --probes ${PROBES} \
    --output out/uce-search \
    --log-path logs
```

Si existe una ejecución anterior:

```bash
rm -rf out/uce-search
```

Ejecutamos:

```bash
sbatch bin/uce_search.slurm
```

---

# 15. Revisar la búsqueda de UCEs

```bash
ls -lh out/uce-search/
```

Encontraremos:

```text
probe.matches.sqlite
${ID}.contigs.lastz
```

El archivo:

```text
probe.matches.sqlite
```

es una base de datos que almacena qué UCE fue identificado en qué contig y en qué muestra.

Podemos revisar el resumen:

```bash
grep -Ri "unique\|match" logs/ | tail -n 20
```

En una ejecución con `SRR40095305` se obtuvo:

```text
SRR40095305: 3145 (11.33%) uniques of 27752 contigs,
0 dupe probe matches,
55 UCE loci removed for matching multiple contigs,
53 contigs removed for matching multiple UCE loci
```

Esto significa que PHYLUCE retuvo **3,145 asociaciones no ambiguas entre contigs y loci UCE**.

---

# 16. Extraer los UCEs de una muestra

Creamos:

```bash
printf "[all]\n%s\n" "${ID}" \
    > metadata/taxon-set.conf
```

Revisamos:

```bash
cat metadata/taxon-set.conf
```

Para una sola muestra podemos obtener los loci con:

```bash
mkdir -p out/loci/log
```

```bash
apptainer exec containers/phyluce_1.6.8.sif \
    phyluce_assembly_get_match_counts \
    --locus-db out/uce-search/probe.matches.sqlite \
    --taxon-list-config metadata/taxon-set.conf \
    --taxon-group all \
    --output out/loci/${ID}-uce-loci.conf
```

Extraemos las secuencias:

```bash
apptainer exec containers/phyluce_1.6.8.sif \
    phyluce_assembly_get_fastas_from_match_counts \
    --contigs out/spades-assemblies/contigs \
    --locus-db out/uce-search/probe.matches.sqlite \
    --match-count-output out/loci/${ID}-uce-loci.conf \
    --output out/loci/${ID}-uce-loci.fasta
```

Contamos:

```bash
grep -c "^>" out/loci/${ID}-uce-loci.fasta
```

Ejemplo para `SRR40095305`:

```text
3145
```

---

# 17. ¿Qué representa cada secuencia?

Revisamos:

```bash
grep "^>" out/loci/${ID}-uce-loci.fasta | head
```

Ejemplo:

```text
>uce-2413_SRR40095305 |uce-2413
>uce-9410_SRR40095305 |uce-9410
>uce-20003818_SRR40095305 |uce-20003818
```

Cada registro corresponde a un **locus UCE diferente** recuperado para la muestra.

Por ejemplo:

```text
uce-2413
    │
    └── secuencia recuperada en SRR40095305
```

Una secuencia recuperada puede contener el núcleo ultraconservado y regiones flanqueantes.

Estas regiones flanqueantes suelen presentar mayor variación y pueden aportar información filogenética.

---

# PARTE II. De una muestra a un análisis filogenómico

Hasta ahora trabajamos con una sola muestra para observar el proceso.

Sin embargo, para realizar un análisis filogenómico necesitamos comparar **el mismo locus entre varias muestras**.

Con cuatro muestras:

```text
SRR40095305
SRR40095304
SRR40095292
SRR40095270
```

cada una se ensambla de manera independiente:

```text
SRR40095305 → SPAdes → contigs
SRR40095304 → SPAdes → contigs
SRR40095292 → SPAdes → contigs
SRR40095270 → SPAdes → contigs
```

Después PHYLUCE trabaja sobre el directorio que contiene los contigs de todas las muestras:

```text
spades-assemblies/
└── contigs/
    ├── SRR40095305.contigs.fasta
    ├── SRR40095304.contigs.fasta
    ├── SRR40095292.contigs.fasta
    └── SRR40095270.contigs.fasta
```

Entonces ejecutamos **una sola búsqueda de probes contra el conjunto completo de contigs**.

El resultado permite construir relaciones como:

```text
uce-2413
├── SRR40095305
├── SRR40095304
├── SRR40095292
└── SRR40095270

uce-9410
├── SRR40095305
├── SRR40095304
└── SRR40095270

uce-5435
├── SRR40095305
├── SRR40095292
└── SRR40095270
```

No todos los loci tienen que estar presentes en todas las muestras.

---

# 18. Matriz completa e incompleta

Una **matriz completa** conserva únicamente loci presentes en todas las muestras.

Por ejemplo:

```text
                uce-1   uce-2   uce-3
muestra A         ✓       ✓       ✓
muestra B         ✓       ✓       ✓
muestra C         ✓       ✓       ✓
muestra D         ✓       ✓       ✓
```

Una **matriz incompleta** permite ausencia de algunos loci:

```text
                uce-1   uce-2   uce-3
muestra A         ✓       ✓       ✓
muestra B         ✓       -       ✓
muestra C         ✓       ✓       ✓
muestra D         -       ✓       ✓
```

Para un primer análisis con PHYLUCE utilizaremos:

```text
--incomplete-matrix
```

porque permite recuperar y alinear loci aunque no estén presentes en las cuatro muestras.

La completitud de la matriz puede filtrarse posteriormente.

---

# 19. Crear una lista de muestras

Para automatizar el flujo utilizaremos:

```text
samples.txt
```

Contenido:

```text
SRR40095305
SRR40095304
SRR40095292
SRR40095270
```

La idea es que el pipeline no dependa de estas cuatro muestras en particular.

Un usuario puede reemplazar el contenido por sus propias muestras:

```text
muestra_01
muestra_02
muestra_03
muestra_04
muestra_05
```

siempre que existan:

```text
muestra_01_1.fastq.gz
muestra_01_2.fastq.gz
```

etc.

También se permiten líneas vacías y líneas que comienzan con `#`.

---

# PARTE III. Pipeline global

Para análisis repetibles separaremos dos procesos:

```text
1. limpieza de reads
2. PHYLUCE
```

Esto permite revisar los reads limpios antes de comenzar el ensamblaje.

Los scripts incluidos son:

```text
01_clean_reads.slurm
02_phyluce_all.slurm
```

---

# 20. Script global de limpieza

El script:

```text
01_clean_reads.slurm
```

lee cada identificador de:

```text
samples.txt
```

y ejecuta `fastp`.

Uso:

```bash
sbatch bin/01_clean_reads.slurm metadata/samples.txt
```

El script espera encontrar:

```text
data/raw-fastq/ID_1.fastq.gz
data/raw-fastq/ID_2.fastq.gz
```

y produce:

```text
data/clean-fastq/ID/ID-READ1.fastq.gz
data/clean-fastq/ID/ID-READ2.fastq.gz
```

Los reportes quedan en:

```text
out/fastp-report/
```

Si una muestra ya tiene ambos FASTQ limpios, el script la omite.

---

# 21. Copiar las cuatro muestras de la práctica

Antes de ejecutar el script global de limpieza podemos copiar todas las muestras listadas:

```bash
while read -r ID
do
    [[ -z "${ID}" || "${ID}" =~ ^# ]] && continue

    cp \
    /srv/bishop/phylogenomics/data/class-fastq/${ID}_1.fastq.gz \
    /srv/bishop/phylogenomics/data/class-fastq/${ID}_2.fastq.gz \
    data/raw-fastq/
done < metadata/samples.txt
```

Comprobamos:

```bash
ls -lh data/raw-fastq/
```

---

# 22. Ejecutar la limpieza global

```bash
sbatch bin/01_clean_reads.slurm metadata/samples.txt
```

Revisamos:

```bash
squeue -u $USER
```

Cuando termine:

```bash
find data/clean-fastq -type f | sort
```

---

# 23. Script global de PHYLUCE

El segundo script:

```text
02_phyluce_all.slurm
```

parte de los FASTQ ya limpios.

Realiza automáticamente:

```text
samples.txt
    │
    ▼
crear assembly-all.conf
    │
    ▼
SPAdes para cada muestra
    │
    ▼
contigs de todas las muestras
    │
    ▼
copiar probes
    │
    ▼
match_contigs_to_probes
    │
    ▼
probe.matches.sqlite
    │
    ▼
crear taxon-set-all.conf
    │
    ▼
get_match_counts
    │
    ▼
get_fastas_from_match_counts
    │
    ▼
all-uce-loci.fasta
    │
    ▼
MAFFT por locus
    │
    ▼
un FASTA alineado por cada UCE
```

Uso para las muestras de la práctica:

```bash
sbatch bin/02_phyluce_all.slurm metadata/samples.txt
```

También puede especificarse otro archivo de probes:

```bash
sbatch bin/02_phyluce_all.slurm \
    metadata/samples.txt \
    /ruta/a/mis-probes.fasta
```

Si no se proporciona un segundo argumento, utiliza por defecto:

```text
/srv/bishop/phylogenomics/data/spiders/RTA-v3-probe-combine-spider-color-DUPE-SCREENED.fasta
```

---

# 24. `assembly-all.conf`

El pipeline genera automáticamente:

```text
metadata/assembly-all.conf
```

Con las cuatro muestras tendrá una estructura semejante a:

```ini
[samples]
SRR40095305:/ruta/proyecto/data/clean-fastq/SRR40095305
SRR40095304:/ruta/proyecto/data/clean-fastq/SRR40095304
SRR40095292:/ruta/proyecto/data/clean-fastq/SRR40095292
SRR40095270:/ruta/proyecto/data/clean-fastq/SRR40095270
```

PHYLUCE recorre estas muestras y ensambla cada una de manera independiente.

---

# 25. Contigs de todas las muestras

El pipeline global utiliza:

```text
out/global/spades-assemblies/
```

y produce:

```text
out/global/spades-assemblies/contigs/
├── SRR40095305.contigs.fasta
├── SRR40095304.contigs.fasta
├── SRR40095292.contigs.fasta
└── SRR40095270.contigs.fasta
```

Este directorio representa el punto donde las muestras dejan de analizarse de forma aislada y comienzan a integrarse en el análisis filogenómico.

---

# 26. Buscar UCEs en todas las muestras

El pipeline ejecuta:

```text
phyluce_assembly_match_contigs_to_probes
```

una sola vez sobre:

```text
out/global/spades-assemblies/contigs/
```

El resultado principal es:

```text
out/global/uce-search/probe.matches.sqlite
```

Esta base contiene los matches de **todas las muestras**.

---

# 27. Crear el conjunto de taxones automáticamente

A partir de `samples.txt`, el pipeline genera:

```text
metadata/taxon-set-all.conf
```

Con:

```ini
[all]
SRR40095305
SRR40095304
SRR40095292
SRR40095270
```

---

# 28. Determinar qué loci están presentes

El pipeline ejecuta:

```text
phyluce_assembly_get_match_counts
```

con:

```text
--incomplete-matrix
```

para generar:

```text
out/global/loci/all-uce-loci.conf
```

Este archivo representa la relación entre:

```text
muestras × loci UCE
```

---

# 29. Extraer todas las secuencias UCE

Después se ejecuta:

```text
phyluce_assembly_get_fastas_from_match_counts
```

para generar:

```text
out/global/loci/all-uce-loci.fasta
```

Este archivo es un **FASTA monolítico** que contiene secuencias de múltiples loci y múltiples muestras.

Por ejemplo:

```text
>uce-2413_SRR40095305
SECUENCIA...

>uce-2413_SRR40095304
SECUENCIA...

>uce-9410_SRR40095305
SECUENCIA...

>uce-9410_SRR40095270
SECUENCIA...
```

Todavía no es un único alineamiento.

Contiene muchos loci diferentes.

---

# 30. Separar y alinear por locus

PHYLUCE utiliza:

```text
phyluce_align_seqcap_align
```

para organizar las secuencias por locus y alinearlas.

Conceptualmente:

```text
all-uce-loci.fasta
        │
        ├── uce-2413
        │      ├── muestra 1
        │      ├── muestra 2
        │      ├── muestra 3
        │      └── muestra 4
        │             │
        │             ▼
        │           MAFFT
        │             │
        │             ▼
        │      uce-2413.fasta
        │
        ├── uce-9410
        │      ├── muestra 1
        │      ├── muestra 2
        │      └── muestra 4
        │             │
        │             ▼
        │           MAFFT
        │             │
        │             ▼
        │      uce-9410.fasta
        │
        └── ...
```

El script usa:

```text
--incomplete-matrix
```

porque no todos los loci necesariamente tienen datos de todas las muestras.

Además utiliza:

```text
--no-trim
```

para que esta práctica separe claramente:

```text
alineamiento
```

de:

```text
trimming / filtrado
```

Así los archivos finales representan los alineamientos producidos por MAFFT sin aplicar el trimming automático de PHYLUCE.

Si se desea utilizar el trimming automático de PHYLUCE, puede eliminarse:

```text
--no-trim
```

---

# 31. Resultados del pipeline global

Los resultados principales se encuentran en:

```text
out/global/
├── spades-assemblies/
│   └── contigs/
│       ├── muestra1.contigs.fasta
│       ├── muestra2.contigs.fasta
│       └── ...
│
├── uce-search/
│   ├── probe.matches.sqlite
│   └── *.lastz
│
├── loci/
│   ├── all-uce-loci.conf
│   ├── all-uce-loci.fasta
│   └── all-uce-loci.incomplete
│
└── alignments/
    ├── uce-....fasta
    ├── uce-....fasta
    ├── uce-....fasta
    └── ...
```

---

# 32. Revisar los alineamientos

Contamos cuántos loci fueron alineados:

```bash
find out/global/alignments \
    -type f -name "*.fasta" | wc -l
```

Podemos seleccionar uno:

```bash
find out/global/alignments \
    -type f -name "*.fasta" | head
```

Y observarlo:

```bash
head out/global/alignments/ARCHIVO.fasta
```

En un locus recuperado para las cuatro muestras esperamos encontrar hasta cuatro secuencias:

```text
uce-X
├── SRR40095305
├── SRR40095304
├── SRR40095292
└── SRR40095270
```

---

# 33. ¿Qué hemos logrado?

Al inicio teníamos:

```text
millones de reads cortos
```

Al final tenemos:

```text
un conjunto de alineamientos homólogos
```

donde cada archivo corresponde a una región UCE:

```text
uce-2413.fasta
uce-9410.fasta
uce-5435.fasta
...
```

y cada alineamiento contiene las muestras en las que se recuperó ese locus.

Este es el punto de partida para etapas posteriores como:

```text
filtrado de alineamientos
        │
        ▼
selección por completitud
        │
        ▼
concatenación o análisis por loci
        │
        ▼
inferencia filogenética
```

---

# 34. Resumen conceptual

```text
MUESTRA 1
FASTQ → fastp → SPAdes → contigs
                              │
MUESTRA 2                     │
FASTQ → fastp → SPAdes → contigs
                              │
MUESTRA 3                     ├──► PHYLUCE + probes
FASTQ → fastp → SPAdes → contigs      │
                              │        ▼
MUESTRA 4                     │   probe.matches.sqlite
FASTQ → fastp → SPAdes → contigs      │
                                       ▼
                              muestra × locus
                                       │
                                       ▼
                              secuencias UCE
                                       │
                                       ▼
                                    MAFFT
                                       │
                                       ▼
                              alineamiento por UCE
```

---

# 35. Archivos incluidos

```text
README.md
metadata/samples.txt
bin/01_clean_reads.slurm
bin/02_phyluce_all.slurm
```

Para utilizar otras muestras, en la mayoría de los casos basta con modificar:

```text
metadata/samples.txt
```

y colocar los FASTQ correspondientes en:

```text
data/raw-fastq/
```

El nombre incluido en `samples.txt` debe coincidir con el prefijo de los archivos:

```text
ID_1.fastq.gz
ID_2.fastq.gz
```

---

# Referencias

- PHYLUCE 1.6.8 documentation: https://phyluce.readthedocs.io/en/v1.6.8/
- UCE processing: https://phyluce.readthedocs.io/en/v1.6.8/uce-processing.html
