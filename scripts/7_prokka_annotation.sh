#!/bin/bash
#SBATCH --job-name=bacterial_assembly_annotation
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8

## Set up directories and variables
PROJECT_ROOT="/work/fauverlab/zachpella/namer_surface_ster_L3_pool_mRNA_transcript_data"
BACTERIAL_READS_DIR="${PROJECT_ROOT}/bacterial_reads_for_interproscan"

# Use unique directory names with timestamp to avoid conflicts
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
BACTERIAL_ASSEMBLY_DIR="${PROJECT_ROOT}/bacterial_assembly_${TIMESTAMP}"
ANNOTATION_OUTPUT_DIR="${PROJECT_ROOT}/bacterial_prokka_annotation_${TIMESTAMP}"

# Input bacterial reads
BACTERIAL_R1="${BACTERIAL_READS_DIR}/NaL3_surfster_mRNA.bacterial.R1.fastq.gz"
BACTERIAL_R2="${BACTERIAL_READS_DIR}/NaL3_surfster_mRNA.bacterial.R2.fastq.gz"

# Assembly output
BACTERIAL_CONTIGS="${BACTERIAL_ASSEMBLY_DIR}/NaL3_surfster_bacterial.contigs.fa"

# Define output directory and prefix for Prokka
PROKKA_OUT_DIR="${ANNOTATION_OUTPUT_DIR}/NaL3_surfster_bacterial_prokka"
PROKKA_PREFIX="NaL3_surfster_bacterial"

# Create the output directories
mkdir -p "${ANNOTATION_OUTPUT_DIR}"

echo "Using assembly directory: ${BACTERIAL_ASSEMBLY_DIR}"
echo "Using annotation directory: ${ANNOTATION_OUTPUT_DIR}"

## Check if input bacterial read files exist
if [ ! -f "${BACTERIAL_R1}" ] || [ ! -f "${BACTERIAL_R2}" ]; then
    echo "Error: Bacterial read files not found: ${BACTERIAL_R1} and ${BACTERIAL_R2}"
    exit 1
fi

## Load modules
module purge
module load megahit  # For assembly
module load prokka
module load signalp
module load perl/5.26
module load bioperl  # Add this for prokka dependencies
module load prodigal

## Step 1: Assemble bacterial reads
echo "Assembling bacterial reads..."
# Note: megahit will create the output directory, so don't create it beforehand
megahit \
    -1 "${BACTERIAL_R1}" \
    -2 "${BACTERIAL_R2}" \
    -o "${BACTERIAL_ASSEMBLY_DIR}" \
    --out-prefix NaL3_surfster_bacterial \
    --force \
    -t "${SLURM_CPUS_PER_TASK}"

# Check if assembly was successful
if [ ! -s "${BACTERIAL_CONTIGS}" ]; then
    echo "Error: Bacterial assembly failed or produced no contigs: ${BACTERIAL_CONTIGS}"
    exit 1
fi

echo "Assembly completed successfully. Found contigs at: ${BACTERIAL_CONTIGS}"

## Step 2: Run Prokka annotation on bacterial contigs
echo "Starting Prokka annotation for bacterial contigs..."
prokka \
    --kingdom Bacteria \
    --genus Bacterial_microbiome \
    --outdir "${PROKKA_OUT_DIR}" \
    --prefix "${PROKKA_PREFIX}" \
    --cpus "${SLURM_CPUS_PER_TASK}" \
    --locustag BAMIC \
    --compliant \
    --force \
    "${BACTERIAL_CONTIGS}"

## Verify output files were created
GFF_OUTPUT_FILE="${PROKKA_OUT_DIR}/${PROKKA_PREFIX}.gff"
if [[ -f "${GFF_OUTPUT_FILE}" && -s "${GFF_OUTPUT_FILE}" ]]; then
    echo "✓ Bacterial assembly and Prokka annotation completed successfully."
    echo "Assembly output: ${BACTERIAL_CONTIGS}"
    echo "Annotation files located in: ${PROKKA_OUT_DIR}"
    echo "Main annotation file: ${GFF_OUTPUT_FILE}"
else
    echo "✗ Error: Prokka annotation output file not created or is empty."
    exit 1
fi
