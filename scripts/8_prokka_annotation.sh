#!/bin/bash
#SBATCH --job-name=prokka_annotation
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8

## Set up directories and variables
PROJECT_ROOT="/work/fauverlab/zachpella/namer_surface_ster_L3_pool_mRNA_transcript_data"
ASSEMBLY_INPUT_DIR="${PROJECT_ROOT}/megahit_assembly"
ANNOTATION_OUTPUT_DIR="${PROJECT_ROOT}/prokka_annotation"

# Define input file
INPUT_CONTIGS="${ASSEMBLY_INPUT_DIR}/NaL3_surfster_mRNA.contigs.fa"

# Define output directory and prefix for Prokka
PROKKA_OUT_DIR="${ANNOTATION_OUTPUT_DIR}/NaL3_surfster_mRNA_prokka"
PROKKA_PREFIX="NaL3_surfster_microbiome"

# Create the output directory
mkdir -p "${PROKKA_OUT_DIR}"

## Check if input contig file exists and is not empty
if [ ! -s "${INPUT_CONTIGS}" ]; then
    echo "Error: Input contigs file not found or is empty: ${INPUT_CONTIGS}"
    exit 1
fi

## Load modules
module purge
# Load any modules required for prokka. Common ones are perl and prokka itself.
# Adjust the module name if it's different on your system.
module load prokka # For example, if your system has a direct 'prokka' module
module load signalp # Prokka needs this dependency to function correctly
module load perl/5.26.0 # Often required as a dependency
module load prodigal # A core dependency for prokka

## Run Prokka annotation
echo "Starting Prokka annotation for ${INPUT_CONTIGS}..."
prokka \
    --kingdom Bacteria \
    --genus Necator_americanus_microbiome \
    --outdir "${PROKKA_OUT_DIR}" \
    --prefix "${PROKKA_PREFIX}" \
    --cpus "${SLURM_CPUS_PER_TASK}" \
    --locustag NAMIC \
    --compliant \
    --force \
    "${INPUT_CONTIGS}"

## Verify output files were created
# A common file to check is the GFF3 annotation file
GFF_OUTPUT_FILE="${PROKKA_OUT_DIR}/${PROKKA_PREFIX}.gff"

if [[ -f "${GFF_OUTPUT_FILE}" && -s "${GFF_OUTPUT_FILE}" ]]; then
    echo "✓ Prokka annotation completed successfully."
    echo "Output files are located in: ${PROKKA_OUT_DIR}"
    echo "Main annotation file: ${GFF_OUTPUT_FILE}"
else
    echo "✗ Error: Prokka annotation output file not created or is empty."
    exit 1
fi
