#!/bin/bash
#SBATCH --job-name=prokka_annotation
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8

## Set up directories and variables
PROJECT_ROOT="/work/fauverlab/zachpella/namer_surface_ster_L3_pool_mRNA_transcript_data"
# The input for Prokka should be the assembly of your *bacterial* reads
BACTERIAL_ASSEMBLY_DIR="${PROJECT_ROOT}/bacterial_assembly"
ANNOTATION_OUTPUT_DIR="${PROJECT_ROOT}/prokka_annotation"

# Define input file to be the assembled contigs of *only* bacterial reads
# This assumes you have already run an assembly tool (e.g., MEGAHIT) on the filtered reads.
INPUT_CONTIGS="${BACTERIAL_ASSEMBLY_DIR}/NaL3_surfster_bacterial_contigs.fa"

# Define output directory and prefix for Prokka
PROKKA_OUT_DIR="${ANNOTATION_OUTPUT_DIR}/NaL3_surfster_prokka_bacterial"
PROKKA_PREFIX="NaL3_surfster_microbiome"

# Create the output directory
mkdir -p "${PROKKA_OUT_DIR}"

## Check if input contig file exists and is not empty
if [ ! -s "${INPUT_CONTIGS}" ]; then
    echo "Error: Input contigs file not found or is empty: ${INPUT_CONTIGS}"
    echo "Make sure to run a de novo assembly on the bacterial reads first."
    exit 1
fi

## Load modules
module purge
module load prokka
module load signalp
module load perl/5.26.0
module load prodigal

## Run Prokka annotation
echo "Starting Prokka annotation for ${INPUT_CONTIGS}..."
prokka \
    --kingdom Bacteria \
    --metagenome \
    --outdir "${PROKKA_OUT_DIR}" \
    --prefix "${PROKKA_PREFIX}" \
    --cpus "${SLURM_CPUS_PER_TASK}" \
    --locustag MICROBIOME \
    --compliant \
    --force \
    "${INPUT_CONTIGS}"

## Verify output files were created
GFF_OUTPUT_FILE="${PROKKA_OUT_DIR}/${PROKKA_PREFIX}.gff"

if [[ -f "${GFF_OUTPUT_FILE}" && -s "${GFF_OUTPUT_FILE}" ]]; then
    echo "✓ Prokka annotation completed successfully."
    echo "Output files are located in: ${PROKKA_OUT_DIR}"
    echo "Main annotation file: ${GFF_OUTPUT_FILE}"
else
    echo "✗ Error: Prokka annotation output file not created or is empty."
    exit 1
fi
