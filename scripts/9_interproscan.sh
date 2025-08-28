#!/bin/bash
#SBATCH --job-name=interproscan_bacterial_annotation
#SBATCH --output=/work/fauverlab/zachpella/namer_surface_ster_L3_pool_mRNA_transcript_data/prokka_annotation/NaL3_surfster_prokka_bacterial/interproscan_output/interproscan_%j.out
#SBATCH --error=/work/fauverlab/zachpella/namer_surface_ster_L3_pool_mRNA_transcript_data/prokka_annotation/NaL3_surfster_prokka_bacterial/interproscan_output/interproscan_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=20G

##
## This script runs InterProScan on the protein sequences predicted by Prokka.
## It is designed to follow the updated Prokka script.
##

# Set up directories and variables based on the previous Prokka script
PROJECT_ROOT="/work/fauverlab/zachpella/namer_surface_ster_L3_pool_mRNA_transcript_data"
# This directory should match the --outdir from the Prokka script
PROKKA_OUT_DIR="${PROJECT_ROOT}/prokka_annotation/NaL3_surfster_prokka_bacterial"
# This prefix should match the --prefix from the Prokka script
PROKKA_PREFIX="NaL3_surfster_microbiome"

# InterProScan installation directory.
INTERPROSCAN_DIR="/work/fauverlab/zachpella/braker_run/braker_output_removed_dups_all_ours/my_interproscan/interproscan-5.73-104.0"

# Define input file
# The .faa file contains the predicted amino acid sequences
INPUT_FILE="${PROKKA_OUT_DIR}/${PROKKA_PREFIX}.faa"

# Define output and temporary directories
OUTPUT_DIR="${PROKKA_OUT_DIR}/interproscan_output"
TEMP_DIR="${PROKKA_OUT_DIR}/interproscan_temp"

# Create directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${TEMP_DIR}"

# Check if input file exists
if [ ! -f "${INPUT_FILE}" ]; then
    echo "Error: Prokka .faa file not found: ${INPUT_FILE}"
    echo "Please ensure the Prokka script completed successfully and the file exists."
    exit 1
fi

echo "Processing input file: ${INPUT_FILE}"

# Create a temporary input file with stop codons removed, as some databases
# within InterProScan do not handle them correctly.
MODIFIED_INPUT="${TEMP_DIR}/$(basename ${INPUT_FILE} .faa)_nostop.faa"
sed 's/\*$//' "${INPUT_FILE}" > "${MODIFIED_INPUT}"

# Load modules
module purge
module load interproscan
module load java
module load python/3.8

# Run InterProScan with full paths and fixed CPU allocation
echo "Starting InterProScan..."
"${INTERPROSCAN_DIR}/interproscan.sh" \
    -i "${MODIFIED_INPUT}" \
    -o "${OUTPUT_DIR}/interproscan_results" \
    -f TSV,XML,GFF3 \
    -T "${TEMP_DIR}" \
    -cpu "${SLURM_CPUS_PER_TASK}" \
    --goterms \
    --pathways

# Check if the job completed successfully
if [ $? -eq 0 ]; then
    echo "✓ InterProScan completed successfully."
else
    echo "✗ Error: InterProScan failed. Please check the log file for details."
    exit 1
fi

echo ""
echo "=== OUTPUT FILES CREATED ==="
echo "Results can be found in: ${OUTPUT_DIR}"
echo "InterProScan TSV report: ${OUTPUT_DIR}/interproscan_results.tsv"
echo "InterProScan XML report: ${OUTPUT_DIR}/interproscan_results.xml"
echo "InterProScan GFF3 report: ${OUTPUT_DIR}/interproscan_results.gff3"
echo ""
echo "✓ Script completed."
