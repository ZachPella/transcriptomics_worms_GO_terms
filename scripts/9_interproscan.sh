#!/bin/bash
#SBATCH --job-name=interproscan_prokka
#SBATCH --output=/work/fauverlab/zachpella/namer_surface_ster_L3_pool_mRNA_transcript_data/prokka_annotation/NaL3_surfster_mRNA_prokka/interproscan_output/interproscan_%j.out
#SBATCH --error=/work/fauverlab/zachpella/namer_surface_ster_L3_pool_mRNA_transcript_data/prokka_annotation/NaL3_surfster_mRNA_prokka/interproscan_output/interproscan_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=20G

# The path has been corrected to include the 'prokka_annotation' subdirectory.
PROKKA_OUTPUT_DIR="/work/fauverlab/zachpella/namer_surface_ster_L3_pool_mRNA_transcript_data/prokka_annotation/NaL3_surfster_mRNA_prokka"

# Name of the input file from Prokka (the amino acid sequences).
INPUT_FILE_NAME="NaL3_surfster_microbiome.faa"

# InterProScan installation directory.
INTERPROSCAN_DIR="/work/fauverlab/zachpella/braker_run/braker_output_removed_dups_all_ours/my_interproscan/interproscan-5.73-104.0"

# Set paths based on configuration
INPUT_FILE="${PROKKA_OUTPUT_DIR}/${INPUT_FILE_NAME}"
OUTPUT_DIR="${PROKKA_OUTPUT_DIR}/interproscan_output"
TEMP_DIR="${PROKKA_OUTPUT_DIR}/interproscan_temp"

# Create directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${TEMP_DIR}"

# Check if input file exists
if [ ! -f "${INPUT_FILE}" ]; then
    echo "Error: Input file not found: ${INPUT_FILE}"
    exit 1
fi

echo "Processing input file: ${INPUT_FILE}"

# Create modified input file with stop codons removed, as some databases
# within InterProScan do not handle them correctly.
MODIFIED_INPUT="${OUTPUT_DIR}/$(basename ${INPUT_FILE_NAME} .faa)_nostop.faa"
cat "${INPUT_FILE}" | sed 's/\*$//' > "${MODIFIED_INPUT}"

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
    -f TSV \
    -T "${TEMP_DIR}" \
    -cpu "${SLURM_CPUS_PER_TASK}" \
    -verbose \
    --goterms \
    --pathways

# Check if the job completed successfully
if [ $? -eq 0 ]; then
    echo "✓ InterProScan completed successfully."
else
    echo "✗ Error: InterProScan failed."
    # Print the last few lines of the log file if it exists
    LOG_FILE="${INTERPROSCAN_DIR}/interproscan.log"
    if [ -f "${LOG_FILE}" ]; then
        echo "Last 50 lines of log file:"
        tail -50 "${LOG_FILE}"
    fi
    exit 1
fi

echo ""
echo "=== OUTPUT FILES CREATED ==="
echo "Results can be found in: ${OUTPUT_DIR}"
echo "InterProScan TSV report: ${OUTPUT_DIR}/interproscan_results.tsv"
echo ""
echo "✓ Script completed."
