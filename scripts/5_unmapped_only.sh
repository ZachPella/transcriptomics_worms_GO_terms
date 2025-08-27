#!/bin/bash
#SBATCH --job-name=unmapped_for_assembly
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mem=15G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1

## Set up directories and variables
PROJECT_ROOT="/work/fauverlab/zachpella/namer_surface_ster_L3_pool_mRNA_transcript_data"
BAM_INPUT_DIR="${PROJECT_ROOT}/bams_aligned"
UNMAPPED_FASTQ_OUTPUT_DIR="${PROJECT_ROOT}/unmapped_fastq"

# Create output directory
mkdir -p "${UNMAPPED_FASTQ_OUTPUT_DIR}"

# Define file names
BAM_INPUT_FILE="NaL3_surfster_mRNA.sorted.bam"
OUTPUT_R1="NaL3_surfster_mRNA.unmapped.R1.fastq.gz"
OUTPUT_R2="NaL3_surfster_mRNA.unmapped.R2.fastq.gz"
OUTPUT_SINGLETONS="NaL3_surfster_mRNA.unmapped.singleton.fastq.gz"

# Construct the full path to the current sample's BAM file
CURRENT_BAM="${BAM_INPUT_DIR}/${BAM_INPUT_FILE}"

echo "Input BAM file: ${CURRENT_BAM}"
echo "Output directory for unmapped FastQ: ${UNMAPPED_FASTQ_OUTPUT_DIR}"

## Check if input BAM file exists
if [ ! -f "${CURRENT_BAM}" ]; then
    echo "Error: Input BAM file not found: ${CURRENT_BAM}"
    exit 1
fi

## Load modules
module purge
module load samtools

# Extract unmapped paired reads and convert to FASTQ
echo "Extracting unmapped paired reads to FASTQ..."
samtools fastq -f 12 -F 256 \
    -1 "${UNMAPPED_FASTQ_OUTPUT_DIR}/${OUTPUT_R1}" \
    -2 "${UNMAPPED_FASTQ_OUTPUT_DIR}/${OUTPUT_R2}" \
    -s "${UNMAPPED_FASTQ_OUTPUT_DIR}/${OUTPUT_SINGLETONS}" \
    "${CURRENT_BAM}"

## Verify output files were created
if [[ -f "${UNMAPPED_FASTQ_OUTPUT_DIR}/${OUTPUT_R1}" && \
      -f "${UNMAPPED_FASTQ_OUTPUT_DIR}/${OUTPUT_R2}" ]]; then
    echo "✓ Successfully extracted unmapped reads to FastQ."
else
    echo "✗ Error: Unmapped FastQ files not created."
    exit 1
fi
