#!/bin/bash
#SBATCH --job-name=megahit_assembly
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8

## Set up directories and variables
PROJECT_ROOT="/work/fauverlab/zachpella/namer_surface_ster_L3_pool_mRNA_transcript_data"
UNMAPPED_FASTQ_DIR="${PROJECT_ROOT}/unmapped_fastq"
ASSEMBLY_OUTPUT_DIR="${PROJECT_ROOT}/megahit_assembly"

# Create output directories
mkdir -p "${ASSEMBLY_OUTPUT_DIR}"

## Set input and output file names
R1_FASTQ="${UNMAPPED_FASTQ_DIR}/NaL3_surfster_mRNA.unmapped.R1.fastq.gz"
R2_FASTQ="${UNMAPPED_FASTQ_DIR}/NaL3_surfster_mRNA.unmapped.R2.fastq.gz"
ASSEMBLY_SAMPLE_DIR="${ASSEMBLY_OUTPUT_DIR}/NaL3_surfster_mRNA_assembly"
FINAL_CONTIGS="${ASSEMBLY_OUTPUT_DIR}/NaL3_surfster_mRNA.contigs.fa"

## Check if input files exist
if [ ! -f "${R1_FASTQ}" ] || [ ! -f "${R2_FASTQ}" ]; then
    echo "Error: Unmapped FastQ files not found."
    exit 1
fi

## Load modules
module purge
module load megahit/1.2

## Run Megahit assembly
echo "Starting MEGAHIT assembly..."
megahit \
    -1 "${R1_FASTQ}" \
    -2 "${R2_FASTQ}" \
    -t "${SLURM_CPUS_PER_TASK}" \
    -o "${ASSEMBLY_SAMPLE_DIR}"

## Copy the final contigs to a central location
echo "Copying final contigs..."
cp "${ASSEMBLY_SAMPLE_DIR}/final.contigs.fa" "${FINAL_CONTIGS}"

## Verify output file was created
if [[ -f "${FINAL_CONTIGS}" && -s "${FINAL_CONTIGS}" ]]; then
    echo "✓ MEGAHIT assembly completed successfully."
    echo "Final contigs are located at: ${FINAL_CONTIGS}"
else
    echo "✗ Error: Final contigs file not created or is empty."
    exit 1
fi
