#!/bin/bash
#SBATCH --job-name=fastp_rna_data
#SBATCH --mail-user=zpella@unmc.edu
#SBATCH --mail-type=ALL
#SBATCH --time=0-06:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8 # Increased CPUs for faster processing of a single large sample
#SBATCH --mem=30G # Increased memory for a single large sample
#SBATCH --partition=batch

## Set up directories and variables
PROJECT_ROOT="/work/fauverlab/zachpella/namer_surface_ster_L3_pool_mRNA_transcript_data"
READSDIR="${PROJECT_ROOT}"
WORKDIR="${PROJECT_ROOT}/trimmed_reads"

# Use quotes for safety with paths
mkdir -p "${WORKDIR}"

## Set file names
READS1="NaL3_surfster_mRNA_S1_L001_R1_001.fastq.gz"
READS2="NaL3_surfster_mRNA_S1_L001_R2_001.fastq.gz"

## Set output file names
READS1_TRIMMED="NaL3_surfster_mRNA_S1_L001_R1_001_trimmed.fastq.gz"
READS2_TRIMMED="NaL3_surfster_mRNA_S1_L001_R2_001_trimmed.fastq.gz"

## Confirm that variables are properly assigned
echo "PROJECT DIRECTORY: ${PROJECT_ROOT}"
echo "FORWARD READS: ${READSDIR}/${READS1}"
echo "REVERSE READS: ${READSDIR}/${READS2}"
echo "TRIMMED FORWARD READS: ${WORKDIR}/${READS1_TRIMMED}"
echo "TRIMMED REVERSE READS: ${WORKDIR}/${READS2_TRIMMED}"
echo "Starting fastp..."

## Check if input files exist
if [ ! -f "${READSDIR}/${READS1}" ] || [ ! -f "${READSDIR}/${READS2}" ]; then
    echo "Error: Input read files not found:"
    echo "  ${READSDIR}/${READS1}"
    echo "  ${READSDIR}/${READS2}"
    exit 1
fi

## Load modules
module purge
module load fastp

## Run fastp
fastp \
    --in1 "${READSDIR}/${READS1}" \
    --in2 "${READSDIR}/${READS2}" \
    --out1 "${WORKDIR}/${READS1_TRIMMED}" \
    --out2 "${WORKDIR}/${READS2_TRIMMED}" \
    --length_required 50 \
    --thread 8 \
    --json "${WORKDIR}/NaL3_surfster_mRNA_fastp.json" \
    --html "${WORKDIR}/NaL3_surfster_mRNA_fastp.html"

## Verify output files were created
if [[ -f "${WORKDIR}/${READS1_TRIMMED}" && -f "${WORKDIR}/${READS2_TRIMMED}" ]]; then
    echo "✓ fastp completed successfully."
    echo "  Output R1: ${WORKDIR}/${READS1_TRIMMED}"
    echo "  Output R2: ${WORKDIR}/${READS2_TRIMMED}"
else
    echo "✗ Error: fastp output files not created."
    exit 1
fi
