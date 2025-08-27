#!/bin/bash
#SBATCH --job-name=bwa_rna_data
#SBATCH --mail-user=zpella@unmc.edu
#SBATCH --mail-type=ALL
#SBATCH --time=6-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8 # Increased CPUs for faster BWA alignment
#SBATCH --mem=45G
#SBATCH --partition=batch

## Set working directory and variables
PROJECT_ROOT="/work/fauverlab/zachpella/namer_surface_ster_L3_pool_mRNA_transcript_data"
WORKDIR="${PROJECT_ROOT}/sam_files"
READSDIR="${PROJECT_ROOT}/trimmed_reads"
REFERENCEDIR="/work/fauverlab/zachpella/braker_run/microbiome_full_pipeline"
REFERENCE="MaSuRCA_config_purged_namericanus_withMito.short.masked.fasta"

mkdir -p "${WORKDIR}"

## Set file names
READS1_TRIMMED="NaL3_surfster_mRNA_S1_L001_R1_001_trimmed.fastq.gz"
READS2_TRIMMED="NaL3_surfster_mRNA_S1_L001_R2_001_trimmed.fastq.gz"
OUTPUT_SAM="NaL3_surfster_mRNA.sam" # Consolidated SAM filename

## Check if input files exist
if [ ! -f "${READSDIR}/${READS1_TRIMMED}" ] || [ ! -f "${READSDIR}/${READS2_TRIMMED}" ]; then
    echo "Error: Trimmed read files not found."
    exit 1
fi

## Check if reference exists
if [ ! -f "${REFERENCEDIR}/${REFERENCE}" ]; then
    echo "Error: Reference file not found: ${REFERENCEDIR}/${REFERENCE}"
    exit 1
fi

## Load modules
module purge
module load bwa

## Run BWA-MEM alignment
echo "Running BWA MEM alignment..."
bwa mem \
    -t ${SLURM_CPUS_PER_TASK} \
    -M \
    "${REFERENCEDIR}/${REFERENCE}" \
    "${READSDIR}/${READS1_TRIMMED}" \
    "${READSDIR}/${READS2_TRIMMED}" > "${WORKDIR}/${OUTPUT_SAM}"

## Verify output file was created
if [[ -f "${WORKDIR}/${OUTPUT_SAM}" && -s "${WORKDIR}/${OUTPUT_SAM}" ]]; then
    echo "✓ BWA alignment completed successfully."
else
    echo "✗ Error: SAM file not created or is empty."
    exit 1
fi
