#!/bin/bash
#SBATCH --job-name=fastqc_rna_data
#SBATCH --mail-user=zpella@unmc.edu
#SBATCH --mail-type=ALL
#SBATCH --time=0-06:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=25G
#SBATCH --partition=batch

## Set working directory and variables
PROJECT_ROOT="/work/fauverlab/zachpella/namer_surface_ster_L3_pool_mRNA_transcript_data"
WORKDIR="${PROJECT_ROOT}/trimmed_reads"
QCDIR="${PROJECT_ROOT}/fastqc_results"

mkdir -p "${QCDIR}"

## Set file names
READS1_TRIMMED="NaL3_surfster_mRNA_S1_L001_R1_001_trimmed.fastq.gz"
READS2_TRIMMED="NaL3_surfster_mRNA_S1_L001_R2_001_trimmed.fastq.gz"

## Check if input files exist
if [[ ! -f "${WORKDIR}/${READS1_TRIMMED}" || ! -f "${WORKDIR}/${READS2_TRIMMED}" ]]; then
    echo "ERROR: Trimmed files not found for ${READS1_TRIMMED} and ${READS2_TRIMMED}"
    exit 1
fi

## Load modules
module purge
module load fastqc/0.12

## Run FastQC on both trimmed read files
echo "Running FastQC on R1..."
fastqc --threads ${SLURM_CPUS_PER_TASK} --outdir="${QCDIR}" "${WORKDIR}/${READS1_TRIMMED}"

echo "Running FastQC on R2..."
fastqc --threads ${SLURM_CPUS_PER_TASK} --outdir="${QCDIR}" "${WORKDIR}/${READS2_TRIMMED}"

## Verify output files were created
if [[ -f "${QCDIR}/NaL3_surfster_mRNA_S1_L001_R1_001_trimmed_fastqc.html" && \
      -f "${QCDIR}/NaL3_surfster_mRNA_S1_L001_R2_001_trimmed_fastqc.html" ]]; then
    echo "✓ FastQC completed successfully."
else
    echo "✗ Error: FastQC output files not created."
    exit 1
fi
