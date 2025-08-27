#!/bin/bash
#SBATCH --job-name=samtools_rna_data
#SBATCH --mail-user=zpella@unmc.edu
#SBATCH --mail-type=ALL
#SBATCH --time=4-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --partition=batch

## Set working directory and variables
PROJECT_ROOT="/work/fauverlab/zachpella/namer_surface_ster_L3_pool_mRNA_transcript_data"
SAMDIR="${PROJECT_ROOT}/sam_files"
WORKDIR="${PROJECT_ROOT}/bams_aligned"
REFERENCEDIR="/work/fauverlab/zachpella/braker_run/microbiome_full_pipeline"
REFERENCE="MaSuRCA_config_purged_namericanus_withMito.short.masked.fasta"
TARGETS="${REFERENCEDIR}/${REFERENCE}.bed"

mkdir -p "${WORKDIR}"

## Set file names
INPUT_SAM="NaL3_surfster_mRNA.sam"
OUTPUT_BAM="NaL3_surfster_mRNA.sorted.bam"
BAM_TMP_PREFIX="${WORKDIR}/NaL3_surfster_mRNA.tmp"

## Check if SAM file and BED file exist
if [ ! -f "${SAMDIR}/${INPUT_SAM}" ]; then
    echo "Error: SAM file not found: ${SAMDIR}/${INPUT_SAM}"
    exit 1
fi
if [ ! -f "${TARGETS}" ]; then
    echo "Error: TARGETS BED file not found: ${TARGETS}"
    exit 1
fi

## Load modules
module purge
module load samtools/1.19

## Convert SAM to BAM and sort
echo "Converting SAM to BAM and sorting..."
samtools sort -@ ${SLURM_CPUS_PER_TASK} -m 10G -o "${WORKDIR}/${OUTPUT_BAM}" -T "${BAM_TMP_PREFIX}" "${SAMDIR}/${INPUT_SAM}"

## Remove intermediate files and verify output
if [ -s "${WORKDIR}/${OUTPUT_BAM}" ]; then
    echo "✓ SAM-to-BAM conversion successful."
    rm "${SAMDIR}/${INPUT_SAM}" # Remove original SAM file to save space
else
    echo "Error: BAM file is empty or conversion failed."
    exit 1
fi

## Index the sorted BAM file
echo "Indexing BAM file..."
samtools index "${WORKDIR}/${OUTPUT_BAM}"

## Generate various statistics
echo "Generating flagstat statistics..."
samtools flagstat "${WORKDIR}/${OUTPUT_BAM}" > "${WORKDIR}/flagstats.NaL3_surfster_mRNA.out"

echo "Generating detailed samtools stats..."
samtools stats --reference "${REFERENCEDIR}/${REFERENCE}" "${WORKDIR}/${OUTPUT_BAM}" > "${WORKDIR}/stats.NaL3_surfster_mRNA.out"

echo "Calculating depth of coverage..."
samtools depth -a "${WORKDIR}/${OUTPUT_BAM}" | awk '{ total += $3; count++ } END { if(count > 0) print total/count; else print "N/A" }' > "${WORKDIR}/averageDOC.NaL3_surfster_mRNA.out"

echo "Counting mapped reads..."
samtools view -c -F 4 "${WORKDIR}/${OUTPUT_BAM}" > "${WORKDIR}/countmappedreads.NaL3_surfster_mRNA.out"

echo "✓ Processing completed successfully."
