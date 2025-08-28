#!/bin/bash
#SBATCH --job-name=kraken2_bracken_unmapped_reads_bacteria_filter
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mem=30G
#SBATCH --time=15:00:00
#SBATCH --cpus-per-task=8
#SBATCH --partition=batch

##
## This script performs Kraken2 and Bracken taxonomic classification on unmapped reads,
## followed by a filtering step using an external Python script.
##

# Set up directories and variables
PROJECT_ROOT="/work/fauverlab/zachpella/namer_surface_ster_L3_pool_mRNA_transcript_data"
UNMAPPED_FASTQ_DIR="${PROJECT_ROOT}/unmapped_fastq"
KRAKEN2_OUTPUT_DIR="${PROJECT_ROOT}/kraken2_output_unmapped_reads"
BRACKEN_OUTPUT_DIR="${PROJECT_ROOT}/bracken_output_unmapped_reads"
FILTERED_OUTPUT_DIR="${PROJECT_ROOT}/bracken_bacteria_filtered"

# Create all necessary output directories
mkdir -p "${KRAKEN2_OUTPUT_DIR}"
mkdir -p "${BRACKEN_OUTPUT_DIR}"
mkdir -p "${FILTERED_OUTPUT_DIR}"

## Sample information
SAMPLE="NaL3_surfster_mRNA"
# Input files are the unmapped reads from the previous step
INPUT_R1="${UNMAPPED_FASTQ_DIR}/${SAMPLE}.unmapped.R1.fastq.gz"
INPUT_R2="${UNMAPPED_FASTQ_DIR}/${SAMPLE}.unmapped.R2.fastq.gz"

## Output file names for Kraken2 and Bracken
KRAKEN2_REPORT="${KRAKEN2_OUTPUT_DIR}/${SAMPLE}.kraken2_unmapped_reads.report"
KRAKEN2_OUTPUT="${KRAKEN2_OUTPUT_DIR}/${SAMPLE}.kraken2_unmapped_reads.output"

## Check if input files exist
if [ ! -f "${INPUT_R1}" ] || [ ! -f "${INPUT_R2}" ]; then
    echo "Error: Unmapped FastQ files not found: ${INPUT_R1} and ${INPUT_R2}"
    exit 1
fi

echo "Processing sample: ${SAMPLE}"
echo "Input reads: ${INPUT_R1} and ${INPUT_R2}"

## Load modules
module purge
module load kraken2/2.0.8-beta
module load bracken

##
## STEP 1: Run Kraken2 taxonomic classification on the reads
##
echo "Running Kraken2 taxonomic classification for ${SAMPLE} on unmapped reads..."
kraken2 --db $KRAKEN2_DB \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --report "${KRAKEN2_REPORT}" \
    --output "${KRAKEN2_OUTPUT}" \
    --paired "${INPUT_R1}" "${INPUT_R2}"

## Check if Kraken2 completed successfully
if [ ! -s "${KRAKEN2_REPORT}" ]; then
    echo "✗ Error: Kraken2 report not found or empty for ${SAMPLE}. Skipping Bracken and filtering."
    exit 1
fi

echo "✓ Kraken2 classification completed successfully."

##
## STEP 2: Run Bracken for abundance estimation at Genus and Species levels
##
echo "Running Bracken abundance estimation for ${SAMPLE} at Genus and Species levels..."

LEVEL_G="G"
BRACKEN_REPORT_G="${BRACKEN_OUTPUT_DIR}/${SAMPLE}.bracken_unmapped_reads.${LEVEL_G}.report"
bracken -d $KRAKEN2_DB \
    -i "${KRAKEN2_REPORT}" \
    -o "${BRACKEN_OUTPUT_DIR}/${SAMPLE}.bracken_unmapped_reads.${LEVEL_G}.output" \
    -w "${BRACKEN_REPORT_G}" \
    -l "${LEVEL_G}"

LEVEL_S="S"
BRACKEN_REPORT_S="${BRACKEN_OUTPUT_DIR}/${SAMPLE}.bracken_unmapped_reads.${LEVEL_S}.report"
bracken -d $KRAKEN2_DB \
    -i "${KRAKEN2_REPORT}" \
    -o "${BRACKEN_OUTPUT_DIR}/${SAMPLE}.bracken_unmapped_reads.${LEVEL_S}.output" \
    -w "${BRACKEN_REPORT_S}" \
    -l "${LEVEL_S}"

echo "✓ Bracken abundance estimation completed at Genus and Species levels."

##
## STEP 3: Filter Bracken reports for bacteria only using a Python script
##
echo "Filtering Bracken reports for bacteria, then excluding specific IDs using Python script..."

# Define the path to the filter_bracken_out.py script
# IMPORTANT: You must ensure this script is accessible at this path.
FILTER_BRACKEN_SCRIPT="/work/fauverlab/zachpella/braker_run/microbiome_full_pipeline/namer_surface_ster_f14_for_nri_raw_illumina_data/scripts/filter_bracken_out.py"

# Define the NCBI Taxonomy IDs to exclude
# Added Homo sapiens (9606) to the list
EXCLUDE_IDS="2759 9605 10239 2157 9606"

# Define the output filenames for the filtered reports
FILTERED_G_OUTPUT="${FILTERED_OUTPUT_DIR}/${SAMPLE}.bracken_unmapped_reads.G.bacteria_only.output"
FILTERED_S_OUTPUT="${FILTERED_OUTPUT_DIR}/${SAMPLE}.bracken_unmapped_reads.S.bacteria_only.output"

# Run the Python script for the Genus report
python "${FILTER_BRACKEN_SCRIPT}" \
    -i "${BRACKEN_OUTPUT_DIR}/${SAMPLE}.bracken_unmapped_reads.G.output" \
    -o "${FILTERED_G_OUTPUT}" \
    --exclude ${EXCLUDE_IDS}

# Run the Python script for the Species report
python "${FILTER_BRACKEN_SCRIPT}" \
    -i "${BRACKEN_OUTPUT_DIR}/${SAMPLE}.bracken_unmapped_reads.S.output" \
    -o "${FILTERED_S_OUTPUT}" \
    --exclude ${EXCLUDE_IDS}

# Convert the filtered outputs to CSV
tr '\t' ',' < "${FILTERED_G_OUTPUT}" > "${FILTERED_G_OUTPUT%.output}.csv"
tr '\t' ',' < "${FILTERED_S_OUTPUT}" > "${FILTERED_S_OUTPUT%.output}.csv"

echo "✓ Filtering completed successfully using the Python script."

##
## FINAL SUMMARY
##
echo "Generating final summary statistics..."

# Count total unmapped reads
TOTAL_READS=$(zcat "${INPUT_R1}" | grep -c "^@")
# Count classified and unclassified reads from Kraken2 output
CLASSIFIED_READS=$(awk '$1 == "C"' "${KRAKEN2_OUTPUT}" | wc -l)
UNCLASSIFIED_READS=$(awk '$1 == "U"' "${KRAKEN2_OUTPUT}" | wc -l)

echo ""
echo "=== CLASSIFICATION SUMMARY ==="
echo "Sample: ${SAMPLE}"
echo "Total unmapped reads: ${TOTAL_READS}"
echo "Classified reads: ${CLASSIFIED_READS}"
echo "Unclassified reads: ${UNCLASSIFIED_READS}"

if [ ${TOTAL_READS} -gt 0 ]; then
    PERCENT_CLASSIFIED=$(awk "BEGIN {printf \"%.2f\", (${CLASSIFIED_READS}/${TOTAL_READS})*100}")
    echo "Percent classified: ${PERCENT_CLASSIFIED}%"
fi

echo ""
echo "=== OUTPUT FILES CREATED ==-"
echo "Kraken2 report: ${KRAKEN2_REPORT}"
echo "Bracken Genus report: ${BRACKEN_REPORT_G}"
echo "Bracken Species report: ${BRACKEN_REPORT_S}"
echo ""
echo "=== BACTERIA-ONLY FILTERED REPORTS ==="
echo "Filtered Genus report (tab-delimited): ${FILTERED_G_OUTPUT}"
echo "Filtered Genus report (CSV): ${FILTERED_G_OUTPUT%.output}.csv"
echo "Filtered Species report (tab-delimited): ${FILTERED_S_OUTPUT}"
echo "Filtered Species report (CSV): ${FILTERED_S_OUTPUT%.output}.csv"


echo ""
echo "✓ Taxonomic classification and abundance estimation completed for ${SAMPLE}"
