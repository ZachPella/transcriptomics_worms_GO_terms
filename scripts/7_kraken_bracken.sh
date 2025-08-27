#!/bin/bash
#SBATCH --job-name=kraken2_bracken_rna_contigs
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mem=50G
#SBATCH --time=15:00:00
#SBATCH --cpus-per-task=8
#SBATCH --partition=batch

## Set up directories and variables for RNA-seq project
PROJECT_ROOT="/work/fauverlab/zachpella/namer_surface_ster_L3_pool_mRNA_transcript_data"
CONTIGS_DIR="${PROJECT_ROOT}/megahit_assembly"
KRAKEN2_OUTPUT_DIR="${PROJECT_ROOT}/kraken2_output"
BRACKEN_OUTPUT_DIR="${PROJECT_ROOT}/bracken_output"

# Create output directories
mkdir -p "${KRAKEN2_OUTPUT_DIR}"
mkdir -p "${BRACKEN_OUTPUT_DIR}"

## Sample information
SAMPLE="NaL3_surfster_mRNA"
ASSEMBLED_CONTIGS_FILE="${CONTIGS_DIR}/${SAMPLE}.contigs.fa"

## Output file names
KRAKEN2_REPORT="${KRAKEN2_OUTPUT_DIR}/${SAMPLE}.kraken2_contigs.report"
KRAKEN2_OUTPUT="${KRAKEN2_OUTPUT_DIR}/${SAMPLE}.kraken2_contigs.output"

## Check if contigs file exists
if [ ! -f "${ASSEMBLED_CONTIGS_FILE}" ]; then
    echo "Error: Contigs file not found: ${ASSEMBLED_CONTIGS_FILE}"
    exit 1
fi

echo "Processing sample: ${SAMPLE}"
echo "Contigs file: ${ASSEMBLED_CONTIGS_FILE}"
echo "Kraken2 report: ${KRAKEN2_REPORT}"
echo "Kraken2 output: ${KRAKEN2_OUTPUT}"

## Load modules
module purge
module load kraken2/2.0.8-beta
module load bracken

## Run Kraken2 taxonomic classification
echo "Running Kraken2 taxonomic classification for ${SAMPLE}..."
kraken2 --db $KRAKEN2_DB \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --report "${KRAKEN2_REPORT}" \
    --output "${KRAKEN2_OUTPUT}" \
    "${ASSEMBLED_CONTIGS_FILE}"

## Check if Kraken2 completed successfully
if [ -f "${KRAKEN2_REPORT}" ] && [ -s "${KRAKEN2_REPORT}" ]; then
    echo "✓ Kraken2 classification completed successfully."

    ## Run Bracken for abundance estimation at Genus level (G)
    echo "Running Bracken abundance estimation for ${SAMPLE} at Genus level..."

    LEVEL_G="G"
    BRACKEN_REPORT_G="${BRACKEN_OUTPUT_DIR}/${SAMPLE}.bracken_contigs.${LEVEL_G}.report"
    bracken -d $KRAKEN2_DB \
        -i "${KRAKEN2_REPORT}" \
        -o "${BRACKEN_OUTPUT_DIR}/${SAMPLE}.bracken_contigs.${LEVEL_G}.output" \
        -w "${BRACKEN_REPORT_G}" \
        -l "${LEVEL_G}"

    echo "✓ Bracken abundance estimation completed at Genus level."

    ## Run Bracken for abundance estimation at Species level (S)
    echo "Running Bracken abundance estimation for ${SAMPLE} at Species level..."

    LEVEL_S="S"
    BRACKEN_REPORT_S="${BRACKEN_OUTPUT_DIR}/${SAMPLE}.bracken_contigs.${LEVEL_S}.report"
    bracken -d $KRAKEN2_DB \
        -i "${KRAKEN2_REPORT}" \
        -o "${BRACKEN_OUTPUT_DIR}/${SAMPLE}.bracken_contigs.${LEVEL_S}.output" \
        -w "${BRACKEN_REPORT_S}" \
        -l "${LEVEL_S}"

    echo "✓ Bracken abundance estimation completed at Species level."


else
    echo "✗ Error: Kraken2 report not found or empty for ${SAMPLE}. Skipping Bracken."
    exit 1
fi

## Generate summary statistics
echo "Generating summary statistics..."

# Count total contigs classified
TOTAL_CONTIGS=$(grep -c "^>" "${ASSEMBLED_CONTIGS_FILE}")
CLASSIFIED_CONTIGS=$(awk '$1 == "C"' "${KRAKEN2_OUTPUT}" | wc -l)
UNCLASSIFIED_CONTIGS=$(awk '$1 == "U"' "${KRAKEN2_OUTPUT}" | wc -l)

echo ""
echo "=== CLASSIFICATION SUMMARY ==="
echo "Sample: ${SAMPLE}"
echo "Total contigs: ${TOTAL_CONTIGS}"
echo "Classified contigs: ${CLASSIFIED_CONTIGS}"
echo "Unclassified contigs: ${UNCLASSIFIED_CONTIGS}"

if [ ${TOTAL_CONTIGS} -gt 0 ]; then
    PERCENT_CLASSIFIED=$(awk "BEGIN {printf \"%.2f\", (${CLASSIFIED_CONTIGS}/${TOTAL_CONTIGS})*100}")
    echo "Percent classified: ${PERCENT_CLASSIFIED}%"
fi

echo ""
echo "=== OUTPUT FILES CREATED ==="
echo "Kraken2 report: ${KRAKEN2_REPORT}"
echo "Kraken2 output: ${KRAKEN2_OUTPUT}"
echo "Bracken Genus report: ${BRACKEN_REPORT_G}"
echo "Bracken Species report: ${BRACKEN_REPORT_S}"

echo ""
echo "✓ Taxonomic classification and abundance estimation completed for ${SAMPLE}"
