#!/bin/bash
#SBATCH --job-name=rna_stats_summary
#SBATCH --mail-user=zpella@unmc.edu
#SBATCH --mail-type=ALL
#SBATCH --time=0-01:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --partition=batch

## Set working directories to match your RNA-seq project paths
PROJECT_ROOT="/work/fauverlab/zachpella/namer_surface_ster_L3_pool_mRNA_transcript_data"
TRIMMED_DIR="${PROJECT_ROOT}/trimmed_reads"
BAM_DIR="${PROJECT_ROOT}/bams_aligned"
UNMAPPED_DIR="${PROJECT_ROOT}/unmapped_fastq"
OUTDIR="${PROJECT_ROOT}"

## Set reference file
REFERENCEDIR="/work/fauverlab/zachpella/braker_run/microbiome_full_pipeline"
REFERENCE="MaSuRCA_config_purged_namericanus_withMito.short.masked.fasta"

## Load samtools module
module load samtools/1.19

## Sample name
SAMPLE="NaL3_surfster_mRNA"

echo "Generating comprehensive RNA-seq statistics for ${SAMPLE}..."

## Create a header for the comprehensive summary file
echo "Sample,Raw_Reads_R1,Raw_Reads_R2,Total_Raw_Reads,Trimmed_Reads_R1,Trimmed_Reads_R2,Total_Trimmed_Reads,Percent_Reads_Kept,Total_Aligned_Reads,Primary_Reads,Primary_Mapped,Perc_Primary_Mapped,Properly_Paired,Perc_Properly_Paired,Mean_Coverage,Total_Unmapped_Reads,Perc_Unmapped" > ${OUTDIR}/rna_seq_stats_summary.csv

## Get raw read counts from original FASTQ files
echo "Extracting raw read counts..."
RAW_R1_FILE="${PROJECT_ROOT}/NaL3_surfster_mRNA_S1_L001_R1_001.fastq.gz"
RAW_R2_FILE="${PROJECT_ROOT}/NaL3_surfster_mRNA_S1_L001_R2_001.fastq.gz"

if [ -f "${RAW_R1_FILE}" ]; then
    RAW_READS_R1=$(zcat "${RAW_R1_FILE}" | wc -l | awk '{print $1/4}')
else
    RAW_READS_R1="NA"
fi

if [ -f "${RAW_R2_FILE}" ]; then
    RAW_READS_R2=$(zcat "${RAW_R2_FILE}" | wc -l | awk '{print $1/4}')
else
    RAW_READS_R2="NA"
fi

# Calculate total raw reads
if [ "$RAW_READS_R1" != "NA" ] && [ "$RAW_READS_R2" != "NA" ]; then
    TOTAL_RAW_READS=$((RAW_READS_R1 + RAW_READS_R2))
else
    TOTAL_RAW_READS="NA"
fi

## Get trimmed read counts
echo "Extracting trimmed read counts..."
TRIMMED_R1_FILE="${TRIMMED_DIR}/NaL3_surfster_mRNA_S1_L001_R1_001_trimmed.fastq.gz"
TRIMMED_R2_FILE="${TRIMMED_DIR}/NaL3_surfster_mRNA_S1_L001_R2_001_trimmed.fastq.gz"

if [ -f "${TRIMMED_R1_FILE}" ]; then
    TRIMMED_READS_R1=$(zcat "${TRIMMED_R1_FILE}" | wc -l | awk '{print $1/4}')
else
    TRIMMED_READS_R1="NA"
fi

if [ -f "${TRIMMED_R2_FILE}" ]; then
    TRIMMED_READS_R2=$(zcat "${TRIMMED_R2_FILE}" | wc -l | awk '{print $1/4}')
else
    TRIMMED_READS_R2="NA"
fi

# Calculate total trimmed reads and percentage kept
if [ "$TRIMMED_READS_R1" != "NA" ] && [ "$TRIMMED_READS_R2" != "NA" ]; then
    TOTAL_TRIMMED_READS=$((TRIMMED_READS_R1 + TRIMMED_READS_R2))
    if [ "$TOTAL_RAW_READS" != "NA" ] && [ "$TOTAL_RAW_READS" -ne 0 ]; then
        PERCENT_READS_KEPT=$(awk "BEGIN {printf \"%.2f\", ($TOTAL_TRIMMED_READS/$TOTAL_RAW_READS)*100}")
    else
        PERCENT_READS_KEPT="NA"
    fi
else
    TOTAL_TRIMMED_READS="NA"
    PERCENT_READS_KEPT="NA"
fi

## Get alignment statistics from flagstat file
echo "Extracting alignment statistics..."
FLAGSTAT_FILE="${BAM_DIR}/flagstats.${SAMPLE}.out"

if [ -f "${FLAGSTAT_FILE}" ]; then
    # Extract alignment metrics
    TOTAL_ALIGNED=$(grep "in total" "${FLAGSTAT_FILE}" | head -1 | awk '{print $1}')
    PRIMARY=$(grep "^[0-9]* + [0-9]* primary$" "${FLAGSTAT_FILE}" | awk '{print $1}')

    # Get primary mapped stats specifically
    PRIMARY_MAPPED_LINE=$(grep "primary mapped" "${FLAGSTAT_FILE}")
    PRIMARY_MAPPED=$(echo "$PRIMARY_MAPPED_LINE" | awk '{print $1}')
    PERC_PRIMARY_MAPPED=$(echo "$PRIMARY_MAPPED_LINE" | awk -F'[()]' '{print $2}' | awk '{print $1}')

    # Get properly paired stats
    PAIRED_LINE=$(grep "properly paired" "${FLAGSTAT_FILE}")
    PROPERLY_PAIRED=$(echo "$PAIRED_LINE" | awk '{print $1}')
    PERC_PROPERLY_PAIRED=$(echo "$PAIRED_LINE" | awk -F'[()]' '{print $2}' | awk '{print $1}')

    # Handle empty values
    if [ -z "$PERC_PRIMARY_MAPPED" ]; then PERC_PRIMARY_MAPPED="0"; fi
    if [ -z "$PERC_PROPERLY_PAIRED" ]; then PERC_PROPERLY_PAIRED="0"; fi
    if [ -z "$PROPERLY_PAIRED" ]; then PROPERLY_PAIRED="0"; fi
    if [ -z "$PRIMARY_MAPPED" ]; then PRIMARY_MAPPED="0"; fi
else
    TOTAL_ALIGNED="NA"
    PRIMARY="NA"
    PRIMARY_MAPPED="NA"
    PERC_PRIMARY_MAPPED="NA"
    PROPERLY_PAIRED="NA"
    PERC_PROPERLY_PAIRED="NA"
fi

## Get mean coverage from averageDOC file
echo "Extracting coverage statistics..."
AVGDOC_FILE="${BAM_DIR}/averageDOC.${SAMPLE}.out"
if [ -f "${AVGDOC_FILE}" ]; then
    MEAN_COVERAGE=$(cat "${AVGDOC_FILE}")
    # Check if we got a valid number and format it
    if [ -n "$MEAN_COVERAGE" ] && [[ "$MEAN_COVERAGE" =~ ^[0-9]+\.?[0-9]*$ ]]; then
        MEAN_COVERAGE=$(awk "BEGIN {printf \"%.2f\", ${MEAN_COVERAGE}}")
    else
        MEAN_COVERAGE="NA"
    fi
else
    MEAN_COVERAGE="NA"
fi

## Get unmapped read counts
echo "Extracting unmapped read counts..."
UNMAPPED_R1_FILE="${UNMAPPED_DIR}/${SAMPLE}.unmapped.R1.fastq.gz"
UNMAPPED_R2_FILE="${UNMAPPED_DIR}/${SAMPLE}.unmapped.R2.fastq.gz"
UNMAPPED_SINGLETON_FILE="${UNMAPPED_DIR}/${SAMPLE}.unmapped.singleton.fastq.gz"

UNMAPPED_R1=0
UNMAPPED_R2=0
UNMAPPED_SINGLETON=0

if [ -f "${UNMAPPED_R1_FILE}" ]; then
    UNMAPPED_R1=$(zcat "${UNMAPPED_R1_FILE}" | wc -l | awk '{print $1/4}')
fi

if [ -f "${UNMAPPED_R2_FILE}" ]; then
    UNMAPPED_R2=$(zcat "${UNMAPPED_R2_FILE}" | wc -l | awk '{print $1/4}')
fi

if [ -f "${UNMAPPED_SINGLETON_FILE}" ]; then
    UNMAPPED_SINGLETON=$(zcat "${UNMAPPED_SINGLETON_FILE}" | wc -l | awk '{print $1/4}')
fi

TOTAL_UNMAPPED_READS=$((UNMAPPED_R1 + UNMAPPED_R2 + UNMAPPED_SINGLETON))

# Calculate percentage unmapped
if [ "$TOTAL_TRIMMED_READS" != "NA" ] && [ "$TOTAL_TRIMMED_READS" -ne 0 ]; then
    PERC_UNMAPPED=$(awk "BEGIN {printf \"%.2f\", ($TOTAL_UNMAPPED_READS/$TOTAL_TRIMMED_READS)*100}")
else
    PERC_UNMAPPED="NA"
fi

## Output comprehensive summary
echo "Creating summary output..."
echo "$SAMPLE,$RAW_READS_R1,$RAW_READS_R2,$TOTAL_RAW_READS,$TRIMMED_READS_R1,$TRIMMED_READS_R2,$TOTAL_TRIMMED_READS,$PERCENT_READS_KEPT,$TOTAL_ALIGNED,$PRIMARY,$PRIMARY_MAPPED,$PERC_PRIMARY_MAPPED,$PROPERLY_PAIRED,$PERC_PROPERLY_PAIRED,$MEAN_COVERAGE,$TOTAL_UNMAPPED_READS,$PERC_UNMAPPED" >> ${OUTDIR}/rna_seq_stats_summary.csv

echo "RNA-seq statistics summary created: ${OUTDIR}/rna_seq_stats_summary.csv"

## Create a tab-separated version for command line viewing
sed 's/,/\t/g' ${OUTDIR}/rna_seq_stats_summary.csv > ${OUTDIR}/rna_seq_stats_summary.tsv
echo "Tab-separated version created: ${OUTDIR}/rna_seq_stats_summary.tsv"

## Create a human-readable report
echo "Creating detailed report..."
cat > ${OUTDIR}/rna_seq_detailed_report.txt << EOF
================================================================================
                    RNA-seq Processing Summary Report
                         Sample: ${SAMPLE}
                    Generated: $(date)
================================================================================

RAW DATA STATISTICS:
--------------------
Raw reads (R1):                    ${RAW_READS_R1}
Raw reads (R2):                    ${RAW_READS_R2}
Total raw reads:                   ${TOTAL_RAW_READS}

QUALITY TRIMMING (fastp):
-------------------------
Trimmed reads (R1):                ${TRIMMED_READS_R1}
Trimmed reads (R2):                ${TRIMMED_READS_R2}
Total trimmed reads:               ${TOTAL_TRIMMED_READS}
Percentage of reads kept:          ${PERCENT_READS_KEPT}%

ALIGNMENT STATISTICS (BWA-MEM):
-------------------------------
Total aligned reads:               ${TOTAL_ALIGNED}
Primary reads:                     ${PRIMARY}
Primary mapped reads:              ${PRIMARY_MAPPED} (${PERC_PRIMARY_MAPPED})
Properly paired reads:             ${PROPERLY_PAIRED} (${PERC_PROPERLY_PAIRED})
Mean coverage depth:               ${MEAN_COVERAGE}x

UNMAPPED READS:
---------------
Total unmapped reads:              ${TOTAL_UNMAPPED_READS}
Percentage unmapped:               ${PERC_UNMAPPED}%
- Unmapped R1 reads:               ${UNMAPPED_R1}
- Unmapped R2 reads:               ${UNMAPPED_R2}
- Unmapped singleton reads:        ${UNMAPPED_SINGLETON}

================================================================================
Notes:
- Unmapped reads have been extracted for potential downstream assembly
- All percentages are calculated relative to appropriate denominators
- Coverage statistics are based on the reference genome
================================================================================
EOF

echo "Detailed report created: ${OUTDIR}/rna_seq_detailed_report.txt"

echo "All statistics files have been generated successfully!"
echo ""
echo "Output files:"
echo "  - ${OUTDIR}/rna_seq_stats_summary.csv (machine-readable)"
echo "  - ${OUTDIR}/rna_seq_stats_summary.tsv (tab-separated)"
echo "  - ${OUTDIR}/rna_seq_detailed_report.txt (human-readable)"
