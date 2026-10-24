#!/bin/bash

# Usage: ./somatic_vc_pipeline_fixed.sh <tumor_id> <normal_id>
# Example: ./somatic_vc_pipeline_fixed.sh SRR1653287 SRR1656612

set -euo pipefail

#=======================================
# Inputs (taken from command line)
#=======================================
TUMOR="${1:?Usage: $0 <tumor_id> <normal_id>}"
NORMAL="${2:?Usage: $0 <tumor_id> <normal_id>}"

#=======================================
# Defining paths
#=======================================
WORKDIR="/home/sayoni/WORKDIR/Samples/Tumor"
REF="/home/sayoni/WORKDIR/Samples/Tumor/REF/Homo_sapiens_assembly38.fasta"
DBSNP="/home/sayoni/WORKDIR/DBSNP/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
INDELS="/home/sayoni/WORKDIR/INDELS/Homo_sapiens_assembly38.known_indels.vcf.gz"
MILLS="/home/sayoni/WORKDIR/MILLS/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
FUNCOTATOR="/home/sayoni/WORKDIR/funcotator_dataSources.v1.8.hg38.20230908g"
PANEL="/home/sayoni/WORKDIR/targets.bed"
ADAPTERS="$WORKDIR/TruSeq3-PE.fa"

# Optional but recommended for Mutect2 (leave unset or missing to skip)
GERMLINE_RES="/home/sayoni/WORKDIR/GNOMAD/af-only-gnomad.hg38.vcf.gz"
PON="/home/sayoni/WORKDIR/PON/pon.vcf.gz"
COMMON_SNPS="/home/sayoni/WORKDIR/GNOMAD/small_exac_common_3.hg38.vcf.gz"

THREADS=8

cd "$WORKDIR"

#=======================================
# FASTQ naming convention
#=======================================
FASTQ_TUMOR1="${TUMOR}_1.fastq.gz"
FASTQ_TUMOR2="${TUMOR}_2.fastq.gz"
FASTQ_NORMAL1="${NORMAL}_1.fastq.gz"
FASTQ_NORMAL2="${NORMAL}_2.fastq.gz"

#=======================================
# Output structure
#=======================================
OUTPUT="$WORKDIR/outputs"
mkdir -p "$OUTPUT"/{qc_reports,trimmed_reads,bam_files,vcf_files,logs}

#=======================================
# Check tools
#=======================================
echo "Checking required WES tools..."
REQUIRED_TOOLS=("fastqc" "trimmomatic" "bwa" "samtools" "gatk" "tabix")

for tool in "${REQUIRED_TOOLS[@]}"; do
    if command -v "$tool" &> /dev/null; then
        echo "[OK] $tool found"
    else
        echo "[ERROR] $tool not found. Please install it or add it to PATH."
        exit 1
    fi
done
echo "All required tools are available."

#=======================================
# Check input files and reference indexes
#=======================================
echo "Checking input files..."
for f in "$FASTQ_TUMOR1" "$FASTQ_TUMOR2" "$FASTQ_NORMAL1" "$FASTQ_NORMAL2" \
         "$REF" "$DBSNP" "$INDELS" "$MILLS" "$PANEL" "$ADAPTERS"; do
    [ -f "$f" ] || { echo "[ERROR] Missing file: $f"; exit 1; }
done
[ -d "$FUNCOTATOR" ] || { echo "[ERROR] Funcotator data sources not found: $FUNCOTATOR"; exit 1; }

# BWA index
if [ ! -f "${REF}.bwt" ]; then
    echo "Indexing reference with BWA..."
    bwa index "$REF"
fi
# samtools .fai
if [ ! -f "${REF}.fai" ]; then
    samtools faidx "$REF"
fi
# GATK sequence dictionary
REF_DICT="${REF%.fasta}.dict"
if [ ! -f "$REF_DICT" ]; then
    gatk CreateSequenceDictionary -R "$REF" -O "$REF_DICT"
fi
# Known-sites indexes (needed by BaseRecalibrator)
for vcf in "$DBSNP" "$INDELS" "$MILLS"; do
    if [ ! -f "${vcf}.tbi" ] && [ ! -f "${vcf}.idx" ]; then
        echo "Indexing $vcf ..."
        gatk IndexFeatureFile -I "$vcf"
    fi
done

#=======================================
# FASTQC (raw reads)
#=======================================
echo "FastQC is running on raw reads.........."
fastqc -t "$THREADS" "$FASTQ_TUMOR1" "$FASTQ_TUMOR2" "$FASTQ_NORMAL1" "$FASTQ_NORMAL2" \
  -o "$OUTPUT/qc_reports"

#=======================================
# TRIMMING
#=======================================
echo "Trimming the reads..."
for SAMPLE in "$TUMOR" "$NORMAL"; do
    trimmomatic PE -threads 4 \
      "${SAMPLE}_1.fastq.gz" "${SAMPLE}_2.fastq.gz" \
      "$OUTPUT/trimmed_reads/${SAMPLE}_R1_paired.fq.gz" "$OUTPUT/trimmed_reads/${SAMPLE}_R1_unpaired.fq.gz" \
      "$OUTPUT/trimmed_reads/${SAMPLE}_R2_paired.fq.gz" "$OUTPUT/trimmed_reads/${SAMPLE}_R2_unpaired.fq.gz" \
      ILLUMINACLIP:"$ADAPTERS":2:30:10 LEADING:3 TRAILING:3 MINLEN:36
done

#=======================================
# FASTQC on trimmed reads
#=======================================
echo "Running FastQC on trimmed reads........."
fastqc -t "$THREADS" \
  "$OUTPUT/trimmed_reads/${TUMOR}_R1_paired.fq.gz" \
  "$OUTPUT/trimmed_reads/${TUMOR}_R2_paired.fq.gz" \
  "$OUTPUT/trimmed_reads/${NORMAL}_R1_paired.fq.gz" \
  "$OUTPUT/trimmed_reads/${NORMAL}_R2_paired.fq.gz" \
  -o "$OUTPUT/qc_reports"

#=======================================
# Aligning with bwa
#=======================================
echo "Aligning reads........"
for SAMPLE in "$TUMOR" "$NORMAL"; do
    bwa mem -t "$THREADS" -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" "$REF" \
      "$OUTPUT/trimmed_reads/${SAMPLE}_R1_paired.fq.gz" \
      "$OUTPUT/trimmed_reads/${SAMPLE}_R2_paired.fq.gz" \
      | samtools sort -@ 4 -T "$OUTPUT/bam_files/${SAMPLE}_tmp" \
          -o "$OUTPUT/bam_files/${SAMPLE}_sorted.bam"
done

#=======================================
# Indexing with samtools
#=======================================
echo "Indexing sorted bams........."
samtools index "$OUTPUT/bam_files/${TUMOR}_sorted.bam"
samtools index "$OUTPUT/bam_files/${NORMAL}_sorted.bam"

#=======================================
# Marking duplicates with gatk
#=======================================
echo "Marking duplicates.............."
for SAMPLE in "$TUMOR" "$NORMAL"; do
    gatk MarkDuplicates \
      -I "$OUTPUT/bam_files/${SAMPLE}_sorted.bam" \
      -O "$OUTPUT/bam_files/${SAMPLE}_dedup.bam" \
      -M "$OUTPUT/logs/${SAMPLE}_metrics.txt"
    samtools index "$OUTPUT/bam_files/${SAMPLE}_dedup.bam"
done

#=======================================
# BQSR
#=======================================
echo "BQSR......................"
for SAMPLE in "$TUMOR" "$NORMAL"; do
    gatk BaseRecalibrator \
      -I "$OUTPUT/bam_files/${SAMPLE}_dedup.bam" -R "$REF" \
      --known-sites "$DBSNP" --known-sites "$INDELS" --known-sites "$MILLS" \
      -L "$PANEL" --interval-padding 100 \
      -O "$OUTPUT/logs/${SAMPLE}_recal.table"

    gatk ApplyBQSR \
      -I "$OUTPUT/bam_files/${SAMPLE}_dedup.bam" -R "$REF" \
      --bqsr-recal-file "$OUTPUT/logs/${SAMPLE}_recal.table" \
      -O "$OUTPUT/bam_files/${SAMPLE}_recal.bam"

    samtools index "$OUTPUT/bam_files/${SAMPLE}_recal.bam"
done

#=======================================
# Paths for VCF outputs
#=======================================
RAW_VCF="$OUTPUT/vcf_files/${TUMOR}_vs_${NORMAL}_somatic_raw.vcf.gz"
F1R2="$OUTPUT/vcf_files/${TUMOR}_f1r2.tar.gz"
ORIENT_MODEL="$OUTPUT/vcf_files/${TUMOR}_read-orientation-model.tar.gz"
PILEUP="$OUTPUT/vcf_files/${TUMOR}_pileups.table"
CONTAM="$OUTPUT/vcf_files/${TUMOR}_contamination.table"
FILTERED_VCF="$OUTPUT/vcf_files/${TUMOR}_vs_${NORMAL}_somatic_filtered.vcf.gz"
ANNOTATED_VCF="$OUTPUT/vcf_files/${TUMOR}_vs_${NORMAL}_somatic_annotated.vcf.gz"

#=============================================
# Calling somatic variants with gatk-Mutect2
#=============================================
echo "Calling somatic variants with Mutect2..."
MUTECT_ARGS=(
  -R "$REF"
  -I "$OUTPUT/bam_files/${TUMOR}_recal.bam"  -tumor  "$TUMOR"
  -I "$OUTPUT/bam_files/${NORMAL}_recal.bam" -normal "$NORMAL"
  -L "$PANEL"
  --f1r2-tar-gz "$F1R2"
  -O "$RAW_VCF"
)
[ -f "$GERMLINE_RES" ] && MUTECT_ARGS+=(--germline-resource "$GERMLINE_RES")
[ -f "$PON" ]          && MUTECT_ARGS+=(--panel-of-normals "$PON")

gatk --java-options "-Xmx8g" Mutect2 "${MUTECT_ARGS[@]}"

#=============================================
# Orientation bias model
#=============================================
gatk LearnReadOrientationModel -I "$F1R2" -O "$ORIENT_MODEL"

#=============================================
# Contamination estimate (needs common-SNP resource)
#=============================================
FILTER_ARGS=(-R "$REF" -V "$RAW_VCF" --ob-priors "$ORIENT_MODEL" -O "$FILTERED_VCF")

if [ -f "$COMMON_SNPS" ]; then
    gatk GetPileupSummaries \
      -I "$OUTPUT/bam_files/${TUMOR}_recal.bam" \
      -V "$COMMON_SNPS" -L "$COMMON_SNPS" \
      -O "$PILEUP"
    gatk CalculateContamination -I "$PILEUP" -O "$CONTAM"
    FILTER_ARGS+=(--contamination-table "$CONTAM")
fi

#========================================================
# Filtering somatic variants with gatk-FilterMutectCalls
#========================================================
echo "Filtering Mutect2 calls..."
gatk FilterMutectCalls "${FILTER_ARGS[@]}"

#========================================================
# Indexing the variants (GATK normally does this already)
#========================================================
if [ ! -f "${FILTERED_VCF}.tbi" ]; then
    tabix -p vcf "$FILTERED_VCF"
fi

#========================================================
# ANNOTATING with funcotator
#========================================================
echo "Annotating with Funcotator......................"
gatk Funcotator \
  --variant "$FILTERED_VCF" \
  --reference "$REF" \
  --ref-version hg38 \
  --data-sources-path "$FUNCOTATOR" \
  --output "$ANNOTATED_VCF" \
  --output-file-format VCF

echo "Somatic Variant Calling Pipeline was completed successfully.....!"
