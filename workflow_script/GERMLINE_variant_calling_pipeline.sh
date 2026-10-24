#!/bin/bash

# Usage: ./germline_vc_pipeline_fixed.sh <sample_id>
# Example: ./germline_vc_pipeline_fixed.sh SRR1656612

set -euo pipefail

#=======================================
# Inputs (taken from command line)
#=======================================
SAMPLE="${1:?Usage: $0 <sample_id>}"

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

THREADS=8
JAVA_MEM="-Xmx14g"

cd "$WORKDIR"

#=======================================
# FASTQ naming convention
#=======================================
FASTQ1="${SAMPLE}_1.fastq.gz"
FASTQ2="${SAMPLE}_2.fastq.gz"

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
for f in "$FASTQ1" "$FASTQ2" "$REF" "$DBSNP" "$INDELS" "$MILLS" "$PANEL" "$ADAPTERS"; do
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
fastqc -t "$THREADS" "$FASTQ1" "$FASTQ2" -o "$OUTPUT/qc_reports"

#=======================================
# TRIMMING
#=======================================
echo "Trimming the reads..."
trimmomatic PE -threads 4 \
  "$FASTQ1" "$FASTQ2" \
  "$OUTPUT/trimmed_reads/${SAMPLE}_R1_paired.fq.gz" "$OUTPUT/trimmed_reads/${SAMPLE}_R1_unpaired.fq.gz" \
  "$OUTPUT/trimmed_reads/${SAMPLE}_R2_paired.fq.gz" "$OUTPUT/trimmed_reads/${SAMPLE}_R2_unpaired.fq.gz" \
  ILLUMINACLIP:"$ADAPTERS":2:30:10 LEADING:3 TRAILING:3 MINLEN:36

#=======================================
# FASTQC on trimmed reads
#=======================================
echo "Running FastQC on trimmed reads........."
fastqc -t "$THREADS" \
  "$OUTPUT/trimmed_reads/${SAMPLE}_R1_paired.fq.gz" \
  "$OUTPUT/trimmed_reads/${SAMPLE}_R2_paired.fq.gz" \
  -o "$OUTPUT/qc_reports"

#=======================================
# Aligning with bwa
#=======================================
echo "Aligning reads........"
bwa mem -t "$THREADS" -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" "$REF" \
  "$OUTPUT/trimmed_reads/${SAMPLE}_R1_paired.fq.gz" \
  "$OUTPUT/trimmed_reads/${SAMPLE}_R2_paired.fq.gz" \
  | samtools sort -@ 4 -T "$OUTPUT/bam_files/${SAMPLE}_tmp" \
      -o "$OUTPUT/bam_files/${SAMPLE}_sorted.bam"

samtools index "$OUTPUT/bam_files/${SAMPLE}_sorted.bam"

#=======================================
# Marking duplicates with gatk
#=======================================
echo "Marking duplicates.............."
gatk MarkDuplicates \
  -I "$OUTPUT/bam_files/${SAMPLE}_sorted.bam" \
  -O "$OUTPUT/bam_files/${SAMPLE}_dedup.bam" \
  -M "$OUTPUT/logs/${SAMPLE}_metrics.txt"

samtools index "$OUTPUT/bam_files/${SAMPLE}_dedup.bam"

#=======================================
# BQSR
#=======================================
echo "BQSR......................"
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

#=======================================
# Paths for VCF outputs
#=======================================
RAW_GVCF="$OUTPUT/vcf_files/${SAMPLE}_raw.g.vcf.gz"
GENOTYPED_VCF="$OUTPUT/vcf_files/${SAMPLE}_genotyped.vcf.gz"
RAW_SNP_VCF="$OUTPUT/vcf_files/${SAMPLE}_raw_snps.vcf.gz"
RAW_INDEL_VCF="$OUTPUT/vcf_files/${SAMPLE}_raw_indels.vcf.gz"
FILTERED_SNP_VCF="$OUTPUT/vcf_files/${SAMPLE}_filtered_snps.vcf.gz"
FILTERED_INDEL_VCF="$OUTPUT/vcf_files/${SAMPLE}_filtered_indels.vcf.gz"
MERGED_FILTERED_VCF="$OUTPUT/vcf_files/${SAMPLE}_filtered_merged.vcf.gz"
ANNOTATED_VCF="$OUTPUT/vcf_files/${SAMPLE}_annotated.vcf.gz"

#====================================================
# Calling germline variants with gatk-HaplotypeCaller
#====================================================
echo "Germline variant calling with HaplotypeCaller.........."
gatk --java-options "$JAVA_MEM" HaplotypeCaller \
  -R "$REF" \
  -I "$OUTPUT/bam_files/${SAMPLE}_recal.bam" \
  -L "$PANEL" --interval-padding 100 \
  -O "$RAW_GVCF" \
  -ERC GVCF

#====================================================
# GENOTYPING GVCF
#====================================================
echo "Genotyping GVCF..."
gatk --java-options "$JAVA_MEM" GenotypeGVCFs \
  -R "$REF" \
  -V "$RAW_GVCF" \
  -L "$PANEL" --interval-padding 100 \
  -D "$DBSNP" \
  -O "$GENOTYPED_VCF"

#====================================================
# FILTERING SNPs (GATK hard-filter recommendations)
#====================================================
echo "FILTERING SNPs........."
gatk SelectVariants \
  -R "$REF" \
  -V "$GENOTYPED_VCF" \
  --select-type-to-include SNP \
  -O "$RAW_SNP_VCF"

gatk VariantFiltration \
  -R "$REF" \
  -V "$RAW_SNP_VCF" \
  --filter-expression "QD < 2.0"              --filter-name "QD2" \
  --filter-expression "FS > 60.0"             --filter-name "FS60" \
  --filter-expression "MQ < 40.0"             --filter-name "MQ40" \
  --filter-expression "SOR > 3.0"             --filter-name "SOR3" \
  --filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  -O "$FILTERED_SNP_VCF"

#====================================================
# FILTERING INDELS
#====================================================
echo "FILTERING INDELS........."
gatk SelectVariants \
  -R "$REF" \
  -V "$GENOTYPED_VCF" \
  --select-type-to-include INDEL \
  -O "$RAW_INDEL_VCF"

gatk VariantFiltration \
  -R "$REF" \
  -V "$RAW_INDEL_VCF" \
  --filter-expression "QD < 2.0"                --filter-name "QD2" \
  --filter-expression "FS > 200.0"              --filter-name "FS200" \
  --filter-expression "SOR > 10.0"              --filter-name "SOR10" \
  --filter-expression "ReadPosRankSum < -20.0"  --filter-name "ReadPosRankSum-20" \
  -O "$FILTERED_INDEL_VCF"

#====================================================
# MERGING FILTERED SNPS AND INDELS
#====================================================
echo "Merging filtered SNPs and INDELs..."
gatk MergeVcfs \
  -I "$FILTERED_SNP_VCF" \
  -I "$FILTERED_INDEL_VCF" \
  -O "$MERGED_FILTERED_VCF"

# GATK normally indexes .vcf.gz outputs already; this is a safety net
if [ ! -f "${MERGED_FILTERED_VCF}.tbi" ]; then
    tabix -p vcf "$MERGED_FILTERED_VCF"
fi

#====================================================
# ANNOTATING VARIANTS WITH FUNCOTATOR
#====================================================
echo "Annotating filtered variants with Funcotator....."
gatk Funcotator \
  --variant "$MERGED_FILTERED_VCF" \
  --reference "$REF" \
  --ref-version hg38 \
  --data-sources-path "$FUNCOTATOR" \
  --output "$ANNOTATED_VCF" \
  --output-file-format VCF

echo "Germline Variant Calling Pipeline was completed successfully.....!" 
