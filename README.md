# MELANOMA-Targeted-Variant-Analysis
Identification and Characterization of Germline and Somatic Variants in Melanoma Through Targeted Sequencing
## Overview

This repository contains an end-to-end bioinformatics workflow for targeted panel sequencing analysis of melanoma samples.

The pipeline performs:

- Quality control of raw FASTQ files
- Adapter and quality trimming
- Read alignment to the human reference genome (GRCh38)
- BAM processing and quality assessment
- Germline variant calling
- Somatic variant calling using matched tumour-normal samples
- Variant filtering
- Functional variant annotation
- Variant interpretation

Workflow
FASTQ
   │
   ▼
Quality Control (FastQC)
   │
   ▼
Read Trimming (fastp / Trimmomatic)
   │
   ▼
Alignment (BWA-MEM)
   │
   ▼
Sorted BAM
   │
   ▼
Duplicate Marking
   │
   ▼
Base Quality Score Recalibration (BQSR)
          │
   ├───────────────┐
   ▼               ▼
Germline         Somatic
HaplotypeCaller  Mutect2
   │               │
   ▼               ▼
     Filtered VCFs
          │
          ▼
  Annotation (Funcotator)
          │
          ▼
 Final Annotated Variants


Targeted Panels
| Panel                       | Description                                      |
| --------------------------- | ------------------------------------------------ |
| Ion AmpliSeq Melanoma Panel | 22 Targeted sequencing of melanoma-associated genes |

Software
| Tool            | Purpose                               |
| --------------- | ------------------------------------- |
| FastQC          | Raw read quality assessment           |
| fastp           | Adapter and quality trimming          |
| BWA-MEM         | Read alignment                        |
| SAMtools        | BAM processing and quality assessment |
| Picard          | Duplicate marking and BAM processing  |
| GATK            | Variant calling and filtering         |
| HaplotypeCaller | Germline variant detection            |
| Mutect2         | Somatic variant detection             |
| Funcotator      | Functional variant annotation         |
| bcftools        | VCF processing and variant filtering  |

Reproducibility
The complete workflows are provided as shell scripts in the workflow_script/ directory.

Citation
If you use this workflow in your research or academic project, please cite this repository.
