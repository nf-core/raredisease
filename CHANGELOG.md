# nf-core/raredisease: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0 - [xxxx-xx-xx]

### `Added`

- Add `public_aws_ecr` profile for using AWS ECR public gallery images
- GATK's ShiftFasta to generate all the files required for mitochondrial analysis
- Feature to calculate CADD scores for indels
- HmtNote to annotate mitochondria
- MT del script to detect mitochondrial deletions
- `public_aws_ecr` profile for using AWS ECR public gallery images

## v1.0.0 - [2023-06-01]

Initial release of nf-core/raredisease, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- FastQC read quality control
- Read mapping with BWAmem2/Sentieon
- Qualimap & Picard tools quality control metrics
- Call repeat expansions with ExpansionHunter and Stranger
- SNV calling with DeepVariant/Sentieon
- SV calling with Manta and TIDDIT
- SNV annotation with bcftools roh, vcfanno, and vep
- SV annotation with SVDB query and vep
- Separate workflow for analysing and annotating mitochondrial variants
- Call copy number variants in the SMN gene using SMNCopyNumberCaller
