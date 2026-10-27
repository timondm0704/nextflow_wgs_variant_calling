# Whole Genome Sequencing Variant Calling Pipeline with Nextflow

A Nextflow pipeline for paired-end whole-genome sequencing (WGS) data, tested on Northeastern University's HPC cluster Explorer.

## Background

This project started as an exercise in BINF 6310 (Northeastern University), where I wrote single-step Nextflow scripts for read trimming with cutadapt and alignment with Bowtie2. I extended it into a multi-step WGS pipeline with a sample sheet, paired-end reads, BWA-MEM alignment, a coverage check, joint variant calling, and multi-sample vcf file.

## Objective

Build a pipeline that takes raw sequencing reads from multiple individuals to a multi-sample VCF of SNPs, with every sample processed the same way:

- **Sample sheet input:** each sample is one row with its FASTQ files and metadata (site, latitude, longitude, depth, date). Adding samples means adding rows, not changing code.
- **Same processing for every sample:** identical trimming, alignment and quality checks.
- **Coverage check before genotyping:** samples below 12x mean coverage are excluded, because low coverage can miss one of the two alleles at heterozygous sites (*heterozygote dropout*).
- **Joint variant calling:** SNPs are called across all passing samples together, producing one VCF.
- **Reproducible and portable:** tools are installed automatically by conda, and the pipeline runs as a SLURM job on an HPC cluster.

## Pipeline

![Pipeline overview](https://timondm0704.github.io/images/nextflow_wgs_pipeline.svg)


### Test data

All sequencing data are simulated with `scripts/simulate_data.py`.

- **Reference:** human *SHH* region from hg38, chr7:155,799,980–155,812,463 (12,484 bp), used as a compact test sequence
- **Samples:** 8 samples, 2 from each of 4 sites (S1–S4)
- **Reads:** paired-end 2 × 100 bp, 0.2% sequencing errors, some low-quality read ends
- **Coverage:** 15x for 7 samples; **6x for `S4_ind2`**, to test the coverage check
- **Truth set:** `data/truth/truth_variants.vcf`, with 39 true SNPs and every sample's genotype

### 1. TRIM

- **Tool:** cutadapt
- **What it does:** trims low-quality bases (below Q20) from the 3′ end of each read and drops read pairs where either read becomes shorter than 50 bp. R1 and R2 are trimmed together so pairs stay matched.
- **Input:** `data/reads/<sample>_R1.fastq.gz`, `data/reads/<sample>_R2.fastq.gz`
- **Output:** `<sample>_R1.trimmed.fastq.gz`, `<sample>_R2.trimmed.fastq.gz` (intermediate, kept in `work/`)

### 2. INDEX

- **Tool:** bwa index, samtools faidx
- **What it does:** builds the search index BWA needs for alignment and the position index bcftools needs for variant calling. Runs once, shared by all samples.
- **Input:** `data/reference/SHH_locus_hg38.fa`
- **Output:** `ref/` folder with `genome.fa` and its index files (intermediate, kept in `work/`)

### 3. ALIGN

- **Tool:** bwa mem, samtools sort, samtools index
- **What it does:** aligns each read pair to the reference, labels reads with the sample name (read group), sorts them by position, and indexes the result.
- **Input:** trimmed FASTQ files from TRIM + reference index from INDEX
- **Output:** `results/bam/<sample>.bam`, `results/bam/<sample>.bam.bai`

### 4. COVERAGE

- **Tool:** samtools coverage
- **What it does:** calculates the mean coverage (average number of reads covering each base) and marks the sample PASS if it is at least 12x, or FAIL if not. Only PASS samples continue.
- **Input:** `results/bam/<sample>.bam` + `.bai`
- **Output:** `results/coverage/<sample>.coverage.txt` (e.g. `S1_ind1 mean_coverage=14.4544 PASS`)

### 5. CALL_VARIANTS

- **Tool:** bcftools mpileup, bcftools call
- **What it does:** looks at all passing samples together and calls SNPs, recording each sample's genotype (0/0, 0/1 or 1/1).
- **Input:** BAM + BAI files of all PASS samples + reference from INDEX
- **Output:** `results/variants/cohort.vcf.gz`


## Folder structure

```
nextflow_wgs_variant_calling/
├── main.nf                 # pipeline
├── nextflow.config         # settings (standard and conda profiles)
├── data/
│   ├── samplesheet.csv     # one row per sample: FASTQ paths + metadata
│   ├── reference/          # reference genome
│   ├── reads/              # 16 FASTQ files (8 samples × R1/R2)
│   └── truth/              # true genotypes
├── scripts/
│   ├── run_pipeline.sbatch # SLURM batch script
│   └── simulate_data.py    # regenerates the test data
├── envs/
│   └── environment.yml     # cutadapt, bwa, samtools, bcftools
├── results/                # pipeline outputs
└── logs/                   # SLURM logs
```

## How to run

Run on a SLURM cluster (e.g. Explorer):

```bash
git clone https://github.com/timondm0704/nextflow_wgs_variant_calling.git
cd nextflow_wgs_variant_calling
sbatch scripts/run_pipeline.sbatch
```

The batch script requests 2 CPUs, 4 GB and 1 hour, loads Nextflow and Anaconda, and runs `nextflow run main.nf -profile conda -resume`. Module names and the partition are specific to Explorer.

**Options:**

- `--samplesheet`: sample sheet (default `data/samplesheet.csv`)
- `--reference`: reference genome (default `data/reference/SHH_locus_hg38.fa`)
- `--min_coverage`: minimum mean coverage to PASS (default `12`)
- `--outdir`: output folder (default `results`)

## Results

| Sample | Mean coverage | Status |
|---|---|---|
| S1_ind1 | 14.45x | PASS |
| S1_ind2 | 14.37x | PASS |
| S2_ind1 | 14.42x | PASS |
| S2_ind2 | 14.38x | PASS |
| S3_ind1 | 14.41x | PASS |
| S3_ind2 | 14.41x | PASS |
| S4_ind1 | 14.44x | PASS |
| S4_ind2 | 5.79x | **FAIL** (excluded) |

- The low-coverage sample was correctly excluded from variant calling. (Samples were simulated at 15x; trimming lowers this to about 14.4x.)
- **40 SNPs** were called across the 7 passing samples, compared with **39 true SNPs** in the simulated data.

To view the first VCF records:

```bash
zcat results/variants/cohort.vcf.gz | grep -v '^##' | head -5
```

Example VCF record:

```
#CHROM     POS  REF  ALT  QUAL     ...  S1_ind2            S4_ind1
SHH_locus  506  C    G    643.762  ...  1/1:255,45,0:0,15  0/1:195,0,119:5,10
```

At position 506, `S1_ind2` is homozygous for G (`1/1`, 15 reads with G) and `S4_ind1` is heterozygous (`0/1`, 5 reads with C and 10 with G).

How to read a VCF record:

- `POS`: position on the reference
- `REF` / `ALT`: the reference letter and the alternate letter
- `QUAL`: confidence that the variant is real (higher is better)
- Each sample column, e.g. `0/1:195,0,119:5,10`, has three parts:
  - `GT` (genotype): `0/0` = two reference copies, `0/1` = one of each, `1/1` = two alternate copies
  - `PL`: likelihoods of 0/0, 0/1 and 1/1 (0 = most likely)
  - `AD`: number of reads supporting the reference and alternate letter

## Design choices

- **Nextflow:** parallel execution, `-resume` after failures, isolated work directories, and the same code locally or on SLURM.
- **BWA-MEM:** a standard aligner for WGS variant calling.
- **12x coverage threshold:** below this, heterozygous genotypes become unreliable.
- **Joint calling with bcftools:** one VCF with consistent sites across all samples.
- **conda environment file:** records and installs the required software for reproducible runs.

## Next steps

- Variant quality filtering (e.g. `bcftools view -i 'QUAL>=30'`)
- Duplicate marking (`samtools markdup` or Picard)
- Read QC reports (FastQC/fastp + MultiQC)
- An automatic accuracy check against the truth set
- A SLURM profile so each step runs as its own cluster job
