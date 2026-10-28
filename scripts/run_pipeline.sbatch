#!/bin/bash
#SBATCH --job-name=nf_wgs
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/nextflow_wgs_%j.out
#SBATCH --error=logs/nextflow_wgs_%j.err

module load nextflow/24.10.3
module load anaconda3/2024.06

# Submit from the project folder:  sbatch scripts/run_pipeline.sbatch
cd $SLURM_SUBMIT_DIR

# runs the main.nf pipeline with Conda environments
nextflow run main.nf -profile conda -resume
