#!/bin/bash
#SBATCH --job-name=BenchModels
#SBATCH --workdir=/scratch/thocking/PeakSegFPOP-paper/jss-reproducible
#SBATCH --output=/scratch/thocking/PeakSegFPOP-paper/jss-reproducible/logs/job_%A_%a.log
#SBATCH --time=72:00:00
#SBATCH --array=4
#SBATCH --mem=2G
#SBATCH --account=rrg-bourqueg-ad
#SBATCH --cpus-per-task=1
set -o errexit
module load r/3.5.0
Rscript jss.bench.models.several.one.R $SLURM_ARRAY_TASK_ID
