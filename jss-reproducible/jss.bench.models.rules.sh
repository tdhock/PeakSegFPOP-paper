#!/bin/bash
#SBATCH --job-name=BenchRules
#SBATCH --workdir=/scratch/thocking/PeakSegFPOP-paper/jss-reproducible
#SBATCH --output=/scratch/thocking/PeakSegFPOP-paper/jss-reproducible/logs/job_%A_%a.log
#SBATCH --time=24:00:00
#SBATCH --array=1
#SBATCH --mem=2G
#SBATCH --account=rrg-bourqueg-ad
#SBATCH --cpus-per-task=1
set -o errexit
module load r/3.5.0
Rscript jss.bench.models.rules.R
