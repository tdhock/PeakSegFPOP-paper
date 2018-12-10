#!/bin/bash
ls jss.bench.models.rules|wc -l
wc -l jss.bench.models.csv
R -q -e "slurm::sjob()"
##echo "sjob($job)" > monitor.R
##Rscript monitor.R
