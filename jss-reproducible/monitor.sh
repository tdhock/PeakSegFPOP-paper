#!/bin/bash
ls jss.bench.models.rules|wc -l
wc -l jss.bench.models.csv
squeue -hru thocking|sed 's/ \+/\t/g'|cut -f5|sort|uniq -c
