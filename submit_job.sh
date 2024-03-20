#!/bin/bash

nohup snakemake --rerun-incomplete -j 1 \
         --latency-wait 30 \
         --use-conda \
         --executor cluster-generic \
         --cluster-generic-submit-cmd 'qsub -pe smp 64 -l h_vmem=30G -V -cwd'