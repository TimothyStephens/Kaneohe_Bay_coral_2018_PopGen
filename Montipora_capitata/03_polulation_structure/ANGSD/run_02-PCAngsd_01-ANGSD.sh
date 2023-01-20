#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate angsd; set -eu

NCPUS=48


#### Start Script
run_cmd "angsd -GL 2 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -nThreads ${NCPUS} -bam bam.filelist -out PCAngsd.angsd"



