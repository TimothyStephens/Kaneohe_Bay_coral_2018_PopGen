#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate angsd; set -eu

NCPUS=48


#### Start Script
run_cmd "angsd -gl 2 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.05 -doGlf 3 -nThreads ${NCPUS} -bam bam.filelist -out NgsRelate.angsd"



