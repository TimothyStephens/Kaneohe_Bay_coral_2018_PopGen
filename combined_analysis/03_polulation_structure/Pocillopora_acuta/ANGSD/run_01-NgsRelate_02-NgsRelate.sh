#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate angsd; set -eu
export PATH="$PATH:/home/timothy/programs/ngsRelate"

NCPUS=48
MAFS="NgsRelate.angsd.mafs.gz"
GLF="NgsRelate.angsd.glf.gz"

#### Start Script
# See https://github.com/ANGSD/NgsRelate
run_cmd "zcat ${MAFS} | cut -f5 |sed 1d > ${MAFS}.freq"
run_cmd "ngsRelate -g ${GLF} -n 151 -f ${MAFS}.freq -z bam.filelist.labels -O NgsRelate.results.tsv -p ${NCPUS}"

