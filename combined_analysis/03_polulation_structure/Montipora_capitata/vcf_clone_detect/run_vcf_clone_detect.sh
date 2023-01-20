#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py38; set -eu

# git clone https://github.com/pimbongaerts/radseq.git
RAD="/home/timothy/programs/pimbongaerts_radseq"

VCF="GVCFall.filtered.recode.vcf.gz"


#### Start Script
run_cmd "gunzip -c ${VCF} > ${VCF%*.gz}"
run_cmd "python3 ${RAD}/vcf_clone_detect.py --vcf ${VCF%.gz} --output ${VCF}.allelic_similarity.csv --threshold 95.0"



