#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

export PATH="$PATH:/home/timothy/programs/vcftools/bin"
VCF="GVCFall.vcf.gz"


#### Start Script
# Filter variants
run_cmd "vcftools \
  --gzvcf ${VCF} --out ${VCF%.vcf.gz}.filtered \
  --remove-indels \
  --min-meanDP 10 --max-missing 1.0 \
  --recode --recode-INFO-all"

# Compress VCF
run_cmd "gzip ${VCF%.vcf.gz}.filtered.recode.vcf"

# Calculate 
run_cmd "vcftools \
  --gzvcf ${VCF%.vcf.gz}.filtered.recode.vcf.gz --out ${VCF%.vcf.gz}.filtered.recode \
  --relatedness2"



