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
  --remove-indv SRR7059699 \
  --remove-indv SRR7062295 \
  --remove-indv SRR7062350 \
  --remove-indv SRR7039808 \
  --remove-indv SRR7040514 \
  --remove-indv SRR7041301 \
  --remove-indv SRR7042978 \
  --remove-indv SRR7043704 \
  --remove-indv SRR7055829 \
  --remove-indv SRR7058378 \
  --remove-indv SRR7058566 \
  --remove-indv SRR7058606 \
  --remove-indv SRR7058616 \
  --remove-indv SRR6951423 \
  --remove-indv SRR6951744 \
  --remove-indv SRR6952431 \
  --remove-indv SRR6963586 \
  --remove-indv SRR6963878 \
  --remove-indv SRR6963891 \
  --remove-indv SRR6964364 \
  --remove-indv SRR6986864 \
  --remove-indv SRR6987146 \
  --remove-indv SRR6914609 \
  --remove-indv SRR6914151 \
  --remove-indv SRR6914908 \
  --remove-indv SRR6934388 \
  --remove-indv SRR6934542 \
  --remove-indv SRR6935629 \
  --remove-indv SRR6942678 \
  --remove-indv SRR6942729 \
  --remove-indv SRR7043013 \
  --remove-indv SRR7046161 \
  --remove-indels \
  --min-meanDP 10 --max-missing 1.0 \
  --recode --recode-INFO-all"

# Compress VCF
run_cmd "gzip ${VCF%.vcf.gz}.filtered.recode.vcf"

# Calculate 
run_cmd "vcftools \
  --gzvcf ${VCF%.vcf.gz}.filtered.recode.vcf.gz --out ${VCF%.vcf.gz}.filtered.recode \
  --relatedness2"



