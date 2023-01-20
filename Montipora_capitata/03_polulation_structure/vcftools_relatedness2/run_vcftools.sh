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
  --remove-indv SRR5453765 \
  --remove-indv SRR5453748 \
  --remove-indv SRR5453751 \
  --remove-indv SRR5453739 \
  --remove-indv SRR5453742 \
  --remove-indv SRR5453741 \
  --remove-indv SRR5453753 \
  --remove-indv SRR5453745 \
  --remove-indv SRR5453747 \
  --remove-indv SRR5453763 \
  --remove-indv SRR5453744 \
  --remove-indv SRR5453743 \
  --remove-indv SRR5453760 \
  --remove-indv SRR5453746 \
  --remove-indv SRR5453757 \
  --remove-indv SRR5453756 \
  --remove-indv SRR5453755 \
  --remove-indv SRR5453740 \
  --remove-indv SRR5453759 \
  --remove-indv SRR5453749 \
  --remove-indv SRR5453758 \
  --remove-indv SRR5453754 \
  --remove-indv SRR5453750 \
  --remove-indv SRR5453762 \
  --remove-indv SRR5453752 \
  --remove-indv SRR5453761 \
  --remove-indv SRR5453764 \
  --remove-indels \
  --min-meanDP 10 --max-missing 1.0 \
  --recode --recode-INFO-all"

# Compress VCF
run_cmd "gzip ${VCF%.vcf.gz}.filtered.recode.vcf"

# Calculate 
run_cmd "vcftools \
  --gzvcf ${VCF%.vcf.gz}.filtered.recode.vcf.gz --out ${VCF%.vcf.gz}.filtered.recode \
  --relatedness2"



