#!/usr/bin/env bash

BAM="${1}"

# Print all info to log file
exec 1> "${BAM%*.bam}.nQuire.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

export PATH="$PATH:/home/timothy/programs/nQuire_ret20210707"
F=${BAM%*.bam}

## Step 01 - create *.bin file
run_cmd "nQuire create -f 0.05 -q 20 -c 20 -x -b $BAM -o $F"

## Step 02 - denoise biallelic sites
run_cmd "nQuire denoise -o $F.denoise $F.bin"

## Step 03 - Extract biallelic sites from denoised data.
run_cmd "nQuire view $F.denoise.bin -a $BAM | gzip -c > $F.denoise.bin.coverage.gz"
# Output format:
# 1. Reference sequence (ID)
# 2. Reference position (0-based)
# 3. Coverage
# 4. Base count 1
# 5. Base count 2

## Step 04 - Print proportion of reads supporting each allele at each bi-allelic site (so each site will have two proportions printed)
run_cmd "zcat $F.denoise.bin.coverage.gz | awk -F'\\t' '{print \$4/\$3\"\\n\"\$5/\$3}' | gzip -c > $F.denoise.bin.coverage.sitesProp.gz"



