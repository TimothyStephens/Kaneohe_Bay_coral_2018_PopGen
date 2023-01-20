#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

export PATH="$PATH:/home/timothy/programs/nQuire_ret20210707"
BAM="SRR16077714.coordsorted.bam"
BED="Pocillopora_acuta_KBHIv2.assembly.fasta.1Mbp_scaffold_features.bed"
F=${BAM%*.bam}

## Step 01
run_cmd "nQuire create -q 20 -c 20 -x -b $BAM -o $F -r $BED"

## Step 02
for BIN in ${F}-*.bin;
do
  B=${BIN%.bin}
  run_cmd "nQuire denoise -o $B.denoise $B.bin > $B.denoise.counts.txt"
done

## Step 03
run_cmd "nQuire lrdmodel -t 2 *.bin > $BED.lrdmodel"



