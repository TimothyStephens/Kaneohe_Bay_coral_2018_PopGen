#!/usr/bin/env bash

BAM="${1}"

# Print all info to log file
exec 1> "${BAM}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

export PATH="$PATH:/home/timothy/programs/nQuire_ret20210707"
F=${BAM%*.bam}

## Step 01
run_cmd "nQuire create -q 20 -c 20 -x -b $BAM -o $F"

## Step 02
run_cmd "nQuire denoise -o $F.denoise $F.bin"

## Step 03
run_cmd "nQuire lrdmodel -t 2 $F.denoise.bin $F.bin > $F.lrdmodel"

## Step 04
run_cmd "nQuire histotest $F.bin 1> $F.histotest"
run_cmd "nQuire histotest $F.denoise.bin 1> $F.denoise.histotest"



