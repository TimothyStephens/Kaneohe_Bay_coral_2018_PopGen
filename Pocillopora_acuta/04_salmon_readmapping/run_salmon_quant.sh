#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

export PATH="$PATH:/home/timothy/programs/salmon-1.6.0_linux_x86_64/bin"
NCPUS=48

## Index info
TRANSCRIPTS="Pocillopora_acuta_KBHIv2.gentrome.fa.gz"
INDEX="$TRANSCRIPTS.salmon.idx"
## Read orientation
# I = inward
# S = stranded
# R = read 1 (or single-end read) comes from the reverse strand
LIB_TYPE="ISR"

#### Start Script
while read LINE;
do
	R1=$(echo $LINE | awk '{print $3}')
	R2=$(echo $LINE | awk '{print $4}')
	OUTPUT=$(echo $LINE | awk '{print $2}')
	run_cmd "salmon quant --validateMappings --seqBias --gcBias --index $INDEX --libType $LIB_TYPE --mates1 $R1 --mates2 $R2 --output $OUTPUT --threads $NCPUS"
done < samples_Pacuta.txt


