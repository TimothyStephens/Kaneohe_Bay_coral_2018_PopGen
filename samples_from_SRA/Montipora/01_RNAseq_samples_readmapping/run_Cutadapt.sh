#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate cutadapt29; set -eu

SRR=$(cut -f1 Montipora_selected_Runs.samples.txt)
INDIR="prefetch_sra_files"
OUTDIR="prefetch_sra_files"
NCPUS=80

#### Start Script
# Similar to other samples except swapped --nextseq-trim for -q and used slightly shorter adapter sequences (See https://cutadapt.readthedocs.io/en/stable/guide.html#illumina-truseq)
mkdir -p "$OUTDIR"
for ID in $SRR;
do
	run_cmd "cutadapt -q 10 --minimum-length 25 --cores $NCPUS -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o ${OUTDIR}/${ID}_R1_Trimmed.fastq.gz -p ${OUTDIR}/${ID}_R2_Trimmed.fastq.gz ${INDIR}/${ID}.sra_1.fastq.gz ${INDIR}/${ID}.sra_2.fastq.gz > ${INDIR}/${ID}.cutadapt.log 2>&1"
done



