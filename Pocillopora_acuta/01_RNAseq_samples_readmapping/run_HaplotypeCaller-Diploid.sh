#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

SAMPLES=""samples_Pacuta.txt
REF="Pocillopora_acuta_KBHIv2.assembly.fasta"
NCPUS=4
NPARALLEL=20

RNAseq_short_variant_analysis() {
	export PATH="$PATH:/home/timothy/programs/gatk-4.2.0.0"
	
	PLOIDY=2
	NCPUS="${1}"; shift
	REF="${1}"; shift
	SAMPLE="${1}"; shift
	
	OUT=$(echo "$SAMPLE" | awk -F'\t' '{print $2}')
	R1=$(echo "$SAMPLE"  | awk -F'\t' '{print $3}')
	R2=$(echo "$SAMPLE"  | awk -F'\t' '{print $4}')
	
	exec 1> "${OUT}.HaplotypeCaller.ploidy${PLOIDY}.log.$(date +%s)" 2>&1
	
	# Workflow based on
	# https://evodify.com/genomic-variant-calling-pipeline/
	# https://evodify.com/gatk-in-non-model-organism/
	
	
	## 08 - HaplotypeCaller
	# Ignore flags:
	#    -L ${interval_list}	Want to call variants from the whole genome
	#    --dbsnp ${dbSNP_vcf}	Dont have a preexisting database of SNPs avaliable
	#    --standard-min-confidence-threshold-for-calling 20		When running in '-ERC GVCF' mode, the confidence threshold is automatically set to 0.
	#
	# Assumes:
	#	--heterozygosity 0.001 (deafult; dont have prior info to update this with)
	log "Step 5 - HaplotypeCaller"
	if [ ! -f "${OUT}.HaplotypeCaller.done" ];
	then
		run_cmd "gatk HaplotypeCaller \
			--reference ${REF} \
			--input ${OUT}.SplitNCigarReads.split.bam \
			--output ${OUT}.HaplotypeCaller.ploidy${PLOIDY}.g.vcf.gz \
			-dont-use-soft-clipped-bases \
			-ERC GVCF \
			--sample-ploidy ${PLOIDY} \
			--native-pair-hmm-threads ${NCPUS}" && \
		run_cmd "touch ${OUT}.HaplotypeCaller.done"
	else
		log "         - Already Fishied!"
	fi
	echo -e '\n\n\n\n\n'
}


export -f RNAseq_short_variant_analysis
parallel -j $NPARALLEL RNAseq_short_variant_analysis "$NCPUS" "$REF" :::: "$SAMPLES"
log "Parallel finished with exit status: $?"

echo -e "## Done running!"

