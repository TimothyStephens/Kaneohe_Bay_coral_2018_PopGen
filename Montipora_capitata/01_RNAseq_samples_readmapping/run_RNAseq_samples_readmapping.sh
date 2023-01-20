#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

SAMPLES="samples_Mcapitata.txt"
REF="Montipora_capitata_KBHIv3.assembly.fasta"
STAR_DB="../00_databases/STAR"
NCPUS=8
NPARALLEL=10

RNAseq_samples_readmapping() {
	export PATH="$PATH:/home/timothy/programs/STAR-2.7.8a/bin/Linux_x86_64_static"
	export PATH="$PATH:/home/timothy/programs/gatk-4.2.0.0"
	export PATH="$PATH:/home/timothy/programs/rgsam/bin"
	export PATH="$PATH:/home/timothy/programs/samtools-1.11/bin"
	
	NCPUS="${1}"; shift
	REF="${1}"; shift
	STAR_DB="${1}"; shift
	SAMPLE="${1}"; shift
	
	OUT=$(echo "$SAMPLE" | awk -F'\t' '{print $2}')
	R1=$(echo "$SAMPLE"  | awk -F'\t' '{print $3}')
	R2=$(echo "$SAMPLE"  | awk -F'\t' '{print $4}')
	
	exec 1> "${OUT}.log.$(date +%s)" 2>&1
	
	# Workflow based on
	# https://evodify.com/genomic-variant-calling-pipeline/
	# https://evodify.com/gatk-in-non-model-organism/
	
	
	
	## 01 - StarAlign
	# Align RNAseq reads to genome using STAR aligner. 
	# Use the 2-pass mode becuase it can make the exon-intron junctions slightly more accurate. 
	# See: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
	log "Step 1 - StarAlign"
	if [ ! -f "${OUT}.STAR.done" ];
	then
		run_cmd "STAR \
			--genomeDir ${STAR_DB} \
			--runThreadN ${NCPUS} \
			--readFilesIn ${R1} ${R2} \
			--readFilesCommand \"gunzip -c\" \
			--sjdbOverhang 149 \
			--outSAMtype BAM SortedByCoordinate \
			--twopassMode Basic \
			--outFileNamePrefix ${OUT}.STAR." && \
		run_cmd "touch ${OUT}.STAR.done"
	else
		log "         - Already Fishied!"
	fi
	echo -e '\n\n\n\n\n'
	
	
	
	## 02 - FastqToSam + collect_RG + MergeBamAlignment
	# Need to add read group info to aligned reads + run the MergeBamAlignment command.
	# The MergeBamAlignment is technically not required if we can add read groups another way but it also 
	# filters the alinged read (e.g. removes secondary alignments) so would be hard to replicate this step
	# exactly using other tools (i.e. easier to to just use it then replace it).
	
	# Convert paired-fastq to BAM file (sorted by read name); no read group info added yet.
	log "Step 2.1 - FastqToSam"
	if [ ! -f "${OUT}.FastqToSam.done" ];
	then
		run_cmd "gatk FastqToSam \
				--FASTQ ${R1} \
				--FASTQ2 ${R2} \
				--OUTPUT ${OUT}.FastqToSam.unmapped.bam \
				--SAMPLE_NAME ${OUT}" && \
		run_cmd "touch ${OUT}.FastqToSam.done"
	else
		log "         - Already Fishied!"
	fi
	
	# Extract read group info for read names then annotate unaligned bam file with this read group info.
	# See: https://github.com/djhshih/rgsam
	log "Step 2.2 - collect_RG"
        if [ ! -f "${OUT}.collect_RG.done" ];
        then
                run_cmd "samtools view ${OUT}.FastqToSam.unmapped.bam \
				| rgsam collect --format sam --qnformat illumina-1.8 -s ${OUT} -o ${OUT}.collect_RG.txt" && \
		run_cmd "samtools view -h ${OUT}.FastqToSam.unmapped.bam \
				| rgsam tag --qnformat illumina-1.8 -r ${OUT}.collect_RG.txt \
				| samtools view -b - \
				> ${OUT}.FastqToSam.unmapped.rg.bam" && \
		run_cmd "touch ${OUT}.collect_RG.done"
        else
                log "         - Already Fishied!"
        fi
	
	# Merge unalinged bam file (now with read group info) with aligned bam file (read group info from unalinged bam is transfered to aligned bam).
	log "Step 2.3 - MergeBamAlignment"
	if [ ! -f "${OUT}.MergeBamAlignment.done" ];
        then
                run_cmd "gatk MergeBamAlignment \
					--REFERENCE_SEQUENCE ${REF} \
					--UNMAPPED_BAM ${OUT}.FastqToSam.unmapped.rg.bam \
					--ALIGNED_BAM ${OUT}.STAR.Aligned.sortedByCoord.out.bam \
					--OUTPUT ${OUT}.MergeBamAlignment.merged.bam \
					--INCLUDE_SECONDARY_ALIGNMENTS false \
					--VALIDATION_STRINGENCY SILENT" && \
		run_cmd "touch ${OUT}.MergeBamAlignment.done"
        else
                log "         - Already Fishied!"
        fi
	echo -e '\n\n\n\n\n'
	
	
	
	## 03 - MarkDuplicates
	log "Step 3 - MarkDuplicates"
	if [ ! -f "${OUT}.MarkDuplicates.done" ];
	then
		run_cmd "gatk MarkDuplicates \
					--INPUT ${OUT}.MergeBamAlignment.merged.bam \
					--OUTPUT ${OUT}.MarkDuplicates.dedupped.bam \
					--CREATE_INDEX true \
					--VALIDATION_STRINGENCY SILENT \
					--METRICS_FILE ${OUT}.MarkDuplicates.metrics" && \
		run_cmd "touch ${OUT}.MarkDuplicates.done"
	else
		log "         - Already Fishied!"
	fi
	echo -e '\n\n\n\n\n'
	
	
	
	## 04 - SplitNCigarReads
	log "Step 4 - SplitNCigarReads"
	if [ ! -f "${OUT}.SplitNCigarReads.done" ];
	then
		run_cmd "gatk SplitNCigarReads \
					-R ${REF} \
					-I ${OUT}.MarkDuplicates.dedupped.bam \
					-O ${OUT}.SplitNCigarReads.split.bam" && \
		run_cmd "touch ${OUT}.SplitNCigarReads.done"
	else
		log "         - Already Fishied!"
	fi
	echo -e '\n\n\n\n\n'
	
}


export -f RNAseq_samples_readmapping
parallel -j $NPARALLEL RNAseq_samples_readmapping "$NCPUS" "$REF" "$STAR_DB" :::: "$SAMPLES"
log "Parallel finished with exit status: $?"

echo -e "## Done running!"

