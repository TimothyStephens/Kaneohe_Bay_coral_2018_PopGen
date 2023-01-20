#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

export PATH="$PATH:/home/timothy/programs/gatk-4.2.0.0"

REF="../../00_databases/Pocillopora_acuta_KBHIv2.assembly.fasta"
OUT="GVCFall"
VARIANTS=""

## Combine *.g.vcf.gz files and call genotypes
while read ID;
do
	VARIANTS="${VARIANTS} --variant ../${ID}.HaplotypeCaller.ploidy2.g.vcf.gz"; 
done < <(cut -f2 ../samples_Pacuta.txt)
while read ID;
do
        VARIANTS="${VARIANTS} --variant ../additional_samples/${ID}.HaplotypeCaller.ploidy2.g.vcf.gz";
done < <(cut -f2 ../additional_samples/Pocillopora_selected_Runs.samples.txt)


run_cmd "gatk --java-options -Xmx4g CombineGVCFs --reference ${REF} --output ${OUT}.g.vcf.gz ${VARIANTS} > ${OUT}.CombineGVCFs.log 2>&1"
run_cmd "gatk --java-options -Xmx4g GenotypeGVCFs --reference ${REF} --output ${OUT}.vcf.gz -V ${OUT}.g.vcf.gz -stand-call-conf 30 --annotation AS_MappingQualityRankSumTest --annotation AS_ReadPosRankSumTest > ${OUT}.GenotypeGVCFs.log 2>&1"


