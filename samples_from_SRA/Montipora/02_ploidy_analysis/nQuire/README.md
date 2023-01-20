# `nQuire` ploidy analysis for additional *M. capitata* RNA-seq samples from SRA

## Setup analysis directory

Link `bam` files of aligned RNA-seq reads.

```bash
while read F; 
do
  ln -s ../../01_RNAseq_samples_readmapping/$F.SplitNCigarReads.split.bam; 
done < ../../../Montipora_selected_Runs.txt
```

Setup bash environment.

```bash
conda activate py27
```

## Run `nQuire`

Run `nQuire` `create`, `denies`, `lrdmodel`, and `histotest` on each `bam` file. 

```bash
for BAM in *.bam;
do
  ./run_nQuire.sh $BAM
done
```

Check that the jobs completed correct and no errors occurred.

```bash
grep -i 'error\|warn\|fault' *.bam.log.*
grep 'ExitStatus:' *.bam.log.* | awk '{print $2}' | sort | uniq -c
```

Last command should give **27** total.

## Extract and reformat `nQuire` results

Generate a list of samples and their estimated ploidy, sort this list by ploidy and sample ID so that we can use it to order our data later on (both in a table format and for plotting).

```bash
echo -e "ID\tPloidy" > sample_order_by_ploidy.txt

cat ../../../Montipora_selected_Runs.samples.annotations.txt \
  | awk 'NR>1 {print $1"\t"$8}' \
  | sort -k2,2n -k1,1 \
  >> sample_order_by_ploidy.txt
```

Generate a table of the delta-logliklihosd values for each sample using both the normal and denoised data (biallelic sites). Use the list of ordered samples that we just created to build and order the results table. 

```bash
echo -e "ID\tPloidy\tbiallelic_sites\tbiallelic_sites_denoised\tDiploid\tTriploid\tTetraploid\tDiploid_denoised\tTriploid_denoised\tTetraploid_denoised" > nQuire_results.tsv

while read LINE;
do
  ID=$(echo "${LINE}" | awk '{print $1}')
  PL=$(echo "${LINE}" | awk '{print $2}')
  BEFORE=$(awk '$1=="Before:"{print $2}' ${ID}.SplitNCigarReads.split.bam.log.*)
  AFTER=$(awk '$1=="After:" {print $2}' ${ID}.SplitNCigarReads.split.bam.log.*)
  echo -e "${ID} ${PL} ${BEFORE} ${AFTER} "\
$(awk -F'\t' 'NR>1 && $1!~".denoise.bin" {print $6" "$7" "$8}' \
    ${ID}.SplitNCigarReads.split.lrdmodel)\
" "\
$(awk -F'\t' 'NR>1 && $1~".denoise.bin"  {print $6" "$7" "$8}' \
    ${ID}.SplitNCigarReads.split.lrdmodel) \
 | sed -e 's/ /\t/g' \
 >> nQuire_results.tsv
done < <(awk 'NR>1' sample_order_by_ploidy.txt)
```

Download results from server

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0022_Coral_Genotype_Analysis/03_Analysis/2022-02-05/samples_from_SRA/Montipora/02_ploidy_analysis/nQuire/"
rsync -zarv --prune-empty-dirs --relative \
 --include="*/" \
 --include="*.sh*" --include="*.txt" \
 --include="*.tsv" \
 --exclude="*" \
 ${WD}/./* \
 . --dry-run
```





