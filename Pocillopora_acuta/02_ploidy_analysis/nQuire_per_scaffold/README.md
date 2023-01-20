# `nQuire` ploidy analysis of *P. acuta* reference genome scaffolds

Run `nQuire` on each of the largest *P. acuta* scaffolds to see if there is a change in ploidy across the genome (similar to aneuploidy).

## Setup analysis directory

Link to `bam` file of aligned DNA-seq reads.

```bash
ln -s ../../../../../../0031_Coral_genomes/03_Analysis/2022-02-03/coral_genomes/Pocillopora_acuta_KBHIv2/11_omics_data/genomic/assembly_read_mapping/DNA-seq_PRJNA761443_Illumina_bowtie2/SRR16077714.coordsorted.bam
```

Link to script for text parsing.

```bash
ln -s ../../../../../02_Scripts/add_value_to_table.py
```

Setup bash environment.

```bash
conda activate py27
```

## Run `nQuire`

Get `bed` features that cover each scaffolds with a length > 1 Mbp.

```bash
awk -F'\t' '$2>1000000 {print $1"\t0\t"$2"\t"$1}' \
    Pocillopora_acuta_KBHIv2.assembly.fasta.fai \
  | sort -k3,3nr | sed -e 's/\tPocillopora_acuta_KBHIv2___/\t/' \
  > Pocillopora_acuta_KBHIv2.assembly.fasta.1Mbp_scaffold_features.bed
```

Run `nQuire` `create`, `denoies`, and  `lrdmodel`on the `bam` file. Will create a separate `.bin` file for each region in the supplied `bed` file.

```bash
./run_nQuire.sh
```

## Extract and reformat `nQuire` results

```bash
F="Pocillopora_acuta_KBHIv2.assembly.fasta.1Mbp_scaffold_features.bed.lrdmodel"
```

Print the estimated ploidy of each scaffold.

```bash
cat "$F" | awk 'BEGIN{OFS=FS="\t"} { if(NR==1){$1="sample"; print $0"\tploidy"} else if($6<$7 && $6<$8){print $0"\tDiploid"} else if($7<$6 && $7<$8){print $0"\tTriploid"} else {print $0"\tTetraploid"} }' > "$F.ploidy"
```

Extract the number of raw and denoised balletic sites that we are working with for each scaffold.

```bash
for FILE in *.coordsorted-*.denoise.counts.txt;
do
  ID=$(echo "${FILE}" | sed -e 's/.denoise.counts.txt//')
  BEFORE=$(awk '$1=="Before:"{print $2}' ${FILE})
  AFTER=$(awk '$1=="After:" {print $2}' ${FILE})
  echo -e "${ID}\t${BEFORE}\t${AFTER}"
done | awk 'BEGIN{print "sample\tsites\tdenoised_sites"}{print}' \
  > "${F}.biallelic_sites_count.tsv"
```

Combine raw and denoised `nQuire` results into a single line. 

NOTE: the results are "raw" then "denoised".

```bash
cat "${F}.biallelic_sites_count.tsv" \
  | ./add_value_to_table.py \
    -a <(grep -v 'denoise' "$F.ploidy" | sed -e 's/.bin//') \
  | ./add_value_to_table.py \
    -a <(grep 'free\|denoise' "$F.ploidy" | sed -e 's/.denoise.bin//') \
  > "$F.nQuire_results.tsv"
```

Download results from server.

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0022_Coral_Genotype_Analysis/03_Analysis/2022-02-05/Pocillopora_acuta/02_ploidy_analysis/nQuire_per_scaffold/"
rsync -zarv --prune-empty-dirs --relative \
 --include="*/" \
 --include="*.sh*" --include="*.txt" \
 --include="*.tsv" \
 --include="*.sitesProp.gz" \
 --exclude="*" \
 ${WD}/./* \
 . --dry-run
```

Create sample list for plotting (i.e., list of SNP allele proportion files and plot titles to use)

```bash
awk 'NR>1' Pocillopora_acuta_KBHIv2.assembly.fasta.1Mbp_scaffold_features.bed.lrdmodel.nQuire_results.tsv | sort -k19,19 -k1,1 | awk '{split($1, a, "-"); print $1".denoise.bin.coverage.sitesProp.gz\t"a[2]" ("$19")"}' | sed -e 's/.coordsorted-/.coordsorted.lowerProp-/' > samples_list.txt 
```

Use the `plot_allele_depths.Rmd` script to plot the results as a multipage PDF.

