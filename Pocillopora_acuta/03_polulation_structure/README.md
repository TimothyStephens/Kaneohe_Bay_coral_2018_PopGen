# Use various approaches to explore the population of the *P. acuta* RNA-Seq samples collected from the bay.

## `vcftools --relatedness2`

Use `vcftools` to calculate relatedness (`--relatedness2`) stats between the samples.

Filter the `*.vcf.gz` file first using `vcftools` so that we have just sites with mean depth values (over all included individuals) greater >= 10 (`--min-meanDP 10`), no missing alleles across any of the samples (`--max-missing 1.0`), and only SNPs not indels (`--remove-indels`).

```bash
./run_vcftools.sh
```

Download results.

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0022_Coral_Genotype_Analysis/03_Analysis/2022-02-05/Pocillopora_acuta/03_polulation_structure/vcftools_relatedness2/"
rsync -zarv --prune-empty-dirs --relative \
  --include="*/" \
  --include="*.relatedness2" \
  --include="*.sh*" \
  --exclude="*" \
  ${WD}/./* \
  . --dry-run
```

Plot results using `plot_vcftools_relatedness2_results.Rmd`

# `vcf_clone_detect.py`

Use `vcf_clone_detect.py` to calculate the number of shared alleles between each pairwise combination of samples.

Run tool (using the below script) on the filtered variant sites produced by `vcftools --min-meanDP 10 --max-missing 1.0 --remove-indels`.

```bash
./run_vcf_clone_detect.sh
```

The `*.allelic_similarity.csv` file produced by the `run_vcf_clone_detect.sh` script lists the pairwise similairty of the SNPs between each sample. This file however only does this in a single direction (A -> B). We need it to list both directions (A -> B and B -> A; will have same values) and self similarity (A -> A and B -> B; will be 100) so that we can create a symmetric similarity matrix out of these results.

```bash
VCF="GVCFall.filtered.recode.vcf.gz"

echo -e "ind1\tind2\tmatch_perc" > ${VCF}.allelic_similarity.full.tsv
awk -F',' 'NR>1{print $1"\t"$2"\t"$7"\n"$2"\t"$1"\t"$7}' \
    ${VCF}.allelic_similarity.csv \
  >> ${VCF}.allelic_similarity.full.tsv
awk 'NR>1{print $1"\t"$1"\t100"}' \
    ../../samples_*.annotations.txt \
  >> ${VCF}.allelic_similarity.full.tsv
```

Download results.

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0022_Coral_Genotype_Analysis/03_Analysis/2022-02-05/Pocillopora_acuta/03_polulation_structure/vcf_clone_detect/"
rsync -zarv --prune-empty-dirs --relative \
  --include="*/" \
  --include="*.allelic_similarity.*" \
  --include="*.sh*" \
  --exclude="*" \
  ${WD}/./* \
  . --dry-run
```

Plot results using `plot_vcf_clone_detect_results.Rmd`

## `ANGSD`

Use `ANGSD` as well as other associated tools to analyze the *P. acuta* RNA-Seq samples generated in this study. 

Get list of input `BAM` files to analyze.

```bash
ls --color=none -1 ../../01_RNAseq_samples_readmapping/*.SplitNCigarReads.split.bam > bam.filelist
sed -e 's@../../01_RNAseq_samples_readmapping/@@' -e 's@.SplitNCigarReads.split.bam@@' bam.filelist > bam.filelist.labels
```

Run `ngsRelate` on MAFs and GLF files produced by `ANGSD`.

```bash
./run_02-PCAngsd_01-ANGSD.sh
./run_02-PCAngsd_02-PCAngsd.sh
```

Download results.

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0022_Coral_Genotype_Analysis/03_Analysis/2022-02-05/Pocillopora_acuta/03_polulation_structure/ANGSD/"
rsync -zarv --prune-empty-dirs --relative \
  --include="*/" \
  --include="bam.filelist*" \
  --include="PCAngsd.angsd.beagle.gz.*" \
  --include="NgsRelate.results.tsv" \
  --include="*.sh*" \
  --exclude="*" \
  ${WD}/./* \
  . --dry-run
```

Plot results using `plot_ANGSD_results.Rmd`
