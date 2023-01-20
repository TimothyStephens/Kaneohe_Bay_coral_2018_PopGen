# Extract bi-allelic sites using `nQuire view` from *Pocillopora* samples downloaded from SRA

## Setup analysis directory

Link `bam` files.

```bash
while read ID;
do
  ln -s "../nQuire/${ID}.SplitNCigarReads.split.bam"
done < <(cut -f1 ../../../Pocillopora_selected_Runs.txt)
```

Setup bash environment.

```bash
conda activate py27
```

## Run `nQuire view`

Run `nQuire view` on denoised `bin` files that were created with a different (lower "min fraction of read support" [-f param]).

```bash
for BAM in *.bam;
do
  echo "./run_nQuire_extract_AP.sh $BAM"
done | parallel -j 10
```

Check that the jobs completed correct and no errors occurred.

```bash
grep -i 'error\|warn\|fault' *.log.*
grep 'ExitStatus:' *.log.* | awk '{print $2}' | sort | uniq -c
```

Last command should give **32** total.

## Download results from server

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0022_Coral_Genotype_Analysis/03_Analysis/2022-02-05/samples_from_SRA/Pocillopora/02_ploidy_analysis/extract_biallelic_sites/"
rsync -zarv --prune-empty-dirs --relative \
 --include="*/" \
 --include="*.sh*" --include="*.txt" \
 --include="*.sitesProp.gz" \
 --exclude="*" \
 ${WD}/./* \
 . --dry-run
```





