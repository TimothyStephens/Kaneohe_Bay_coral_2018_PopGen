# Salmon read mapping

Align RNA-seq reads against CDS of predicted protein-coding genes using `salmon`.

## Setup analysis directory

Link to CDS and genome fasta files.

```bash
ln -s ../00_databases/Pocillopora_acuta_KBHIv2.genes.cds.fna
ln -s ../00_databases/Pocillopora_acuta_KBHIv2.assembly.fasta
```

Setup bash environment.

```bash
conda activate py27

export PATH="$PATH:/home/timothy/programs/salmon-1.6.0_linux_x86_64/bin"
```

## Prepare decoy sequences

Prepare the decoy sequences. In this case we are using the whole genome as a decoy (See https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode for more information). All we have to do is extract the scaffold names (so salmon knows which sequences are decoys, and combine the genome+CDS sequences ready for `salmon index`)

```bash
grep "^>" Pocillopora_acuta_KBHIv2.assembly.fasta | sed -e 's/>//g' > decoys.txt

cat \
  Pocillopora_acuta_KBHIv2.genes.cds.fna \
  Pocillopora_acuta_KBHIv2.assembly.fasta \
  > Pocillopora_acuta_KBHIv2.gentrome.fa

gzip Pocillopora_acuta_KBHIv2.gentrome.fa
```

## Run `salmon`

Index the "gentrome" file using `salmon index`

```bash
./run_salmon_index.sh
```

> 471 transcripts were removed because they were duplicates (i.e., identical to other transcripts in the file)

Run `salmon quant`

```bash
./run_salmon_quant.sh
```

Check for any unexpected errors or warning.

```bash
grep -i 'error\|warn\|fault' *.log.*
```

## Combine `salmon` results from each sample

Use `salmon quant merge` to combine the results from across all samples into a big matrix. 

Generate matrices of TPM and "numread" values from across the samples. Compress then once done.

```bash
F="Pocillopora_acuta_KBHIv2.gentrome.fa.gz"
DIR_LIST=$(awk -F'\t' '{printf " %s", $2}' samples_Pacuta.txt)

salmon quantmerge --column numreads --output $F.salmon.numreads.matrix --quants $DIR_LIST
salmon quantmerge --column tpm --output $F.salmon.tpm.matrix --quants $DIR_LIST
salmon quantmerge --column elen --output $F.salmon.elen.matrix --quants $DIR_LIST

gzip $F.salmon.numreads.matrix
gzip $F.salmon.tpm.matrix
gzip $F.salmon.elen.matrix
```

## `salmon` run stats

```bash
echo -e "sample_id\tSalmon_mapping_rate" > stats.txt
for D in $(awk -F'\t' '{printf " %s", $2}' samples_Pacuta.txt); 
do 
  echo -e "$D\t"$(grep 'Mapping rate' $D/logs/salmon_quant.log | awk '{print $8}')
done >> stats.txt
```

## Cleanup analysis directory

Tar-zip the salmon analysis directory created for each sample since we don't need it anymore. 

```bash
for D in $(awk -F'\t' '{printf " %s", $2}' samples_Pacuta.txt); do echo "$D"; tar -zcf $D.tar.gz $D && rm -r $D; done
```

## Results

### Download results from server

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0022_Coral_Genotype_Analysis/03_Analysis/2022-02-05/Pocillopora_acuta/04_salmon_readmapping"
rsync -zarv --prune-empty-dirs --relative \
 --include="*/" --include="*.sh*" --include="*.py" \
 --include="*.matrix.gz" --include="duplicate_clusters.tsv" \
 --exclude="*" \
 ${WD}/./ \
 . --dry-run
```

