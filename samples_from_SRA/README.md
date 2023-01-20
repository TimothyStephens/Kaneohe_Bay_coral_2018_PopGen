# Download Scleractinia samples list from SRA

Download list of Scleractinia coral sequencing samples from SRA. Filter results so we only have a list of samples that can be used for genotyping. 

## Setup analysis directory

Link to custom scripts for parsing text files.

```bash
ln -s ../../../02_Scripts/grepf_column.py
```

Setup bash environment.

```bash
conda activate py27
```

## Download and sort SRA Runs

Download all SRA runs that are from Scleractinia species. NCBI SRA search term: `Scleractinia[Organism] `

Downloaded 18,774 entries on the 15th June 2022 (file called: `SraRunInfo-Scleractinia-20220615.csv`). For some reason there are 19,050 entires in downloaded `*.cvs` file; won't matter after filtering.

Columns in file:

| Index | SraRunInfo column     |
| ----- | --------------------- |
| 1     | Run                   |
| 2     | ReleaseDate           |
| 3     | LoadDate              |
| 4     | spots                 |
| 5     | bases                 |
| 6     | spots_with_mates      |
| 7     | avgLength             |
| 6     | size_MB               |
| 9     | AssemblyName          |
| 10    | download_path         |
| 11    | Experiment            |
| 12    | LibraryName           |
| 13    | LibraryStrategy       |
| 14    | LibrarySelection      |
| 15    | LibrarySource         |
| 16    | LibraryLayout         |
| 17    | InsertSize            |
| 18    | InsertDev             |
| 19    | Platform              |
| 20    | Model                 |
| 21    | SRAStudy              |
| 22    | BioProject            |
| 23    | Study_Pubmed_id       |
| 24    | ProjectID             |
| 25    | Sample                |
| 26    | BioSample             |
| 27    | SampleType            |
| 28    | TaxID                 |
| 29    | ScientificName        |
| 30    | SampleName            |
| 31    | g1k_pop_code          |
| 32    | source                |
| 33    | g1k_analysis_group    |
| 34    | Subject_ID            |
| 35    | Sex                   |
| 36    | Disease               |
| 37    | Tumor                 |
| 38    | Affection_Status      |
| 39    | Analyte_Type          |
| 40    | Histological_Type     |
| 41    | Body_Site             |
| 42    | CenterName            |
| 43    | Submission            |
| 44    | dbgap_study_accession |
| 45    | Consent               |
| 46    | RunHash               |
| 47    | ReadHash              |

**LibraryStrategy**

| Tag                                     | Explanation                                                  |      |
| --------------------------------------- | ------------------------------------------------------------ | ---- |
| WGA                                     | Random  sequencing of the whole genome following non-pcr amplification |      |
| WGS                                     | Random  sequencing of the whole genome                       |      |
| WXS                                     | Random  sequencing of exonic regions selected from the genome |      |
| RNA-Seq                                 | Random  sequencing of whole transcriptome                    |      |
| miRNA-Seq                               | Random  sequencing of small miRNAs                           |      |
| WCS                                     | Random  sequencing of a whole chromosome or other replicon isolated from a genome |      |
| CLONE                                   | Genomic  clone based (hierarchical) sequencing               |      |
| POOLCLONE                               | Shotgun of  pooled clones (usually BACs and Fosmids)         |      |
| AMPLICON                                | Sequencing  of overlapping or distinct PCR or RT-PCR products |      |
| CLONEEND                                | Clone end  (5', 3', or both) sequencing                      |      |
| FINISHING                               | Sequencing  intended to finish (close) gaps in existing coverage |      |
| ChIP-Seq                                | Direct  sequencing of chromatin immunoprecipitates           |      |
| MNase-Seq                               | Direct  sequencing following MNase digestion                 |      |
| DNase-Hypersensitivity                  | Sequencing  of hypersensitive sites, or segments of open chromatin that are more readily  cleaved by DNaseI |      |
| Bisulfite-Seq                           | Sequencing  following treatment of DNA with bisulfite to convert cytosine residues to  uracil depending on methylation status |      |
| Tn-Seq                                  | Sequencing  from transposon insertion sites                  |      |
| EST                                     | Single pass  sequencing of cDNA templates                    |      |
| FL-cDNA                                 | Full-length  sequencing of cDNA templates                    |      |
| CTS                                     | Concatenated  Tag Sequencing                                 |      |
| MRE-Seq                                 | Methylation-Sensitive  Restriction Enzyme Sequencing strategy |      |
| MeDIP-Seq                               | Methylated  DNA Immunoprecipitation Sequencing strategy      |      |
| MBD-Seq                                 | Direct  sequencing of methylated fractions sequencing strategy |      |
| Synthetic-Long-Read                     |                                                              |      |
| ATAC-seq                                | Assay for  Transposase-Accessible Chromatin (ATAC) strategy is used to study genome-wide  chromatin accessibility. alternative method to DNase-seq that uses an  engineered Tn5 transposase to cleave DNA and to integrate primer DNA  sequences into the cleaved genomic DNA |      |
| ChIA-PET                                | Direct  sequencing of proximity-ligated chromatin immunoprecipitates. |      |
| FAIRE-seq                               | Formaldehyde  Assisted Isolation of Regulatory Elements. reveals regions of open chromatin |      |
| Hi-C                                    | Chromosome  Conformation Capture technique where a biotin-labeled nucleotide is  incorporated at the ligation junction, enabling selective purification of  chimeric DNA ligation junctions followed by deep sequencing |      |
| ncRNA-Seq                               | Capture of  other non-coding RNA types, including post-translation modification types  such as snRNA (small nuclear RNA) or snoRNA (small nucleolar RNA), or  expression regulation types such as siRNA (small interfering RNA) or  piRNA/piwi/RNA (piwi-interacting RNA). |      |
| RAD-Seq                                 |                                                              |      |
| RIP-Seq                                 | Direct  sequencing of RNA immunoprecipitates (includes CLIP-Seq, HITS-CLIP and  PAR-CLIP). |      |
| SELEX                                   | Systematic  Evolution of Ligands by EXponential enrichment   |      |
| ssRNA-seq                               | strand-specific  RNA sequencing                              |      |
| Targeted-Capture                        |                                                              |      |
| Tethered Chromatin Conformation Capture |                                                              |      |
| OTHER                                   | Library  strategy not listed (please include additional info in the “design  description”) |      |

**LibrarySelection**

| Tag                                     | Explanation                                                  |      |
| --------------------------------------- | ------------------------------------------------------------ | ---- |
| RANDOM                                  | Random  selection by shearing or other method                |      |
| PCR                                     | Source  material was selected by designed primers            |      |
| RANDOM PCR                              | Source  material was selected by randomly generated primers  |      |
| RT-PCR                                  | Source  material was selected by reverse transcription PCR   |      |
| HMPR                                    | Hypo-methylated  partial restriction digest                  |      |
| MF                                      | Methyl  Filtrated                                            |      |
| CF-S                                    | Cot-filtered  single/low-copy genomic DNA                    |      |
| CF-M                                    | Cot-filtered  moderately repetitive genomic DNA              |      |
| CF-H                                    | Cot-filtered  highly repetitive genomic DNA                  |      |
| CF-T                                    | Cot-filtered  theoretical single-copy genomic DNA            |      |
| MDA                                     | Multiple  displacement amplification                         |      |
| MSLL                                    | Methylation  Spanning Linking Library                        |      |
| cDNA                                    | complementary  DNA                                           |      |
| ChIP                                    | Chromatin  immunoprecipitation                               |      |
| MNase                                   | Micrococcal  Nuclease (MNase) digestion                      |      |
| DNAse                                   | Deoxyribonuclease  (MNase) digestion                         |      |
| Hybrid Selection                        | Selection  by hybridization in array or solution             |      |
| Reduced  Representation                 | Reproducible  genomic subsets, often generated by restriction fragment size selection,  containing a manageable number of loci to facilitate re-sampling |      |
| Restriction Digest                      | DNA  fractionation using restriction enzymes                 |      |
| 5-methylcytidine  antibody              | Selection  of methylated DNA fragments using an antibody raised against 5-methylcytosine  or 5-methylcytidine (m5C) |      |
| MBD2 protein  methyl-CpG binding domain | Enrichment  by methyl-CpG binding domain                     |      |
| CAGE                                    | Cap-analysis  gene expression                                |      |
| RACE                                    | Rapid  Amplification of cDNA Ends                            |      |
| size fractionation                      | Physical  selection of size appropriate targets              |      |
| Padlock probes  capture method          | Circularized  oligonucleotide probes                         |      |
| other                                   | Other  library enrichment, screening, or selection process (please include  additional info in the “design description”) |      |
| unspecified                             | Library  enrichment, screening, or selection is not specified (please include  additional info in the “design description”) |      |
| cDNA_oligo_dT                           |                                                              |      |
| cDNA_randomPriming                      |                                                              |      |
| Inverse rRNA                            | depletion of  ribosomal RNA by oligo hybridization.          |      |
| Oligo-dT                                | enrichment of  messenger RNA (mRNA) by hybridization to Oligo-dT. |      |
| PolyA                                   | PolyA  selection or enrichment for messenger RNA (mRNA); should replace cDNA  enumeration. |      |
| repeat  fractionation                   | Selection for  less repetitive (and more gene rich) sequence through Cot filtration (CF) or  other fractionation techniques based on DNA kinetics. |      |



Filter Runs by selecting:

- Platform == "ILLUMINA"
- LibraryStrategy == "WGA" || LibraryStrategy == "WGS" || LibraryStrategy == "RNA-Seq"

```bash
awk -F',' '{ if(NR==1) {print $0} else { if($19=="ILLUMINA" && ($13 == "WGA" || $13 == "WGS" || $13 == "RNA-Seq")){print $0} } }' SraRunInfo-Scleractinia-20220615.csv > SraRunInfo-Scleractinia-20220615-filtered.csv
```

We end up with **9,590** samples.

Print the number of samples we have from each genus.

```bash
cat SraRunInfo-Scleractinia-20220615-filtered.csv | awk -F',' '{print $29}' | awk '{print $1}' | sort | uniq -c | sort -k1,1n | awk '{print $2"\t"$1}'
```

---

Split filtered samples into subsets based on species; group species by genus. Remove samples that are on SRA that were generated by this study (listed in `samples_from_this_study.txt`).



Link to text parsing scripts.

```bash
ln -s ../../../02_Scripts/grepf_column.py
ln -s ../../../02_Scripts/add_value_to_table.py
```



### ***Montipora***

Extract just *Montipora* samples.

```bash
awk -F',' '{ if(NR==1) {print $0} else if ($29~"^Montipora") {print $0} }' \
    SraRunInfo-Scleractinia-20220615-filtered.csv \
  > Montipora_SraRunInfo.csv
cat Montipora_SraRunInfo.csv \
  | grep -v -f samples_from_this_study.txt \
  > Montipora_SraRunInfo.notThisStudy.csv
```

Create directory `BioSample`. From here we will extract and download the Biosample metadata associated with each SRR sample that we want to run.

**Extract the SRA metadata for each SRR ID and also fetch the sample location from the associated BioSample using `efetch`.**

Download BioSample information for each SRR ID in list.

```bash
mkdir -p BioSample
cat Montipora_SraRunInfo.notThisStudy.csv \
| awk -F',' '$1!="Run" {print $26}' | sort | uniq | \
  while read BIOSAMPLE; 
  do 
    echo $BIOSAMPLE; 
    efetch -db biosample -id $BIOSAMPLE > BioSamples/$BIOSAMPLE.info; 
  done
```

Extract just the "geographic location" location information from the BioSample metadata we just downloaded.

```bash
for FILE in BioSamples/SAM*.info;
do 
  BIOSAMPLE=${FILE%*.info}
  ORG=$(grep '^Organism:' $FILE | sed -e 's/.*: //' -e 's/"//g')
  SAM=$(grep '/sample name=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  ISO=$(grep '/isolate=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  ISS=$(grep '/isolation source=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  COL=$(grep '/collection date=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  LOC=$(grep '/geographic location=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  TIS=$(grep '/tissue=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  DEV=$(grep '/development stage=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  echo -e "$BIOSAMPLE\t$ORG\t$SAM\t$ISO\t$ISS\t$COL\t$LOC\t$TIS\t$DEV"
done | awk 'BEGIN{print "BioSample\torganism\tsample_name\tisolate\tisolation_source\tcollection_date\tgeographic_location\ttissue\tdevelopment_stage"}{print}' \
  | sed -e 's@BioSamples/@@' \
  > Montipora_biosamples.info
```

Join SRA RunInfo metadata file with BioSample metadata list.

```bash
sed -e 's/,/\t/g' Montipora_SraRunInfo.notThisStudy.csv \
| ./add_value_to_table.py -c 26 \
    -a Montipora_biosamples.info \
    -o Montipora_SraRunInfo_BioSample.notThisStudy.tsv
```



### ***Pocillopora***

Extract just *Pocillopora* samples

```bash
awk -F',' '{ if(NR==1) {print $0} else if ($29~"^Pocillopora") {print $0} }' \
    SraRunInfo-Scleractinia-20220615-filtered.csv \
  > Pocillopora_SraRunInfo.csv
cat Pocillopora_SraRunInfo.csv \
  | grep -v -f samples_from_this_study.txt \
  > Pocillopora_SraRunInfo.notThisStudy.csv
```

Create directory `BioSamples`. From here we will extract and download the Biosample metadata associated with each SRR sample that we want to run.

**Extract the SRA metadata for each SRR ID and also fetch the sample location from the associated BioSample using `efetch`.**

Download BioSample information for each SRR ID in list.

```bash
cat Pocillopora_SraRunInfo.notThisStudy.csv \
| awk -F',' '$1!="Run" {print $26}' | sort | uniq | \
  while read BIOSAMPLE; 
  do 
    echo $BIOSAMPLE; 
    efetch -db biosample -id $BIOSAMPLE > BioSamples/$BIOSAMPLE.info; 
  done
```

Extract just the "geographic location" location information from the BioSample metadata we just downloaded.

```bash
for FILE in BioSamples/SAM*.info;
do 
  BIOSAMPLE=${FILE%*.info}
  ORG=$(grep '^Organism:' $FILE | sed -e 's/.*: //' -e 's/"//g')
  SAM=$(grep '/sample name=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  ISO=$(grep '/isolate=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  ISS=$(grep '/isolation source=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  COL=$(grep '/collection date=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  LOC=$(grep '/geographic location=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  TIS=$(grep '/tissue=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  DEV=$(grep '/development stage=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  echo -e "$BIOSAMPLE\t$ORG\t$SAM\t$ISO\t$ISS\t$COL\t$LOC\t$TIS\t$DEV"
done | awk 'BEGIN{print "BioSample\torganism\tsample_name\tisolate\tisolation_source\tcollection_date\tgeographic_location\ttissue\tdevelopment_stage"}{print}' \
  | sed -e 's@BioSamples/@@' \
  > Pocillopora_biosamples.info
```

Join SRA RunInfo metadata file with BioSample metadata list.

```bash
sed -e 's/,/\t/g' Pocillopora_SraRunInfo.notThisStudy.csv \
| ./add_value_to_table.py -c 26 \
    -a Pocillopora_biosamples.info \
    -o Pocillopora_SraRunInfo_BioSample.notThisStudy.tsv
```

---

---

---

---

---

Modified close from above to grab out extra info for identifying how each sample relates to each other. For example, for *P. acute* it was helpful to have the lat/long coords to figure out which reef they were from, and for *M. capitata* the  sample name identifier in the file actually gave more info then the sample name field, so extracting that info helped with those samples.

```bash
for FILE in BioSamples/SAM*.info;
do 
  BIOSAMPLE=${FILE%*.info}
  ORG=$(grep '^Organism:' $FILE | sed -e 's/.*: //' -e 's/"//g')
  SAM=$(grep '/sample name=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  ISO=$(grep '/isolate=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  ISS=$(grep '/isolation source=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  COL=$(grep '/collection date=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  LOC=$(grep '/geographic location=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  TIS=$(grep '/tissue=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  DEV=$(grep '/development stage=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  TRE=$(grep '/treatment=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  LAT=$(grep '/latitude and longitude=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  echo -e "$BIOSAMPLE\t$ORG\t$SAM\t$ISO\t$ISS\t$COL\t$LOC\t$TIS\t$DEV\t$TRE\t$LAT"
done | awk 'BEGIN{print "BioSample\torganism\tsample_name\tisolate\tisolation_source\tcollection_date\tgeographic_location\ttissue\tdevelopment_stage\ttreatment\tlat_long"}{print}' \
  | sed -e 's@BioSamples/@@' \
  > Pocillopora_biosamples.info_v2


sed -e 's/,/\t/g' Pocillopora_SraRunInfo.notThisStudy.csv \
| ./add_value_to_table.py -c 26 \
    -a Pocillopora_biosamples.info_v2 \
    -o Pocillopora_SraRunInfo_BioSample.notThisStudy_v2.tsv
```

```bash
for FILE in BioSamples/SAM*.info;
do 
  BIOSAMPLE=${FILE%*.info}
  ORG=$(grep '^Organism:' $FILE | sed -e 's/.*: //' -e 's/"//g')
  SAM=$(grep '/sample name=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  ISO=$(grep '/isolate=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  ISS=$(grep '/isolation source=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  COL=$(grep '/collection date=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  LOC=$(grep '/geographic location=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  TIS=$(grep '/tissue=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  DEV=$(grep '/development stage=' $FILE | sed -e 's/.*=//' -e 's/"//g')
  SM2=$(grep 'Sample name:' $FILE | sed -e 's/.*Sample name: //' -e 's/;.*//g')
  echo -e "$BIOSAMPLE\t$ORG\t$SAM\t$ISO\t$ISS\t$COL\t$LOC\t$TIS\t$DEV\t$SM2"
done | awk 'BEGIN{print "BioSample\torganism\tsample_name\tisolate\tisolation_source\tcollection_date\tgeographic_location\ttissue\tdevelopment_stage\tsample_name_2"}{print}' \
  | sed -e 's@BioSamples/@@' \
  > Montipora_biosamples.info_v2


sed -e 's/,/\t/g' Montipora_SraRunInfo.notThisStudy.csv \
| ./add_value_to_table.py -c 26 \
    -a Montipora_biosamples.info_v2 \
    -o Montipora_SraRunInfo_BioSample.notThisStudy_v2.tsv
```

