# another thing for a classification comparison:
Divoll noted the R.bold function to accesss the BOLD API directly.
  - See: https://github.com/tdivoll/Bat-Diet-Metabarcoding/blob/master/BatDietWorkflow5.md#3-query-bold-for-otu-matches
Would be a useful comparison among methods.
```
library('bold')
## data import is just a 2-column fasta...
mydata <- {import a 2 column dataset where col1 == asvID, col2==nucleotide string}
mydata2 <- as.list(setNames(mydata$seqs, mydata$seqID))  
output <- bold_identify(sequences = mydata2, db = "COX1", response=FALSE)
```

# Motivations
We were interested in evaluating how certain filtering decisions affected the representation of taxa in the database we construct. Some questions pertained to the raw sequences obtained prior to filtering, while others pertained to our filtered and dereplicated dataset.

## Alternative database setup
We also decided to explore similar questions using a pair of databases available online formatted by other labs:

- Jon Palmer's [amptk](https://amptk.readthedocs.io/en/latest/). Following installation, you'll need to install the database and then convert it to a fasta file. Note that in the case of Jon's database, this is strictly the BOLD records pulled using his pipeline - the actual classification incorporates additional UTAX and SINTAX databases which are trained from these records. The point of using the BOLD only records is simply to understand how different filtering parameters from the same resource (BOLD) can retain distinct taxa.
```
## assumes Conda install
amptk install -i COI
cd ~/.conda/envs/amptk/lib/python3.6/site-packages/amptk/DB
vsearch --udb2fasta COI.udb --output amptk.coi.fasta
```

- Terri Porter's Eukaryote CO1 mtDNA [datasets](https://github.com/terrimporter/CO1Classifier). We can download this fasta file directly from their v3.2 release [here](https://github.com/terrimporter/CO1Classifier/releases/tag/v3.2-ref).
```
wget https://github.com/terrimporter/CO1Classifier/releases/download/v3.2-ref/CO1v3_2_training.tar.gz
tar xzf CO1v3_2_training.tar.gz
```

Note that both the Palmer and Porter databases contain non-arthropod COI sequences. In the case of the Palmer dataset, it contains chordate and athropod recrds, while the Porter dataset contains chordate, arthropod, and microeukaryote COI sequences. For the purposes of our evaluations, we restricted our comparisons to records identified as Arthropod and ignored all non-arthropod records. In addition, the Porter dataset is _not_ dereplicated, thus we created both a dereplicated and non-dereplicated dataset for different downstream analyses (see `database_composition.R` script). The Palmer data is already dereplicated.

Filtering the Palmer dataset:
```
## see `database_construction.md` for details on virtual environment setup
source activate arrR
seqkit grep -r -p 'p:Arthropoda' amptk.coi.fasta -w 0 > palmer.arthCOI.fasta
```
This reduces the initial **1,617,885** COI records down to **1,565,831** sequences.

Filtering and dereplicating the Porter dataset:
```
grep ';Arthropoda;' -A 1 mytrainseq.fasta |  sed '/--/d' > porter.arthCOI.fasta
vsearch --derep_fulllength porter.arthCOI.fasta --output porter.derep.fasta --relabel_keep --threads 4 --fasta_width 0 --notrunclabels
```
Filtering for Arthropod recoreds reduces the initial **1,280,577** COI records down to **883,979** sequences. Dereplication further reduces these records down to **515,780** disintct sequences.

## Database questions to address
Using the arthropod-selected Palmer and Porter databases, the raw BOLD COI records we retained with the `bold_datapull.R` script, and  the resulting filtered dataset (see `database_construction.md`), we addressed the following questions:
How do these datasets compare with respect to
  - distribution of sequence length?
  - number of unique taxa per level
  - missingness/completeness?

Specific to our BOLD datasets, we evaluated:
- Does specifying records be obtained/associated with a certain country/region effect the representation of taxonomic information?
- What institutions are providing most COI sequences among BOLD arthropod data?

Specific to our filtered BOLD dataset:
- How does clustering a database effect the representation of taxonomic information going into classification?

# Analyses
## Distribution of sequence lengths:
We evaluated the distribution of sequence lengths among each of the four datasets. Note that for the Porter dataset we only used the dereplicated data so as to not double count identical sequences.
```
seqkit fx2tab --length --name --header-line palmer.arthCOI.fasta | cut -f 4 | gzip --best > palmer.lengths.txt.gz
seqkit fx2tab --length --name --header-line porter.derep.fasta | cut -f 4 | gzip --best > porter.lengths.txt.gz
seqkit fx2tab --length --name --header-line boldCOI.derep.fasta | cut -f 4 | gzip --best > derep.lengths.txt.gz
zcat boldCustom.allArth.seqNtaxa.csv.gz | cut -f 3 -d ',' | awk 'NR > 1 { print length }' | gzip --best > raw.lengths.txt.gz
```
These files were used to generate the plot with the R script `database_lengths.R`

We find that the plurality of sequences are represented by a single length of 658bp for all databases, though the proportion of these ranges from the lowest frequency in the Palmer database (24.2%) to the highest in the Raw dataset (36.5%).
**Why is this the case? Likely because of the dominant arthropodo records having full length COI at that length**
The next most frequently observed length for most datasets ranges is observed for 588 bp amplicons, though there is a substantial range of values near this length that are frequently observed.
**What is significant of 588?**
 However, the Palmer dataset contains a unique outlier among our databases having a ~ 13% of observed sequences at lengths of 318-319bp.
**Why are these in his and not other dbs?** For one, Porter requires a minimum at 500bp, so we won't see theirs.


## Distribution of taxonomic composition and taxonomic completeness
Proportions of unique taxa represented from Class to Species levels were assessed for each of the four databases. In addition, taxonomic completeness was evaluated by determining the number of instances in which no data was present at a given taxonomic Level (Class through Species). Finally, we also evaluated the most abundant (top 20) taxa for Class through Genus Levels in the three dereplicated datasets (Palmer, Porter, and ours), though we also retained the full Porter dataset and our own non-dereplicated dataset to illustrate how some records are represented hundreds or thousands of times in raw data.
We included records that contained any information in the Species designation including an arbitrary "sp." marker, thus these results represent a liberal estimate of Species information. Of note, because the `bold_datapull.R` script used to download data specified either Class or Order names to pull data, we are potentially discarding records within BOLD that contain a COI marker but no Phylum or Class information - we made no attempts at determining how pervasive this level of missingness could be. In addition, the Porter method queried NCBI with a Perl script that specified that "species" [RANK] be included for the Arthropod records, so there should not be any taxonomic information missing from this dataset (however, by requiring species-level information, certain taxa may no longer be represented at higher levels). These discrepancies are evaluated in the `classification_analyses.md` document.

We pulled the taxa strings from the Palmer and Porter (dereplicated and non-dereplicated) datasets. Only the dereplicated Porter dataset was used for analysis of the taxa, :
```
cat palmer.arthCOI.fasta | grep '^>' | cut -d ',' -f2- | awk -F ',' '{print $2,$3,$4,$5,$6}' OFS=',' | gzip --best > palmer.taxa.txt.gz
cat porter.arthCOI.fasta | grep '^>' | cut -f5- -d ';' | tr ';' ',' | gzip --best > porter.all.taxa.txt.gz
cat porter.derep.fasta | grep '^>' | cut -f5- -d ';' | tr ';' ',' | gzip --best > porter.derep.taxa.txt.gz
```

We modified our dereplicated and raw datasets, `boldCOI.derep.txt` and `boldCustom.allArth.meta.txt.gz`, to reduce disk space by pulling just the taxa strings from Class through Species:
```
zcat boldCOI.derep.txt.gz | cut -d ';' -f 3- | gzip --best > derep.taxa.txt.gz
zcat boldCustom.allArth.meta.txt.gz | cut -d ';' -f 8-13 | gzip --best > raw.taxa.txt.gz
```
The `*lengths.txt` and `*taxa.txt.gz` files were used to generate the tables with the R script `database_composition.R`.

## Impact of clustering on taxonomic completeness
In addition to the comparisons evaluating taxonomic completeness based on the database being evaluated, we analyzed how clustering a dataset can effect the remaining taxa. We clustered the `boldCOI.derep.fasta` file at three values: 99%, 97%, and 95% identity using Vsearch as described in the `database_construction.md` document.

The remaining sequenceID records for each clustered fasta were obtained from each fasta, then imported into the `database_clustering.R` script to create the plot.

```
cat boldCOI.clust99.fasta | grep "^>" | sed 's/>//' | gzip --best > clust99.names.txt.gz
cat boldCOI.clust99.fasta | grep "^>" | sed 's/>//' | gzip --best > clust97.names.txt.gz
cat boldCOI.clust99.fasta | grep "^>" | sed 's/>//' | gzip --best > clust95.names.txt.gz
```

Taxonomic info for each record in the `clust*.gz` files were subset from the `boldCOI.derep.txt.gz` file.

## Impact of geographic specificity
We restricted these analyses to our own datasets - the raw and dereplicated records. The `boldCustom.allArth.meta.txt` served as the raw input, with the Sequence ID's present from the `boldCOI.derep.txt` file acting as a filter for comparison. To understand the impact on geographic specificity on taxonomic composition of a database, we filtered our BOLD records by requiring that the `country` column in the BOLD dataset contain either `United States` or `Canada` records. We then investigated how this effects the taxonomic composition at the Class and Order levels, as well as the proportional representation of sequence records from the top 10 sources listed in the `institution_storing` field of the BOLD dataset.

The **name** R script was used to genrate the associated plots.
**what is the name of the R script**
**what is the name of the R script**
**what is the name of the R script**
**what is the name of the R script**
```
countrylist <- c("United States", "Canada")
x.meta %>%
  filter(country %in% countrylist) %>%
  group_by(institution_storing) %>%
  summarise(nRecords = n()) %>%
  mutate(pRecords = nRecords / sum(nRecords)) %>%
  arrange(., -nRecords) -> tmp_inst_df
```
We find that the vast majority of data comes from one place: Centre for Biodiversity Genomics. They possess about 88% of all records attributed to an institution (269,875 of 305,440).

However, if we didn't require any particular filtered country, the proportion of data from that institution is comparably smaller, despite still containing the most records per institution.
```
x.meta %>%
  #filter(country %in% countrylist) %>%
  group_by(institution_storing) %>%
  summarise(nRecords = n()) %>%
  mutate(pRecords = nRecords / sum(nRecords)) %>%
  arrange(., -nRecords) -> tmp_inst_df
```
We find that about 10% of these records come from GenBank!

1. If we filter by just Genbank vs. CBG, what fraction of data are we retaining overall?
2. If that's okay, then what Orders represent the largest fraction?
3. Of those fractions, how are the Genbank/CBG records distributed?