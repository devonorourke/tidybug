# Motivations
I was interested in evaluating how certain filtering decisions used to generate a COI database affected the representation of taxa available to classify representative sequences. In other words, how do things like the source of COI sequences, the filtering strategies, or the clustering strategies affect the subsequent references available?  

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

Note that both the Palmer and Porter databases contain non-arthropod COI sequences. In the case of the Palmer dataset, it contains chordate and arthropod records, while the Porter dataset contains chordate, arthropod, and microeukaryote COI sequences. For the purposes of our evaluations, we restricted our comparisons to records identified as Arthropod and ignored all non-arthropod records. In addition, the Porter dataset is _not_ dereplicated, thus we created both a dereplicated and non-dereplicated dataset for different downstream analyses (see `database_composition.R` script). The Palmer data is already dereplicated.

### Palmer database filtering
Filtering the Palmer dataset:
```
## see `database_construction.md` for details on virtual environment setup
source activate arrR
seqkit grep -r -p 'p:Arthropoda' amptk.coi.fasta -w 0 > palmer.arthCOI.fasta
```
This reduces the initial **1,617,885** COI records down to **1,565,831** sequences.

### Porter database filtering
This takes a bit more work: first, remove all the non-arthropod records:
```
grep ';Arthropoda;' -A 1 mytrainseq.fasta |  sed '/--/d' > porter.arthCOI.fasta
```

Next, generate a QIIME-styled OTU table for all arthropod records:
```
pick_otus.py -i porter.arthCOI.fasta -o porter_pid100_otus --similarity 1.0 --threads 24
```

Create a fasta file with the headers stripped of any taxonomy information:
```
cat porter.arthCOI.fasta | sed '/^>/s/\s.*$//' > porter.arthCOI_notax.fasta.tmp
```

then create a taxonomy mapping file:
```
cat porter.arthCOI.fasta | grep '^>' | sed 's/^>//' | cut -d ' ' -f 1 > tmp.left
cat porter.arthCOI.fasta | grep '^>' | sed 's/^>//' | cut -d ' ' -f 2- | cut -d ';' -f 4- > tmp.right
paste -d '\t' tmp.left tmp.right > porter.taxa.tmp
rm tmp.left tmp.right
```

and run the LCA algorithm to create a taxonomy string that is appropriately classified:
```
python create_consensus_taxonomy.py porter.taxa.tmp porter.arthCOI_notax.fasta.tmp ./porter_pid100_otus/porter.arthCOI_otus.txt porter_outmap.txt
```

Finally, we create an updated fasta file by dereplicating the sequences with VSEARCH:
```
vsearch --derep_fulllength porter.arthCOI_notax.fasta.tmp --output porter.derep.fasta --relabel_keep --threads 4 --fasta_width 0 --notrunclabels
```

Filtering for Arthropod records reduces the initial **1,280,577** COI records down to **883,979** sequences. Dereplication further reduces these records down to **515,780** distinct sequences.

The remaining sequence records following dereplication that match the `porter_outmap.txt` taxonomy records are used for compositional comparisons.

## Database questions to address
Using the arthropod-selected Palmer and Porter databases, the raw BOLD COI records we retained with the [bold_datapull_arths.R](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/bold_datapull_arths.R) R script, and the resulting filtered dataset (see [database_construction.md](https://github.com/devonorourke/tidybug/blob/master/docs/database_construction.md)), we addressed the following questions:  

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
We evaluated the distribution of sequence lengths among each of the four datasets.
```
seqkit fx2tab --length --name --header-line palmer.arthCOI.fasta | cut -f 4 | gzip --best > palmer.lengths.txt.gz
seqkit fx2tab --length --name --header-line porter.derep.fasta | cut -f 4 | gzip --best > porter.lengths.txt.gz
seqkit fx2tab --length --name --header-line boldCOI.derep.fasta | cut -f 4 | gzip --best > derep.lengths.txt.gz
zcat boldCustom.allArth.seqNtaxa.csv.gz | cut -f 3 -d ',' | awk 'NR > 1 { print length }' | gzip --best > raw.lengths.txt.gz
```
These files were used to generate the plot with the R script [database_lengths.R](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/database_lengths.R).

We find that the plurality of sequences are represented by a single length of 658bp for all databases, though the proportion of these ranges from the lowest frequency in the Palmer database (24.2%) to the highest in the Raw dataset (36.5%).
**Why is this the case? Likely because of the dominant arthropod records having full length COI at that length**
The next most frequently observed length for most datasets ranges is observed for 588 bp amplicons, though there is a substantial range of values near this length that are frequently observed.
**What is significant of 588?**
 However, the Palmer dataset contains a unique outlier among our databases having a ~ 13% of observed sequences at lengths of 318-319bp.
**Why are these in his and not other dbs?** For one, Porter requires a minimum at 500bp, so we won't see theirs.


## Porportions of shared taxonomic names and differences in completeness of taxonomic records
Proportions of shared and distinct represented from Class to Species levels were assessed for each of the three databases. These are presented in [Figure 4B]((https://github.com/devonorourke/tidybug/blob/master/figures/figure4_upsetPlot.png)) of the manuscript. In addition, taxonomic completeness was evaluated by determining the number of instances in which no data was present at a given taxonomic Level (Class through Species) - these are represented in [Figure 4A]((https://github.com/devonorourke/tidybug/blob/master/figures/figure4_upsetPlot.png)).  

We included records that contained any information in the Species designation including an arbitrary "sp." marker, thus these results represent a liberal estimate of Species information. Of note, because the [bold_datapull_arths.R](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/bold_datapull_arths.R) script used to download data specified either Class or Order names to pull data, we are potentially discarding records within BOLD that contain a COI marker but no Phylum or Class information - we made no attempts at determining how pervasive this level of missingness could be. In addition, the Porter method queried NCBI with a Perl script that specified that "species" [RANK] be included for the Arthropod records, so there should not be any taxonomic information missing from this dataset (however, by requiring species-level information, certain taxa may no longer be represented at higher levels). These discrepancies are evaluated [Figure 4](https://github.com/devonorourke/tidybug/blob/master/figures/figure4_upsetPlot.png) of the manuscript; these data were generated using the [figure4_database_comparisons.R](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/figure4_database_comparisons.R) R script.

We performed a minor amount of data filtering to keep delimiter stucture and prefix for taxonomy levels consistent across datasets. In addition, we retained only Class through Species-level information (discarding any kingdom/phylum information because it was redundant for this comparison).

For Palmer's dataset:
```
cat palmer.arthCOI.fasta | grep '^>' | cut -d ',' -f2- | awk -F ',' '{print $2,$3,$4,$5,$6}' OFS=',' | gzip --best > palmer.taxa.txt.gz
```

For the Porter dataset we used the dereplicated data only:
```
zcat porter.derep.fasta.gz | grep "^>" | sed 's/>//' > seqids.tmp
awk 'FNR==NR {hash[$1]; next} $1 in hash;' seqids.tmp porter_outmap.txt > match.tmp
rm seqids.tmp
cat match.tmp | cut -f 1 > left.tmp
cat match.tmp | cut -f 2- > right.tmp
rm match.tmp
paste -d ';' left.tmp right.tmp | awk -F ';' '{print $1, $3, $4, $5, $6, $7}' OFS=';' | gzip --best > porter.taxa.txt.gz
rm left.tmp right.tmp
```

For our datasets:
```
zcat boldCOI.derep.txt.gz | cut -d ';' -f 3- | gzip --best > derep.taxa.txt.gz
zcat boldCustom.allArth.meta.txt.gz | cut -d ';' -f 8-13 | gzip --best > raw.taxa.txt.gz
```
> Note that we modified our dereplicated and raw datasets, `boldCOI.derep.txt` and `boldCustom.allArth.meta.txt.gz`, to reduce disk space by pulling just the taxa strings from Class through Species:

The `*lengths.txt` and `*taxa.txt.gz` files were used to generate the tables with the R script `database_composition.R`.

## Impact of clustering on taxonomic completeness
In addition to the comparisons evaluating taxonomic completeness based on the database being evaluated, we analyzed how clustering a dataset can effect the remaining taxa. We clustered the `boldCOI.derep.fasta` file at three values: 99%, 97%, and 95% identity using Vsearch as described in the `database_construction.md` document. Clustering reduced the total number of unique sequences as follows:
- `boldCOI.derep.fasta` contained **1,841,946** unique sequences
- `boldCOI.clust99.fasta` contained **407,356** unique sequences of which **216,866** were singletons
- `boldCOI.clust97.fasta` contained **265,885** unique sequences of which **119,070** were singletons
- `boldCOI.clust95.fasta` contained **215,055** unique sequences of which **87,930** were singletons

These data were presented in [Figure5](https://github.com/devonorourke/tidybug/blob/master/figures/figure5_clustSelectOrders.png) of the mansucript. The [database_clustering.R](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/figure5_database_clustering.R) R script processed the `clust9*.names.txt.gz` files to understand how taxonomic completeness varies as a function of clustering percent identity.

## Impact of geographic specificity
We restricted these analyses to our own datasets - the raw and dereplicated records. The `boldCustom.allArth.meta.txt` served as the raw input, though we filtered the necessary fields to reduce files size:
```
zcat boldCustom.allArth.meta.txt.gz | cut -d ';' -f 1,5-6,8-12 | gzip --best > reduced.allArth.meta.txt.gz
```
The Sequence ID's present from the `boldCOI.derep.txt` file served as input for the dereplicated dataset.

To understand the impact on geographic specificity on taxonomic composition of a database, we filtered our BOLD records by requiring that the `country` column in the BOLD dataset contain either `United States` or `Canada` records. We then investigated how this effects the taxonomic composition at the Class and Order levels, as well as the proportional representation of sequence records from the top 10 sources listed in the `institution_storing` field of the BOLD dataset.

The [database_geography.R](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/database_geography.R) script was used to generate the associated plots (not shown in manuscript, but available in [this directory of the repo containing unpublished images](https://github.com/devonorourke/tidybug/blob/master/figures).  

We find that the vast majority of data comes from one place: Centre for Biodiversity Genomics (CBG). While GenBank (GBK) records also contribute a significant fraction, they are more significant when data is filtered to US or Canada country codes. For the raw BOLD records without any geographic filtering, CBG possesses 62.7% of all sequences, while GBK contains 8.9%, and all other institutions contain about 2.8%. For dereplicated data without any geograhpic filtering the results are similar: 62.7% from CBG, 10.2% from GBK, and 2.7% from others. Yet once a record is filtered as being identified from the United States or Canada, the proportion of records from CBG increases substantially: 89.5% from CBG, less than 1% from GBK, and 9.6% from other institutions. A similar trend is observed for the dereplicated data: 89.9% CBG, less than 1% for GBK, and about 9.1% for other institutions. See [this image](https://github.com/devonorourke/tidybug/blob/master/figures/unused_images/db_6_geography.png) for a plot of the comparisons of data by source and [this image](https://github.com/devonorourke/tidybug/blob/master/figures/unused_images/db_7a-alt_MissingnessByGeography.png) comparing missingness among databases by country.

Nevertheless, filtering does not appear to impact the overall taxonomic composition between Class through Genus levels; that is, when applying the same filtering strategy and limiting only records with "US" or "Canada" information in the BOLD metadata, the most abundant taxa are proportional across Levels. 

Additonal comparisons between pairs of countries were also performed. Interestingly, we find evidence that particular countries - even those that share geographic boarders - often have different distributions of taxa. This may suggest that filtering by a particular country may bias towards particular arthropod orders. See [this image](https://github.com/devonorourke/tidybug/blob/master/figures/unused_images/db_12_databaseCompositionByCountry.png) for details.  

## Next steps
The remaining work related to this manuscript concerns the comparisons of classifiers. See the [classification comparisons](https://github.com/devonorourke/tidybug/blob/master/docs/classification_analyses.md) document for full details. 
