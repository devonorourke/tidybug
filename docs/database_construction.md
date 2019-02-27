# Database curation
I've found that there are a handful of filtering steps that are useful when creating an arthropod COI database from BOLD records. To accomplish this, we used a few virtual environments for the work; `arrR` was used for data acquisition from BOLD, as well as filtering the resulting BOLD records,  `dev_qiime1` was used to apply further filtering parameters using some QIIME1 scripts and related Python scripts, and `qiime2-2018.11` was used to import the final filtered COI dataset. See the `sequence_filtering.md` document for QIIME2 information; the additional virtual environments were created as follows:

```
conda create -n arrR python=3.6
conda install -c r r-base
conda install -c bioconda seqkit
conda install -c bioconda csvtk
source activate arrR
```

However, we also required a few QIIME1 scripts as part of the dereplication process:
```
conda create -n dev_qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda
```
> - I'd recommend testing proper installation, as I've encountered times where certain python packages need to be reinstalled
> - try launching the **dev_qiime1** environment, then type: `print_qiime_config.py -t`
> - if you get an error, uninstall and reinstall the offending python package (for example, I've had to do this for scipy and biom-format before)


## Obtaining data
We followed the vignette provided by the [bold R package](https://github.com/ropensci/bold) to download data from the [Barcode of Life Database](http://v4.boldsystems.org/). See the R script `bold_datapull.R` for full details on how the raw BOLD data was obtained and filtered. In brief, we selected all COI records from BOLD by pulling all records matching the Arthropod taxon; data was filtered to require the `markercode` column matching "COI-5P", and the `phylum_name` column matching "Arthropoda". The output of this script contains a single `.csv` file that contains sequence information, taxonomic information, and associated metadata including the `sequenceID`, `processid`, `bin_uri`, `genbank_accession`, `country`, `institution_storing` records.

Because we're pulling nearly two million records from BOLD, we'll query the servers iteratively by generating a list of groups to pull from. I found that you can avoid any server complications by keeping the list size under about 2 million records, so I divied up all Arthropods into the following groups:
1. all non-Insect arthropods
2. Dipterans
3. Coleopterans
4. Hymenopterans
5. Lepidopterans
6. all other Insects

See the R script `bold_datapull.R` for implementation on how data was queried; notably, I did not follow the vignette exactly. One particular sticking point was how the names are derived when using the `taxize` library - I found that it can erronously _miss_ some taxa because of different names being applied to a group. One such instance I identified occurred between the NCBI name of `Psocoptera` while BOLD uses the term `Psocodea`; the distinction is that NCBI things `Psocodea` is a superOrder instead of just an Order, and BOLD doesn't do the super/supra level distinctions. Because of this discrepency, you'd completely miss all `Psocodea` records if you followed the vignette exactly.

The retrieved records were split into two files: the first is a file containing sequence information with a particular format of taxonomic information, while the second contains metadata with taxonomy information in a different format. The first file taxonomic info is formatted in a way to produce the expected delimiters used with QIIME2 and Vsearch tools want (in case you want to dereplicate, for example), while the metadata file contained the various attributes for each Sequence ID, as well as the Phylum through Species taxa names in their original format (ie. not concatenated together). The latter allows us to quickly parse and tally the data as we impose certain filters.

## Filtering BOLD data
We filtered our raw BOLD sequence records by removing short sequences, sequences with gap characters, sequences with DNA characters that were not part of the standard IUPAC code, and sequences with extremely sparse taxonomic information.

### Filter1 - removing gaps and filtering for sequence length
The initial data set contains **3,175,054** records that passed our initial filtering requiring that a `markercode` column in the BOLD metadata be assgined with `COI-5P`. However, we identified a few hundred records in the initial BOLD data that contain gaps. We removed any gap character then selected only sequences longer than 100 bp to produce a single fasta file formatted for further processing in Vsearch. There were just **552** sequences that were discarded, so the vast  majority of the records (**3,174,502**) were retained:

```
grep -v 'sequenceID,taxon,nucleotides' boldCustom.allArth.seqNtaxa.csv | \
awk -F ',' '{OFS="\t"};{print $1";tax="$2, $3}' | \
sed 's/^/>/g' | tr '\t' '\n' | \
seqkit seq - --min-len 100 --remove-gaps -w 0 > boldCustom.allArth.filt1.fasta
```

### Filter2 - identifying any non-ATCG sequence characters
The BOLD sequences do not exlcusively contain non-ambiguous DNA characters. To see the distribution of which characters are present:
```
seqkit fx2tab -n -i -a boldCustom.allArth.filt1.fasta | \
csvtk -H -t grep -f 4 -r -i -p "[^ACGT]" | \
awk '{print $2}' | sort | uniq -c | sort -k1nr
```

About 90% of all records contain just ATCG characters, but the remaining contain either `N` or other ambiguous characters. Dig deeply into that list and you'll find that a single sequence contains a character that is not part of the standard DNA alphabet: `I` (it's the one with `ACGIT` characters). We need to remove that sequence entirely. To identify which sequence:
```
seqkit fx2tab -n -i -a boldCustom.allArth.filt1.fasta | \
csvtk -H -t grep -f 4 -r -i -p "I" | cut -f 1 > idsWithI.txt
```
> The result is record: `8982958;tax=k__Animalia;p__Arthropoda;c__Insecta;o__Diptera;f__Ephydridae;g__Hydrellia;s__Hydrellia			ACGIT`. This Order, Family, and Genus is well represented in the database, so dropping this one sequence is not going to affect our classification

We can discard that single sequence to create our fasta file filtered by length, gaps, and DNA alphabet:
```
seqkit grep --pattern-file idsWithI.txt --invert-match boldCustom.allArth.filt1.fasta > boldCustom.allArth.filt2.fasta
```

### Filter3 - removing sequences with taxonomic incompleteness
Prior to dereplicating data, we considered whether or not to impose an additional filtering parameter relating to taxonomic completeness - should we discard any sequences that contain no information at a particular taxonomic level. There are **3,175,053** total sequences in `boldCustom.allArth.filt2.fasta`, and while all of these sequences contain Phylum and Class information, **1,856** of these are missing Order information, and **122,709** were missing Family level information. The majority of sequences missing Order level information were associated with the _Arachnida_ Class (1,346), with _Collembola_ (281) and _Malacostraca_ (102) also missing several sequences.  Among the missing Family level data, the _Insecta_ class was missing the most sequences (62854), with the _Arachnida_ (29688), _Collembola_ (16749) and _Malacostraca_ (9622) again missing the most information among non-Insect arthropods. The Orders represented within the _Insecta_ Class with missing Family information were unequally distributed, with _Hempitera_ (35679) missing more than any other Order, though notably each Order has a unique number of total records.

We opted to create an filtered dataset to be used in dereplication in which all sequences that lacked Order or Family information were discarded. This ensured that dereplicating data with our lowest common ancestor (LCA) approach would not suffer from information loss in the situation where two identical sequences containing different taxonomic records - one without, say, Order level information, and another with complete information - was reduced solely because of different levels of information completeness (as opposed to distinct information at equivalent levels). Filtering follows a similar approach in identifying the sequence headers without Family info, then removing them from the fasta file:

```
seqkit fx2tab -n -i -a boldCustom.allArth.filt2.fasta | csvtk -H -t grep -f 1 -r -i -p "f__;" | cut -f 1 > idsNoFam.txt

seqkit -w 0 grep --pattern-file idsNoFam.txt --invert-match boldCustom.allArth.filt2.fasta > boldCustom.allArth.filt3.fasta
```

The resulting `boldCustom.allArth.filt3.fasta` now contains **3,051,814** Arthropod records with at least Family-level information, the proper DNA alphabet, and sequences longer than 100 bp (about 96% of what we started with).
> A `missingFamilyCounts.txt` documents the number of missing sequence records, grouped by taxonomic Orders. See the tidybug/data/databases

## Dereplicating with LCA
Databases can contain redundant sequences; dereplicating datasets is one solution to this problem, however, the default dereplication tools used in programs like Vsearch will select the first record when multiple sequence matches exist. This can be problematic if the two records contain non-identical taxonomic information; this generally can create one of two problems. First, if the sequences contain non-identical records but equally complete levels of taxonomic information, there are two generally adopted strategies to picking a "winner" - majority, or consensus. A majority approach will take the most abundant of classifications, while a consensus approach invokes a least common ancestor (LCA) algorithm which retains only taxonomic information at the level where the matching sequences are equal. The second problem that can occur is if the two sequences contain unequal levels of information - a majority or consensus approach can again be invoked, but in this case you will always lose information with a consensus approach to whichever record contains the least amount of taxonomic information.
We adapted the [instructions](https://github.com/mikerobeson/Misc_Code/tree/master/SILVA_to_RDP) written by Mike Robeson for formatting a SILVA database, and incorporated the [consensus approach to reclassifying taxa](https://gist.github.com/walterst/9ddb926fece4b7c0e12c) script written by William (Tony) Walters. This resolved the problem of ties whereby identical sequences with distinct taxa are resolved to a common taxonomic level. This also required installing QIIME1 - we did so by creating a new virtual environment; see [here for instructions](http://qiime.org/install/install.html) on installing QIIME1. The actions are to first generate the appropriate file structures to perform the LCA algorithm on all data, then, with those "ties" now resolved, we can dereplicate the data knowing that any potential duplicate sequence selected will have the appropriate taxonomy assigned.

First, create the taxonomy mapping file; it's a two-column record of sequenceID and taxonomy string.
```
## Taxonomy mapping file:
cat boldCustom.allArth.filt3.fasta | grep '^>' | sed 's/^>//' | cut -d ';' -f 1 > tmp.left
cat boldCustom.allArth.filt3.fasta | grep '^>' | sed 's/^>//' | cut -d ';' -f2- | sed 's/tax=k__//' | sed 's/p__//' | sed 's/c__//' | sed 's/o__//' | sed 's/f__//' | sed 's/g__//' | sed 's/s__//' > tmp.right
paste -d '\t' tmp.left tmp.right > tmp_nolabels.taxa
rm tmp.left tmp.right
```

Next, create a reduced fasta file without taxa info in the headers to save disk space when writing the subsequent outfiles:
```
## reduced fasta file:
cat boldCustom.allArth.filt3.fasta | grep '^>' | cut -d ';' -f 1 > tmp.left
cat boldCustom.allArth.filt3.fasta | grep -v '^>' > tmp.seqs
paste -d '\t' tmp.left tmp.seqs | tr '\t' '\n' > tmp_nolabels.fasta
rm tmp.left tmp.seqs
```

We'll use that `tmp_nolabels.fasta` file as input for the [pick_otus.py](http://qiime.org/scripts/pick_otus.html) script. We switch the conda environments from the ``
```
conda activate dev_qiime1
pick_otus.py -i tmp_nolabels.fasta -o filt3_otus --similarity 1.0 --threads 24
```
This generates a `tmp_nolabels_otus.txt` text file within a newly created `filt3_otus` directory that is used in the next [create_consensus_taxonomy.py](https://gist.github.com/walterst/bd69a19e75748f79efeb) script to generate the mapping file that will apply the LCA algorithm.

We next apply three inputs (`tmp_nolabels.taxa`, `tmp_nolabels.fasta`, and `tmp_nolabels_otus.txt`) to generate a consensus mapping file output (`outmap.txt`) for our data:
```
python create_consensus_taxonomy.py tmp_nolabels.taxa tmp_nolabels.fasta ./pikd_otus/tmp_nolabels_otus.txt outmap.txt
```

The `outmap.txt`file is a 2-column file with the sequenceID and taxonomy string:
```
5333265	Animalia;Arthropoda;Insecta;Diptera;Cecidomyiidae;;
5333264	Animalia;Arthropoda;Insecta;Hemiptera;Aphididae;Aphis;Ambiguous_taxa
5333267	Animalia;Arthropoda;Insecta;Thysanoptera;Thripidae;;
```

Recall that the `tmp_nolabels.taxa` taxonomy mapping file (created from our filtered `tmp_nolabels.fasta` file) contained **3,051,814** records - this includes potentially redundant sequences. The `outmap.txt` file contains **3,051,814** records also, but we have applied the LCA algorithm to these records. All that remains to do is dereplicate these records; because it doesn't matter which record is chosen among redundant sequences now (they all have similar taxa records) we can dereplicate the `tmp_nolabels.fasta` directly:
```
vsearch \
--derep_fulllength tmp_nolabels.fasta \
--output boldCOI.derep.fasta \
--relabel_keep --threads 4 --fasta_width 0 --notrunclabels
```

The dereplicated data contain **1,841,946** unique records (about 60% of the original data); the Vsearch output provided the following brief summary:
```
Dereplicating file tmp_nolabels.fasta 100%  
1873901062 nt in 3051814 seqs, min 100, max 1982, avg 614
Sorting 100%
1841946 unique sequences, avg cluster 1.7, median 1, max 3547
Writing output file 100%  
```

We'll use this dataset for one final purpose prior to converting into the appropriate QIIME format- clustering this dereplicated set.

## Clustering data
I wanted to understand how clustering databases prior to classification can reduce taxonomic information in a dataset. I suspect the motivation behind this approach is partly due to improving computational speed (fewer sequences means faster searching) while reducing computing hardware (smaller datasets require less memory to compute and less disk space to write to), as well as a biological inclination that clustering sequences with similar sequence identities produces fewer unique sequence variants the resulting dataset can contain. The former argument is certainly important if you want to do this work on a laptop, but not a necessity with the nearly ubiquitous and cheap access to cloud computing resources. The latter argument is just unnecessary - if computing resources are no issue, we shouldn't cluster before hand - we should cluster after the fact if so desired because you can't do the reverse if your database is clustered to begin with.

We'll cluster at 99%, 97%, and 95%; this command produces three output fasta files, one per cluster %id:
```
declare -a arr=("99" "97" "95")
for i in "${arr[@]}"; do
vsearch --cluster_size boldCOI.derep.fasta --threads 24 --strand both --fasta_width 0 --notrunclabels --relabel_keep \
--id 0."$i" --centroids boldCOI.clust"$i".fasta
done
```

These clustered `.fasta` records serve as inputs for subsequent analyses described in the `database_analyses.md` document.

## Formatting data for QIIME import
Two files are required for importing into QIIME2 to perform classification: (1) a taxonomy file, and (2) a fasta file. The taxonomy file uses the same 2 column format as the `outmap.txt` file, except we're only going to retain the records present in the dereplicated dataset. The `boldCOI.derep.fasta` file contains the correct format for import.

We make a temporary list of all the sequenceID values from the headers of the `boldCOI.derep.fasta` file, then use that as a list to query the matches in the `outmap.txt` file to generate the taxonomy file we want:
```
grep '^>' boldCOI.derep.fasta | sed 's/>//' > derep.seqid.tmp
awk 'FNR==NR {hash[$1]; next} $1 in hash' derep.seqid.tmp outmap.txt > boldCOI.derep.txt
rm derep.seqid.tmp
```

The `boldCOI.derep.fasta` and `boldCOI.derep.txt` files are imported into QIIME2:

```
## fasta file import
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path boldCOI.derep.fasta \
  --output-path boldCOI.derep.seqs.qza

## taxonomy file import
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path boldCOI.derep.txt \
  --output-path boldCOI.derep.tax.qza
```
