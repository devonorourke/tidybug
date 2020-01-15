# Data collection
Sequences are available for download at the SRA; see http://www.ncbi.nlm.nih.gov/bioproject/518082.

# Raw sequence data processing
Raw sequences were filtered within a Conda environment installed with QIIME2 version 2018.11. See https://docs.qiime2.org/ for full details about installation and complete documentation. Here we show only representative code and descriptions for each of the main processes used in generating the data underlying figures presented in this manuscript; see github.com/devonorourke/tidybug for complete shell scripts used.

## Data import
Within an active Conda environment with QIIME2 installed, data was imported by creating a manifest file as described [in QIIME's import documentation](https://docs.qiime2.org/2019.1/tutorials/importing/?highlight=manifest#fastq-manifest-formats). Because our dataset contained several independent sequence runs worth of fastq files, these were processed in batches by sequencing run. This is also a requirement for one of the downstream denoising applications, Dada2, which builds error models on a per-sample basis.
See `seqfilter.import_example` for template shell script.

> Note in the code below that `$LIB` signifies the path to a directory containing a single sequencing run's worth of fastq files

```
## each $LIB represents a directory containing
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$LIB".manifest.file \
  --output-path "$LIB".demux.qza \
  --input-format PairedEndFastqManifestPhred33
```

The resulting `$LIB.demux.qza` artifact is then used as input for read trimming. This artifact contains unjoined fastq files.

## Primer trimming and visualization
Forward and reverse fastq file pairs are trimmed separately with Cutadapt within the QIIME environment. Note that adapter trimming was performed within the previous [seqfilter.import_example.sh](https://github.com/devonorourke/tidybug/blob/master/scripts/shell_scripts/seqfilter.import_example.sh) shell script.

```
qiime cutadapt trim-paired \
--i-demultiplexed-sequences "$LIB".demux.qza \
--p-cores 24 \
--p-adapter-f GGTCAACAAATCATAAAGATATTGG \
--p-adapter-r GGATTTGGAAATTGATTAGTWCC \
--o-trimmed-sequences "$LIB".trimd.qza
```
The output is a set of trimmed fastq files which remain separate with respect to the forward and reverse fastq pairs. Per-base error profiles are visualized with a QIIME summary tool which demonstrated that the amplicons produced are of sufficiently high quality across the full length of the expected COI read. This information was used to define the truncation length in the subsequent denoising step.
```
qiime demux summarize \
  --i-data "$LIB".trimd.qza \
  --p-n 500000 \
  --o-visualization "$LIB".trimd.qzv
```

## Denoising
Three separate pipelines were used for further filtering sequence data. Each process ultimately generates a table of dereplicated sequences (OTUs or ASVs) containing per-sample counts of unique sequence variants. These three processes were Vsearch, Dada2, or Deblur. Dada2 and Deblur incorporate unique error models to identify and act upon proposed sequence mistakes including chimeric sequences, while Vsearch only corrects for chimeric sequences without any error detection on a per-base level. Each process required a separate set of commands; notably, Dada2 and Denoise accept the trimmed unpaired output from Cutadapt, while Vsearch requires reads to be joined prior to beginning it's process. We changed a few default parameters from the QIIME implementations of Deblur and VSEARCH in attempts to standardize their filtering assumptions:
- all methods removed per-sample features (ASVs) that contained an abundance of sequence counts < 2 (discarded singleton sequences per sample)
- all methods retained per-dataset features (ASVs) _including_ singleton ASVs

### Dada2
Dada 2 implementation below. See [seqfilter.dada2_example.sh](https://github.com/devonorourke/tidybug/blob/master/scripts/shell_scripts/seqfilter.dada2_example.sh) for example shell script used.
```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$LIB".trimd.qza \
  --p-trunc-len-f 175 \
  --p-trunc-len-r 175 \
  --p-n-threads 24 \
  --o-table "$LIB".table.qza \
  --o-representative-sequences "$LIB".repSeqs.qza \
  --o-denoising-stats "$LIB".denoisingStats.qza
```

DADA2 automatically discards singelton seqences, but retains doubletons, tripletons, etc. It also retains singleton ASVs by default in QIIME. Note that chimera filtering completed in DADA2 is done on a per-sample basis rather than by pooling the entire sample; this follows the recommendation of the developer (though the per-pool chimera detection option is available in QIIME).

### Deblur
Unliked Dada2, Deblur uses a reference database of known sequences to model the error of the sequence data. The database used for denoising here is the same as described in the `database_workflow.md` document - a set of select arthropod sequences derived from the Barcode of Life Database (BOLD). See [seqfilter.deblur_example.sh](https://github.com/devonorourke/tidybug/blob/master/scripts/shell_scripts/seqfilter.deblur_example.sh) for example shell script used.  
> The `$REF` environmental variable specifies the full path to the QIIME artifact that represents the curated arthropod BOLD database

```
qiime deblur denoise-other \
--i-demultiplexed-seqs ../"$LIB".qfiltd.seqs.qza \
--i-reference-seqs "$REF" \
--p-trim-length 181 \
--p-sample-stats \
--p-min-reads 2 \
--p-min-size 1 \
--p-jobs-to-start 24 \
--o-table "$LIB".dblr.table.qza \
--o-representative-sequences "$LIB".dblr.repseqs.qza \
--o-stats "$LIB".dblr.stats.qza
```

Unlike DADA2, Deblur by default discards singleton ASVs, and discards per-sample ASVs with less than 10 reads. We modified these default parameters to match DADA2. The chimera filtering performed in this analysis uses a VSEARCH implementation of Uchime-denovo and is thus independent of a reference library. This requires that chimera filtering is done on a pool of sequences, rather than on a per-sample level.

### Vsearch
The parameters chosen in this implementation mostly follow those described by Vsearch authors in [their Wiki documentation](https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline). The exception is that their vignette removes singleton features, while we are retaining them; however, we are going to add an additional step in which per-sample ASVs with < 2 reads are discarded to match the DADA2 and Deblur filtering settings. An example Slurm script is available: see [seqfilter.vsearch_example.sh](https://github.com/devonorourke/tidybug/blob/master/scripts/shell_scripts/seqfilter.vsearch_example.sh).  

Within the QIIME installation, Vsearch filtering first required the read pairs to be joined first; we then applied a basic quality filter to the sequences:
```
## merge paired end data
## '--p-allowmergestagger' necessary because of the large overhang with our 300bp seqs and 180bp amplicon
qiime vsearch join-pairs \
--p-allowmergestagger \
--i-demultiplexed-seqs "$LIB".trimd.qza \
--o-joined-sequences "$LIB".joind.seqs.qza

## quality filter seqs: used for both deblur and vsearch piplelines
qiime quality-filter q-score-joined \
 --i-demux "$LIB".joind.seqs.qza \
 --o-filtered-sequences "$LIB".qfiltd.seqs.qza \
 --o-filter-stats "$LIB".qfiltd-stats.qza
 ```

 We then dereplicate the merged reads, cluster at 98% before chimera detection, perform the de novo chimera filtering, discard the chimeras, then cluster once more at 97% identity to produce the OTU table and representative sequence variants.

 ```
 ## dereplicate the merged dataset
qiime vsearch dereplicate-sequences \
--i-sequences ../"$LIB".qfiltd.seqs.qza \
--o-dereplicated-table "$LIB".vsrch.derep.table.qza \
--o-dereplicated-sequences "$LIB".vsrch.derep.seqs.qza

## pre cluster at 98% before chimera detection
qiime vsearch cluster-features-de-novo \
--i-sequences "$LIB".vsrch.derep.seqs.qza \
--i-table "$LIB".vsrch.derep.table.qza \
--p-perc-identity 0.98 \
--p-threads 24 \
--o-clustered-table "$LIB".vsrch.preChimera.table.qza \
--o-clustered-sequences "$LIB".vsrch.preChimera.seqs.qza

## perform de novo chimera filtering
qiime vsearch uchime-denovo \
--i-sequences "$LIB".vsrch.preChimera.seqs.qza \
--i-table "$LIB".vsrch.preChimera.table.qza \
--o-chimeras "$LIB".vsrch.chimera.seqs.qza \
--o-nonchimeras "$LIB".vsrch.nonchimera.seqs.qza \
--o-stats "$LIB".chimerafilt-stats.qza

## generate stats file from chimera filtering
qiime metadata tabulate \
  --m-input-file "$LIB".chimerafilt-stats.qza \
  --o-visualization "$LIB".chimerafilt-stats.qzv

## exclude chimeras and borderline chimeras:
qiime feature-table filter-features \
  --i-table "$LIB".vsrch.preChimera.table.qza \
  --m-metadata-file "$LIB".vsrch.nonchimera.seqs.qza \
  --o-filtered-table "$LIB".vsrch.postChimera.table.qza

## cluster at 97% identity
qiime vsearch cluster-features-de-novo \
--i-sequences "$LIB".vsrch.nonchimera.seqs.qza \
--i-table "$LIB".vsrch.postChimera.table.qza \
--p-perc-identity 0.97 \
--p-threads 24 \
--o-clustered-table "$LIB".vsrch.postChimera.p97clust.table.qza \
--o-clustered-sequences "$LIB".vsrch.postChimera.p97clust.seqs.qza

qiime tools export --input-path "$LIB".vsrch.postChimera.p97clust.seqs.qza --output-path "$LIB".p97clust
qiime tools export --input-path "$LIB".chimerafilt-stats.qzv --output-path "$LIB".chimerafilt-stats.txt
 ```

However, this creates a problem when comparing across filtering pipelines: vsearch by default in QIIME's implementation creates a sha1-hashed label, whereas dada2 and deblur use an md5 hashing algorithm to label sequences. As a result, you can't compare identical sequences. To fix this, we reformatted the resulting `"$LIB".vsrch.postChimera.p97clust*.qza` artifacts by (1) appending the sequences with the md5 hash using a native Vsearch function; generating a list of the equivalent md5 and sha1 hashIDs for each sequence; using these hashID pairs to relabel the original .qza feature table. Rather than doing this four times (one for each library), we're going to first combine the four libraries into a single dataset and execute the script just once.

## Combining datasets
Following these three basic filtering pipelines, the per-`$LIB` tables and representative sequences were joined into individual `.qza` artifacts. See the `seqfilter.combineNfilter.sh` shell script for full details. As an example, this is how the Dada2-filtered sequencing batches of tables and sequences were combined (the same would apply for deblur and vsearch):

```
## merge the features for all dada2-produced tables and sequences
# tables
qiime feature-table merge \
  --i-tables /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p41/qiime/reads/dada2/p41.dada2.table.qza \
  --i-tables /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p42/qiime/reads/dada2/p42.dada2.table.qza \
  --i-tables /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p71/qiime/reads/dada2/p71.dada2.table.qza \
  --i-tables /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p72/qiime/reads/dada2/p72.dada2.table.qza \
  --o-merged-table dada2.all.raw.table.qza

# sequences
qiime feature-table merge-seqs \
  --i-data /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p41/qiime/reads/dada2/p41.dada2.repSeqs.qza \
  --i-data /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p42/qiime/reads/dada2/p42.dada2.repSeqs.qza \
  --i-data /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p71/qiime/reads/dada2/p71.dada2.repSeqs.qza \
  --i-data /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p72/qiime/reads/dada2/p72.dada2.repSeqs.qza \
  --o-merged-data dada2.all.raw.seqs.qza
```

## Amending the Vsearch data
The combined `vsearch.all.raw.*.qza` files needed to be relabled from their sha1 format to the md5 hash format, as described above. We did this by creating a test file `vsearch.all.headers.txt` that contained the sha1 and md5 headers in two columns:
```
## reformat the Vsearch fasta file with md5 hash labels, while retaining the original sha1 hash ID
qiime tools export --input-path vsearch.all.raw.seqs.qza --output-path vtmpseqs
vsearch -derep_fulllength ./vtmpseqs/dna-sequences.fasta -relabel_md5 -relabel_keep -xsize -output vsearch.all.md5.seqs.fasta -fasta_width 0
## make a list of these headers
echo "md5 sha1" > vsearch.all.headers.txt
grep "^>" vsearch.all.md5.seqs.fasta | sed 's/>//' >> vsearch.all.headers.txt
rm -r vtmpseqs
```

We then run the following R script which converts the .qza artifact with the sha1 hashIDs into the md5 version, using the "$LIB".headers.txt information provided above:
```
# required: library(devtools)
# runonce: install_github("jbisanz/qiime2R", force=TRUE)
library(tidyverse)
library(reshape2)
library(qiime2R)
library(phyloseq)

featuretable <- read_qza("vsearch.all.raw.table.qza")   ## convert .qza to matrix, then convert wide-format matrix to long-format data.frame object
mat.tmp <- featuretable$data
df.tmp <- as.data.frame(mat.tmp)
df.tmp$OTUid <- rownames(df.tmp)
rownames(df.tmp) <- NULL
vsearch.tmp <- melt(df.tmp, id = "OTUid") %>% filter(value != 0) %>% mutate(Method = "vsearch")
colnames(vsearch.tmp) <- c("sha1", "SeqID", "Reads", "Method")

hashtable <- read.table(gzfile("vsearch.all.headers.txt"), header = TRUE)   ## import hashID lists and resolve sha1 vs MD5-hash'd strings in the headers
colnames(hashtable) <- c("md5", "sha1")
df <- merge(vsearch.tmp, hashtable, all.x = TRUE)       ## merge with dataframe 'df' object
df <- df[c(2,3,5)]
df$md5 <- as.character(df$md5)
df.mat <- dcast(df, md5 ~ SeqID, value.var = "Reads", fill = 0)   ## reformat to matrix
colnames(df.mat)[1] <- "#OTU ID"
write.table(df.mat, "vsearch.all.raw.md5.table.tsv", quote = FALSE, row.names = FALSE)    ## write file to disk
```

We then convert that text table into a QIIME .qza file along with the reformatted .fasta file.
```
## convert fasta back into proper .qza object
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path vsearch.all.md5.seqs.fasta \
  --output-path vsearch.all.raw.seqs.qza

## convert .tsv table into .biom, then .biom back into proper .qza object
biom convert -i vsearch.all.raw.md5.table.tsv -o vsearch.all.raw.md5.table.biom --table-type="OTU table" --to-hdf5

qiime tools import \
  --input-path vsearch.all.raw.md5.table.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path vsearch.all.raw.table.qza
```
> Note that the `vsearch.all.raw.*.qza` files were overwriting the original files used as input into the whole relabeleing process

We removed host sequences from this table first before removing the remaining singleton sequences.

## Filtering potential host sequences
It's possible that following quality filtering the remaining data includes host (bat) COI sequences. Failure to remove these sequences can result in misidentifying the bat COI sequence as some arthropod sequence (because the database that sequences are assigned taxonomy consists exclusively of arthropod information). Because we captured our data within New Hampshire there are a finite list of possible bat hosts. Representative sequences for each bat species were obtained from two Genbank projects:

1. PopSet 726974368: MYSO, MYSE, and MYLE
2. PopSet 301344216: MYAU, MYLU, PESU, EPFU, NYHU, LABO, COTO, LACI, LANO

The resulting fasta files were downloaded for each species. Note that the fasta files needed to be reformatted from original format to work within QIIME, which expects a USEARCH-style file pair. This results in two files created from a single fasta file:
1. fasta file containg header lacking taxonomic info plus sequence data
2. text file containing fasta header and taxonomy information

The resulting datasets were refomratted and are available as `nh.bat.hosts.fasta` and `nh.bat.hosts.txt`, respectively. These are available within the /data/databases section of the Github repo. The `nh.bat.hosts.*` files were imported into qiime for use in filtering out host reads:
```
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path nh.bat.hosts.fasta \
  --output-path NHbat_seqs.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path nh.bat.hosts.txt \
  --output-path NHbat_tax.qza
```

These file pairs were used in the `seqfilter.combineNfilter.sh` script to remove bat COI sequences. In brief, we first determine which sequences align to our host reference (bat) representaive sequences:
> Note here that `$METHOD` indicates the path to the directory containing each filtered table and sequence set (from Dada2, deblur, and vsearch), while `$REF` is the path to the **bat** reference database (not the arthropod BOLD database)

```
  qiime quality-control exclude-seqs \
    --i-query-sequences "$METHOD".all.raw.seqs.qza \
    --i-reference-sequences "$REF" \
    --p-perc-identity 0.95 \
    --o-sequence-hits "$METHOD".hostseqs.qza \
    --o-sequence-misses "$METHOD".arthseqs.qza \
    --p-method vsearch;

  qiime tools export --input-path "$METHOD".hostseqs.qza --output-path "$METHOD"_hostseqs;

  grep "^>" "$METHOD"_hostseqs/dna-sequences.fasta | sed 's/>//' | sed '1 i\#OTU ID' > "$METHOD".droplist;
```

The resulting `"$METHOD".hostseqs.fasta` file was searched with NCBI blast against nr database to determine if putative host hits were mismatches. Nearly all matches regardless of filtering method returned the same host for putative host-associated sequence variants: the Little Brown bat, _Myotis lucifugus_. There were two exceptions: one sequence variant was matched to the Eastern Small-footed bat, _Myotis leibii_ in each of the three filtering methods, while the Vsearch and DADA2 filtering strategies matched a sequence variant to the Evening bat, _Nycticeius humeralis_ (this variant was not present in the DADA2-filtered output).

The `*.droplist` file is then used to exlcude the ASVs we suspect are derived from bat hosts; any further sequence variants are assumed to be derived from arthropod (or otherwise will remain unclassified):
```
qiime feature-table filter-features \
  --i-table "$METHOD".all.raw.table.qza \
  --m-metadata-file "$METHOD".droplist \
  --o-filtered-table "$METHOD".arthtable.qza \
  --p-exclude-ids
```

For VSEARCH specifically, the `vsearch.arthtable.qza` file is filtered to removed singleton sequences per sample:
```
qiime feature-table filter-features \
--i-table vsearch.arthtable.qza \
--p-min-frequency 2 \
--o-filtered-table vsearch.arthtable_nosingles.qza

mv vsearch.arthtable_nosingles.qza vsearch.arthtable.qza
```

We also then update the sequences based on this reduced table (which contains fewer ASVs):
```
qiime feature-table filter-seqs \
--i-data vsearch.arthseqs.qza \
--i-table vsearch.arthtable.qza \
--o-filtered-data vsearch.arthseqs-filtd.qza
mv vsearch.arthseqs-filtd.qza vsearch.arthseqs.qza
```

The `*.arthtable.qza` files produced represent our "basic" host-filtered datasets, one per filtering program (Dada2, Deblur, or Vsearch). These artifacts serve as the input for potentially additional sample and sequence variant-based filtering described next. Files were uploaded to the Github repo (/tidybug/data/qiime/tables).

## Filtering datasets
The "standard" and "extra" tables are created from the `*arthseqs.qza` and `*arthtable.qza` objects (from dada2, deblur, and vsearch-filtered inputs). We opted to use a combination of QIIME-supported functions as well as a custom R script `sequence_filtering.R` to generate the "standard" and "extra" datasets because of the need for generating the figures used in this manuscript and to allow for additional flexibility in statistical tests not supported with QIIME2 version 2018.11 (note that additional packages like the Adonis plugin for PERMANOVA testing of group differences were added to QIIME2 v.2019.1, for example). We provide context for how the "standard" and "extra" filters are designed below, and have commented within the R script for additional clarification.

### Standard filtering approach
The "standard" filter will apply two criteria for inclusion: (1) samples must have at least 5000 sequence reads (from the `*arthtable.qza` artifact), and (2) keeps only those sequence variants observed in at least 2 samples across the entire dataset. The rationale behind dropping low-abundance samples is that sequence errors are more likely in samples with low abundances of reads. We perform the read-minimum filtering on samples first, then apply the sequence variant-minimum filter. As a result, a sample following a the "standard" filter may have less than 5000 reads if one of the OTUs happened to be dropped from that sample. This dataset was created by importing the `*arthtable.qza` file for each pipeline (dada2, deblur, or vsearch) into an R environment and applying the `sequence_filtering.R` script.

### Extra filtering approach
The "extra" filter uses the same filtering criteria as the "standard" filter, but further requires that a fixed-count subtraction is applied to all reads per sequence variant per sample. This value is determined by first identifying which ASVs in a given mock sample are unexpected: those samples that have less than 97% identity to our known mock sequences. Once we have a list of the unexpected, or "miss" sequences, we identify the sequence variant among those "miss" variatns with the highest read count. That maximum value serves as the integer with which all sequence variants are reduced by.

A reference QIIME artifact consisting of the known (expected) mock sequences was generated from the fasta file of known community members and sequences (`CFMR_insect_mock4.fasta`) and uploaded into QIIME format:

```
qiime tools import \
  --input-path CFMR_insect_mock4.fasta \
  --output-path mock.ExpectedSeqs.qza \
  --type 'FeatureData[Sequence]'
```  
> See a similar `CMFR_insect_mock4_wtax.fasta` file for expected taxonomic information.

To create a feature table consisting of just the four mock samples, we first create a small metadata file listing the 4 mock samples of interest:
```
echo SampleID > mocklist.txt
echo mockIM4p4L1 >> mocklist.txt
echo mockIM4p4L2 >> mocklist.txt
echo mockIM4p7L1 >> mocklist.txt
echo mockIM4p7L2 >> mocklist.txt
```

Because all filtered samples contain the same mock sample names, we can apply that file to filter each feature table to include just the mock samples (making a mock-only feature table). We apply this to each of the 3 `*arthtable.qza` files (example shown here for dada2):
> Note the `$MOCKLIST` variable shown below describes the full path to the `mocklist.txt` file

```
qiime feature-table filter-samples \
  --i-table dada2.arthtable.qza \
  --m-metadata-file "$MOCKLIST" \
  --o-filtered-table dada2.mock.table.qza
```

We then use that table as the input to filter out just the sequences we need - those sequences identified in the mock samples:
```
qiime feature-table filter-seqs \
  --i-data dada2.arthseqs.qza \
  --i-table dada2.mock.table.qza \
  --o-filtered-data dada2.mock.seqs.qza
```

This mock sample file is then used to identify the "exact", "partial", and "miss" sequence variants among each of the four mock samples. We used the QIIME quality-control plugin to identify exact and partial matches via vsearch's global alignment algorithms (example shown here for dada2):
> The `$MOCKREF` variable shown below refers to the full path to the `mock.ExpectedSeqs.qza` file created previously.

```
## use 100% identity for the exact matches; default coverage length is 97%
qiime quality-control exclude-seqs --p-method vsearch \
--i-query-sequences dada2.mock.seqs.qza \
--i-reference-sequences $MOCKREF \
--p-perc-identity 1.0 \
--o-sequence-hits dada2.mock.exactHits.qza \
--o-sequence-misses dada2.mock.exactMisses.qza

## use 97% identity for partial matches
qiime quality-control exclude-seqs --p-method vsearch \
--i-query-sequences dada2.mock.exactMisses.qza \
--i-reference-sequences $MOCKREF \
--p-perc-identity 0.97 \
--o-sequence-hits dada2.mock.partialHits.qza \
--o-sequence-misses dada2.mock.partialMisses.qza

rm dada2.mock.exactMisses.qza
```

We then export the each of the `*.mock.qza` artifacts and save the HashIDs for each of the mock samples.
```
qiime tools export --input-path dada2.mock.exactHits.qza --output-path exact_seqs
qiime tools export --input-path dada2.mock.partialHits.qza --output-path partial_seqs
qiime tools export --input-path dada2.mock.partialMisses.qza --output-path miss_seqs

grep "^>" ./exact_seqs/dna-sequences.fasta | sed 's/>//' > dada2.exactseqs.txt
grep "^>" ./partial_seqs/dna-sequences.fasta | sed 's/>//' > dada2.partialseqs.txt
grep "^>" ./miss_seqs/dna-sequences.fasta | sed 's/>//' > dada2.missseqs.txt

```

These `*.txt` files are imported into the custom R script `sequence_filtering.R` to determine what the maximum value is per Library, per pipeline. Each Library of full guano data (or mock data) is then filtered according to that value.

# Additional notes
The output of this script and the additional `sequence_filtering.R` script generates all of the needed data structures for the section of the manuscript concerning filtering pipelines and parameters. Additional R scripts were applied to generate the datasets and figures concerning alpha and beta diversity metrics. See the `database_filtering.md` document for details of bioinformatic processes used in the database-related section of the manscript. Additional R scripts are generated for each figure as needed; each figure has it's own R script for clarity.
