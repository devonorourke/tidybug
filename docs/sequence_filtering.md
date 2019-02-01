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
Forward and reverse fastq file pairs are trimmed separately with Cutadapt within the QIIME environment. Note that adapter trimming was performed within the previous `seqfilter.import_example.sh` shell script.

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
Three separate pipelines were used for further filtering sequence data. Each process ultimately generates a table of dereplicated sequences (OTUs or ASVs) containing per-sample counts of unique sequence variants. These three processes were Vsearch, Dada2, or Deblur. Dada2 and Deblur incorporate unique error models to identify and act upon proposed sequence mistakes including chimeric sequences, while Vsearch only corrects for chimeric sequences without any error detection on a per-base level. Each process required a separate set of commands; notably, Dada2 and Denoise accept the trimmed unpaired output from Cutadapt, while Vsearch requires reads to be joined prior to beginning it's process.

### Dada2
Dada 2 implementation below. See `seqfilter.dada2_example.sh` for example shell script used.
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

### Deblur
Unliked Dada2, Deblur uses a reference database of known sequences to model the error of the sequence data. The database used for denoising here is the same as described in the `database_workflow.md` document - a set of select arthropod sequences derived from the Barcode of Life Database (BOLD).
> The `$REF` environmental variable specifies the full path to the QIIME artifact that represents the curated arthropod BOLD database

```
qiime deblur denoise-other \
--i-demultiplexed-seqs ../"$LIB".qfiltd.seqs.qza \
--i-reference-seqs "$REF" \
--p-trim-length 181 \
--p-sample-stats \
--p-jobs-to-start 24 \
--o-table "$LIB".dblr.table.qza \
--o-representative-sequences "$LIB".dblr.repseqs.qza \
--o-stats "$LIB".dblr.stats.qza
```

### Vsearch
The parameters chosen in this implementation follow those described by Vsearch authors in [their Wiki documentation](https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline). Within the QIIME installation, Vsearch filtering first required the read pairs to be joined first; we then applied a basic quality filter to the sequences:
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

qiime tools export --input-path "$LIB".vsrch.nonchimera.seqs.qza --output-path "$LIB".p100clust.seqs
qiime tools export --input-path "$LIB".vsrch.postChimera.p97clust.seqs.qza --output-path "$LIB".p97clust.seqs
qiime tools export --input-path "$LIB".chimerafilt-stats.qzv --output-path "$LIB".chimerafilt-stats.txt
 ```

## Combining datasets
Following these three filtering pipelines, the per-`$LIB` tables and representative sequences were joined into individual `.qza` artifacts. See the `seqfilter.combineNfilter.sh` shell script for full details. As an example, this is how the Dada2-filtered sequencing batches of tables and sequences were combined (the same would apply for deblur and vsearch, respectively):

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

The resulting `"$METHOD".hostseqs.fasta` file was searched with NCBI blast against nr database to determine if putative host hits were mismatches. Nearly all matches regardless of filtering method returned the same host for putative host-associated sequence variants: the Little Brown bat, _Myotis lucifugus_. There were two exceptions: one sequence variant was matched to the Eastern Small-footed bat, _Myotis leibii_ in each of the three filtering methods, while the Vsearch filtering strategy matched a sequence variant to the Evening bat, _Nycticeius humeralis_ (this variant was not present in either the DADA2 or Deblur pipeline).

The `*.droplist` file is then used to exlcude the ASVs we suspect are derived from bat hosts; any further sequence variants are assumed to be derived from arthropod (or otherwise will remain unclassified):
```
  qiime feature-table filter-features \
  --i-table "$METHOD".all.raw.table.qza \
  --m-metadata-file "$METHOD".droplist \
  --o-filtered-table "$METHOD".arthtable.qza \
  --p-exclude-ids
```

The `*.arthtable.qza` files produced represent our host-filtered datasets, one per filtering program (Dada2, Deblur, or Vsearch). Each `.qza` file was uploaded to the Github repo.

# QIIME artifact filtering and HashID amendments
In the midst of the sequence variant assignment process it was discovered that the QIIME implementation of the hashing algorithms used to assign each unique sequence some string of alphanumeric characters differed between Vsearch and Dada2/Deblur. We reran the vsearch algorithm natively and specified the `--relable_md5` hashID parameter to generate comparable sequence variants across platforms. See [this QIIME forum post](https://forum.qiime2.org/t/where-does-hashid-get-assigned-in-dada2-deblur-and-vsearch/7738) highlighting the observation.

We manually recreated the data tables by deconstructing the QIIME2 artifact file in R, amended the labels, and recreated the QIIME file. See the `qiime2R_datawrangling.R` script within the data/r_scripts section of the repo. Completion of that R script serves as the input for figures in the manuscript -  a single `.csv` file that is essentially a 'wide to long' transformation of the filtering program matricies, with each row representing one element of an OTU table (a per sample, per sequence variant count) which includes the associated metadata.
