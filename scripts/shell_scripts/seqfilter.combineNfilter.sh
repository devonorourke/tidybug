## script combines the 4 libraries worth of sequences and tables per filtering method
## script takes these per-filter-method tables/seqs and filters bat host DNA
## resulting output used as input for R scripted figures used in publication
## short R script used to generate a `.tsv` format of the filtered tables imported to generate figures shown here also
## input created from `seqfilter.deblur_examples.sh`, `seqfilter.vsearch_examples.sh`, and `seqfilter.dada2_examples.sh`

#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/Data/methods_paper
#SBATCH --partition=macmanes,shared
#SBATCH --job-name="mergeNfilt"

module purge
module load anaconda/colsa
source activate qiime2-2018.11

##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ #####
## part 1: merging tables and sequences
##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ #####

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


## repeat merging for deblur tables and seqs
# tables
qiime feature-table merge \
  --i-tables /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p41/qiime/reads/deblur/p41.dblr.table.qza \
  --i-tables /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p42/qiime/reads/deblur/p42.dblr.table.qza \
  --i-tables /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p71/qiime/reads/deblur/p71.dblr.table.qza \
  --i-tables /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p72/qiime/reads/deblur/p72.dblr.table.qza \
  --o-merged-table deblur.all.raw.table.qza
# sequences
qiime feature-table merge-seqs \
  --i-data /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p41/qiime/reads/deblur/p41.dblr.repseqs.qza \
  --i-data /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p42/qiime/reads/deblur/p42.dblr.repseqs.qza \
  --i-data /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p71/qiime/reads/deblur/p71.dblr.repseqs.qza \
  --i-data /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p72/qiime/reads/deblur/p72.dblr.repseqs.qza \
  --o-merged-data deblur.all.raw.seqs.qza

## repeat merging tables for vsearch tables and seqs
# tables
qiime feature-table merge \
  --i-tables /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p41/qiime/reads/vsearch/p41.vsrch.postChimera.p97clust.table.qza \
  --i-tables /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p42/qiime/reads/vsearch/p42.vsrch.postChimera.p97clust.table.qza \
  --i-tables /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p71/qiime/reads/vsearch/p71.vsrch.postChimera.p97clust.table.qza \
  --i-tables /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p72/qiime/reads/vsearch/p72.vsrch.postChimera.p97clust.table.qza \
  --o-merged-table vsearch.all.raw.table.qza
# sequences
qiime feature-table merge-seqs \
  --i-data /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p41/qiime/reads/vsearch/p41.vsrch.postChimera.p97clust.seqs.qza \
  --i-data /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p42/qiime/reads/vsearch/p42.vsrch.postChimera.p97clust.seqs.qza \
  --i-data /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p71/qiime/reads/vsearch/p71.vsrch.postChimera.p97clust.seqs.qza \
  --i-data /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p72/qiime/reads/vsearch/p72.vsrch.postChimera.p97clust.seqs.qza \
  --o-merged-data vsearch.all.raw.seqs.qza

##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ #####
## part 2: remove bat (host) sequences from table/seq files
##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ #####

# Filtering host reads
## See 'host_seq_explanation.md' file for information on how the `$REF` sequences and taxonomies were created and imported

REF=/mnt/lustre/macmaneslab/devon/guano/BOLDdb/NHbat_seqs.qza

for METHOD in $(ls -1 | cut -f 1 -d '.' | sort -u); do
  qiime quality-control exclude-seqs \
    --i-query-sequences "$METHOD".all.raw.seqs.qza \
    --i-reference-sequences "$REF" \
    --p-perc-identity 0.95 \
    --o-sequence-hits "$METHOD".hostseqs.qza \
    --o-sequence-misses "$METHOD".arthseqs.qza \
    --p-method vsearch;

  qiime tools export --input-path "$METHOD".hostseqs.qza --output-path "$METHOD"_hostseqs;

  grep "^>" "$METHOD"_hostseqs/dna-sequences.fasta | sed 's/>//' | sed '1 i\#OTU ID' > "$METHOD".droplist;

  qiime feature-table filter-features \
  --i-table "$METHOD".all.raw.table.qza \
  --m-metadata-file "$METHOD".droplist \
  --o-filtered-table "$METHOD".arthtable.qza \
  --p-exclude-ids

  qiime tools export --input-path "$METHOD".hostseqs.qza --output-path "$METHOD"_hostseqs;
  cd "$METHOD"_hostseqs;
  mv dna-sequences.fasta "$METHOD".hostseqs.fasta;
  mv "$METHOD".hostseqs.fasta ..;
  cd ..;
  rm -r "$METHOD"_hostseqs;
done


## Each `"$METHOD".hostseqs.fasta` file was searched with NCBI blast. See `host_seq_explanation.md` file for further details.
## The `"$METHOD".arthtable.qza` files represent the inputs to convert from QIIME-format into that used in R to generate figures in this publication


##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ #####
## part 3: convert the "$METHOD".arthtable.qza files into .tsv files to load into R for plotting
##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ ##### ~ #####

see `qiime2R_datawrangling.R` script for complete details
