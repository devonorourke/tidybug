## example shell script used to create the Deblur-filtered output used in the analyses of quality filtering pipelines
## note this script was modified for each library

#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p81/qiime/reads
#SBATCH --partition=macmanes,shared
#SBATCH --job-name="p81.dblr"
#SBATCH --output=p81.deblur.log
#SBATCH --cpus-per-task=24

module purge
module load anaconda/colsa
source activate qiime2-2018.11

LIB=$(pwd | cut -d '/' -f 9)

mkdir deblur
cd deblur

## denoise with Deblur using custom BOLD database of arthropod-only records
## note: additional scripts available explaining how database was formed

REF=/mnt/lustre/macmaneslab/devon/guano/BOLDdb/ref_seqs.qza
qiime deblur denoise-other \
--i-demultiplexed-seqs ../"$LIB".qfiltd.seqs.qza \
--i-reference-seqs "$REF" \
--p-trim-length 181 \
--p-sample-stats \
--p-jobs-to-start 24 \
--o-table "$LIB".dblr.table.qza \
--o-representative-sequences "$LIB".dblr.repseqs.qza \
--o-stats "$LIB".dblr.stats.qza

## export data in text format:
qiime deblur visualize-stats --i-deblur-stats "$LIB".dblr.stats.qza --o-visualization "$LIB".dblr.stats.qzv
qiime tools export --input-path "$LIB".dblr.stats.qzv --output-path "$LIB".dblr.output
qiime tools export --input-path "$LIB".dblr.repseqs.qza --output-path "$LIB".dblr.output
qiime tools export --input-path "$LIB".dblr.table.qza --output-path "$LIB".dblr.output

## all resulting `"$LIB".dblr.table.qza` and `"$LIB".dblr.repseqs.qza` files (1 pair for each $LIB) combined..
## ..for processing in R script