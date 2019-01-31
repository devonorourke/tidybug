## example shell script used to create the filtered reads used in the Vsearch and Deblur analyses of quality filtering pipelines
## note this script was modified for each library

#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p81/qiime/reads
#SBATCH --partition=macmanes,shared
#SBATCH --job-name="p81.qfilt"
#SBATCH --output=p81.qfilt.log

module purge
module load anaconda/colsa
source activate qiime2-2018.11

LIB=$(pwd | cut -d '/' -f 9)

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

## "$LIB".qfiltd.seqs.qza used for both Deblur and Vsearch inputs
## see `seqfilter.deblur_example.sh` and `seqfilter.vsearch_example.sh` for further processing
