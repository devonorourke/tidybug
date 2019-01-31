## example shell script used to create the Vsearch-filtered output used in the analyses of quality filtering pipelines
## note this script was modified for each library

#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/Data/individLibs/p81/qiime/reads
#SBATCH --partition=macmanes,shared
#SBATCH --job-name="p81.vsrch"
#SBATCH --output=p81.vsrch.log
#SBATCH --cpus-per-task=24

module purge
module load anaconda/colsa
source activate qiime2-2018.11

LIB=$(pwd | cut -d '/' -f 9)

mkdir vsearch
cd vsearch

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
