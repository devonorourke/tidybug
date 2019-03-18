# Background
We classified the same dataset using four different approaches and three related database constructed under different parameters (all were derived from BOLD records). First, we used our dereplicated database described in the `database_analyses.md` document and applied a Vsearch and Blast approach to classification in QIIME2. In addition, we classified sequences using Jon Palmer's AMPTK program with the full database included in his program; this database includes chordate and arthropod sequences and is not the same database used in the comparisons of taxonomic completeness, etc. described in the `database_analyses.md` document. Finally, we classified our sequences using the BOLD API through the `bold` R package. Filtering parametres that deviated from default settings are discussed for each program below.

We selected the DADA2-filtered guano data presented in the `sequence_filtering.md` script because the mock community analyses demonstrated denoising approaches having lower false positive error rates than the cluster-based approach; in addition the DADA2 approach yielded more unique sequences (13,407) than the other denoising program, Deblur (7,963), which gives us a greater sample size (unique sequences) to compare across classification platforms. We did not apply any additional filters described in the `sequence_filtering.md` document. The `dada2.arthseqs.qza` file served as input for the Vsearch and Blast-based classification in QIIME2; this file was exported into fasta format for use in Jon Palmer's pipeline as follows:
```
## run: source activate qiime2-2018.11
qiime tools export --input-path dada2.arthseqs.qza --output-path dada2.arthASVs.fasta
```

I further converted the fasta file into a 2-column text file for import into R:
```
zcat dada2.arthASVs.fasta.gz | paste - - | sed 's/>//' | gzip --best > dada2.arthASVs.txt.gz
```

# Classification
## QIIME2 classifications
We used QIIME2 classification tools to assign taxonomy to our dataset using BLAST and VSEARCH programs as follows. The reference database files used in classification are described in the `database_construction.md` document, consisting of a file of sequences `boldCOI.derep.seqs.qza` and taxonomic information `boldCOI.derep.tax.qza`. We amended two default parameters: the minimum alignment  coverage and identity to be considered a match. We found that the default values are much to low for COI sequences; in our experience sequences across arthropod taxonomic Orders routinely share similarities above 80%. Because QIIME2's classification scripts incorporate a lowest common ancestor (LCA) approach, the effect of having a very low threshold for percent identity (or coverage) results in very few sequences retaining information at the Species, Genus, or even Family level. We increased the specificity in defining our successful matches; the tradeoff is that more sequences are discarded as "Unassigned", though the ones that are classified have a more complete degree of information per record.

For BLAST:
```
## note file paths reflect our servers
REFPATH=/mnt/lustre/macmaneslab/devon/guano/BOLDdb/custom/allBOLD
FASTAPATH=/mnt/lustre/macmaneslab/devon/guano/Data/tidy_paper/qiime_tables/dada2

qiime feature-classifier classify-consensus-blast \
--i-query "$FASTAPATH"/dada2.arthseqs.qza \
--i-reference-reads "$REFPATH"/boldCOI.derep.seqs.qza \
--i-reference-taxonomy "$REFPATH"/boldCOI.derep.tax.qza \
--p-strand both \
--o-classification derep.blast.qza
```

For VSEARCH:
```
## note file paths reflect our servers
REFPATH=/mnt/lustre/macmaneslab/devon/guano/BOLDdb/custom/allBOLD
FASTAPATH=/mnt/lustre/macmaneslab/devon/guano/Data/tidy_paper/qiime_tables/dada2

qiime feature-classifier classify-consensus-vsearch \
--i-query "$FASTAPATH"/dada2.arthseqs.qza \
--i-reference-reads "$REFPATH"/boldCOI.derep.seqs.qza \
--i-reference-taxonomy "$REFPATH"/boldCOI.derep.tax.qza \
--p-perc-identity 0.97 \
--p-query-cov 0.94 \
--p-strand both \
--p-threads 24 \
--o-classification derep.vsearch.qza
```

The resulting `.qza` files were exported as text files `derep.vsearch.pid97.gz` and `derep.blast.pid97.gz`.
```
qiime tools export --input-path derep.vsearch.qza --output-path classified_seqs
cd classified_seqs
mv taxonomy.tsv derep.vsearch.pid97
gzip --best derep.vsearch.pid97
```
> repeat for `blast` output

## AMPTK classification
Jon's default classification strategy is described on the [amptk site](https://amptk.readthedocs.io/en/latest/taxonomy.html), and uses a hybrid method involving global alignment, followed by a similar LCA implementation for matches above a defined threshold for sequences that surpass the similarity threshold (default is 0.7). Then, the same dataset is classified by running two Bayesian classifers trained on the same COI database, and the most complete taxonomy record produced from either classifier is retained. Classified sequences with less than 97% identity are assigned to the Bayesian classifier record, while those greater than 97% identity are assigned the global alignment taxonomy string. A final LCA value is then applied on each record retained among the remaining instance in which multiple "winners" are present for each group to produce the final dataset.
We included these Bayesian steps, but increased the default identities required for inclusion from the default 0.8 to the maximum allowed 0.9. In addition, We increased the default threshold for global alignment from 0.7 to match that used in the QIIME classifiers at 0.97. Notably, there is no query coverage option available in AMPTK.

We had previously installed the amptk Database - see his [documentation](https://amptk.readthedocs.io/en/latest/taxonomy.html#taxonomy-databases) for details on it's construction.

```
source activate amptk
# database installation command ==> 'amptk install -i COI'

## temporarily unzip:
gzip -d dada2.arthASVs.fasta.gz

amptk taxonomy \
--fasta dada2.arthASVs.fasta.gz \
--utax_cutoff 0.9 --usearch_cutoff 0.97 --sintax_cutoff 0.9 \
--db COI --method hybrid --cpus 24

## zip up:
gzip --best dada2.arthASVs.fasta
```

The fasta headers in the classified fasta file were reformatted for import into R:
```
grep '^>' dada2.arthASVs.otus.taxonomy.fa | sed 's/>//' | cut -f 1 -d ' ' > tmp.left
grep '^>' dada2.arthASVs.otus.taxonomy.fa | sed 's/>//' | cut -f 2- -d ' ' > tmp.right
paste tmp.left tmp.right | gzip --best > amptk.taxa.txt.gz
```

## BOLD-API classification
Querying records through BOLD directly is perhaps the most frequent means with which arthropod COI data are classified because it doesn't require any database construction. Several programs are available to access the BOLD API to query the representative sequences for potential BOLD reference matches through R (ex. see the "bold" package) and Python (ex. see "bold-retriever") scripts. We used the R `bold` package to perform this task in the `boldAPI_classification.R` script.


## Subsequent VSEARCH filtering - changing query coverage
To test whether the addition of a query coverage term altered the number of ASVs that would ultimately be classified we ran standalone VSEARCH classification with three pairs of percent identity and query coverage parameters:

- (p97c94) 97% identity, 94% coverage; matches the parameters used in QIIME2's VSEARCH and BLAST classifiers we employed
- (p97c19) 97% identity, 10% coverage; those included here but not in p97c94 would represent the additional taxa identified by AMPTK but not VSEARCH or BLAST
- (p90c10) 90% identity, 10% coverage; we lowered the percent identity to classify a larger fraction of ASVs that fell below the 97% threshold - these are likely to be passed to UTAX or SINTAX in the AMPTK method

We applied a similar function for each pair, changing the `--id` and `--query_cov` values as appropriate. Here's an example for p97c94:
```
vsearch --usearch_global $FASTA --db $REF --threads 24 \
--id 0.90 --query_cov 0.1 --strand both --maxaccepts 0 --maxrejects 0 --maxhits 1 \
--output_no_hits --blast6out vsearch.p90c10out.txt
```
> `$FASTA` input is the DADA2-processed fasta file `dada2.arthASVs.fasta`; `$REF` represents the dereplicated BOLD dataset `boldCOI.derep.fasta.gz`

Outputs were imported into `classifier_analysesPalmer-Vsearch_BayesianClassify_sub90pid.R` for analyses.


## Chimera filtering
Chimera's were evaluated with Vsearch:
```
## $REPFASTA is the dada2-processed arthropod sequence file "dada2.arthASVs.fasta"
## $REFFASTA is the dereplicated BOLD COI database fasta "boldCOI.derep.fasta.gz"

vsearch --threads 24 \
    --uchime_ref "$REPFASTA" \
    --db $REFFASTA \
    --sizein --sizeout --fasta_score --fasta_width 0 \
    --nonchimeras derep.ref.nonchimeras.fasta \
    --chimeras derep.ref.chimeras.fasta
```

The resulting reference-based chimera filtering was exported into R:
```
grep "^>" derep.ref.chimeras.fasta | sed 's/>//' > chimera.ref.headers.txt
```

# Output analyses:
See the `classification_analyses.R`, `classifier_UnassignedASVsFrequencies`, and `classifier_analysesPalmer-Vsearch_BayesianClassify_sub90pid.R` scripts for details on tables and plots.
