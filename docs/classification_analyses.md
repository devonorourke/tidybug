# Classifiers and databases setup
We compared the taxonomic information provided by five different classifiers. The _tidybug_ database described in the [database_construction.md](https://github.com/devonorourke/tidybug/blob/master/docs/database_construction.md) document served as the reference source to classify mock and guano data for four of the five classifiers:
1. QIIME 2 implementation of VSEARCH (global aligner)
2. QIIME 2 implementation of BLAST (local aligner)
3. QIIME 2 implementation of Naive Bayes (kmer-based machine learning)
4. VSEARCH implementation of SINTAX (kmer-based)

However, note that the input files varied depending on the classifiers used:
A. BLAST and VSEARCH used the QIIME 2 formatted `boldCOI.derep.seqs.qza` and `boldCOI.derep.tax.qza` files.
B. The Naive Bayes classifier required those `boldCOI.derep.seqs.qza` and `boldCOI.derep.tax.qza` files to serve as inputs to generate a trained classifier. The code used to train the classifier is documented below, and resulted in a `nbClassifer_allboldDB.qza` file.
C. The SINTAX classifier used a modified version of the `boldCOI.derep.fasta` file; namely, we needed to append taxonomy information to the header using the information matching for each header ID in the `boldCOI.arth_derep.txt` file. The resulting `tmp_forSINTAX.fasta.gz` file was used for that classifier. Details for generating this file are provided below.

All of these files are available at the [OSF repo for this project](https://osf.io/k3eh6/files/).

The fifth classifier relied on the BOLD API [Taxonomy Engine](http://v4.boldsystems.org/index.php/resources/api?type=taxonomy), which does not provide details about specific parameters invoked during classification, nor does it specify the content of the database used for it's classification. Thus, while the tidybug database was derived from BOLD references, it's unclear whether the contents of our database were identical to those references that the BOLD API used. We classified our sequences using the BOLD API through the `bold` R package. Notably, the output of the Taxonomy API is not a single match like with the other classifiers, rather, it is a list of up to 100 potential "hits". We therefore and applied specific filtering parameters to this list to best match those of the other classifiers when possible. Specifically, we retained matches only above either 95, 97, or 99% identity. Because query coverage was not a parameter that was provided in the output, we applied an LCA process to those remaining hits. See the [boldAPI_classification_mockSamples.R](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/boldAPI_classification_mockSamples.R) and [boldAPI_classification.R](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/boldAPI_classification.R) R scripts for details about filters applied to the mock and bat guano data, respectively.

# Classification Parameters
Mock samples and guano samples were processed separately using the parameters defined below. In both cases we used a QIIME 2 environment to classify samples with VSEARCH, BLAST, and Naive Bayes programs. The BOLD API required a separate workflow described below. The SINTAX classifier was implemented using VSEARCH-v2.13.4.

## Mock samples

As mentioned in the methods section of the manuscript, we amended added species names to two of the 24 mock specimen which lacked information based on sequence alignment searches using both GenBank and BOLD databases. The remaining specimen that did not contain high percent alignment identity and coverage was removed from the analysis for classification - this modified file is available with taxonomic information in the fasta header [at this link](https://github.com/devonorourke/tidybug/blob/master/data/classify_comps/mock_comps/CFMR_insect_mock4_wtax_noIM52.fasta). The fasta format was converted into a tab delimited file containing the updated names and is available [at this link](https://github.com/devonorourke/tidybug/tree/master/data/classify_comps/mock_comps/expected_mockData/mock_expected_taxa.txt).

The mock fasta was converted into a format useful for QIIME 2. We stripped out the taxonomic names from the headers of the fasta file first, then imported that file:
```
## strip headers
cat CFMR_insect_mock4_wtax_noIM52.fasta | paste - - | cut -f 1,4 | sed 's/;//g' | tr '\t' '\n' > CFMR_insect_mock4_noIM52.fa

## import to QIIME 2 format
qiime tools import \
  --input-file CFMR_insect_mock4.fa --output-file mockIM4_noIM42_seqs.qza \
  --input-formate FeatureData[Sequence]
```

The `.qza` file was used as input for each of the QIIME 2-based classifiers, while the `.fa` version was used in the SINTAX and BOLD API inputs.

---

**BLAST**: Modified the number of accepted hits from default (10) to 1000. Would have set to 0, except BLAST doesn't allow for this parameter in the QIIME implementation. Percent query coverage (0.89) was maintained across all possible parameters, while the percent identity parameter was explored using three values: 0.95, 0.97, and 0.99. Thus, we altered the `--p-perc-identity` argument below:

```
qiime feature-classifier classify-consensus-blast \
  --i-query mockIM4_seqs.qza \
  --i-reference-reads boldCOI.derep.seqs.qza \
  --i-reference-taxonomy boldCOI.derep.tax.qza \
  --p-maxaccepts 1000 \
  --p-perc-identity 0.97 \
  --p-query-cov 0.89 \
  --p-strand both \
  --o-classification mockIM4.blast_out.qza
```

**VSEARCH**: The parameters used here matched the BLAST settings for similar number of top hits to retain and same percent query coverage. Likewise, we explored the same set of percent identities.
```
qiime feature-classifier classify-consensus-vsearch \
  --i-query mockIM4_seqs.qza \
  --i-reference-reads boldCOI.derep.seqs.qza \
  --i-reference-taxonomy boldCOI.derep.tax.qza \
  --p-maxaccepts 1000 \
  --p-perc-identity 0.97 \
  --p-query-cov 0.89 \
  --p-strand both \
  --p-threads 12 \
  --o-classification mockIM4.vsearch_out.qza
```

**VSEARCH with "top hits"**: We added one other parameter to the VSEARCH classification by invoking a `--top-hits` parameter using the same code as above. This parameter retains only the top alignment score among those hits that are above the assigned percent identity and coverage. In the event of ties an LCA process will be applied.

**Naive Bayes**: Another kmer-based classification method used the SciKit Learn program wrapped in QIIME 2's `Naive Bayes` classifier. This required first training the database (the same dataset used for the VSEARCH and BLAST classifiers).
```
READS=mockIM4_seqs.qza
CLSSFYR=nbClassifer_allboldDB.qza

qiime feature-classifier classify-sklearn \
  --i-reads "$READS" --i-classifier "$CLSSFYR" \
  --p-n-jobs 1 --p-reads-per-batch 2000 \
  --o-classification mockIM4.sklearn_out.qza \
  --verbose
```

Once the classifier (`nbClassifer_allboldDB.qza`) was trained, we then invoked the Naive Bayes process to the mock dataset. We amended a confidence parameter similar to what was invoked with SINTAX (see below) using the same range of values: 0.3, 0.5, 0.7, 0.8, and 0.9. Thus, we altered the `--p-confidence` as needed (value of 0.3 adjusted in subsequent runs):
```
qiime feature-classifier classify-sklearn \
  --i-reads mockIM4_seqs.qza --i-classifier nbClassifer_allboldDB.qza \
  --p-n-jobs 1 --p-reads-per-batch 2000 \
  --o-classification mockIM4.nbayes_out.qza \
  --p-confidence 0.3 \
  --verbose
```

**SINTAX**: Unlike alignment-based classifiers, SINTAX does not have any percent identity or query coverage parameters. The one change to default parameters was setting the bootstrap confidence value. We explored a range of values including 0.3, 0.5, 0.7, 0.8 (default), and 0.9. Note that we're using the VSEARCH implementation of SINTAX, not the version created by Robert Edgar. The database used for taxonomic evaluation is the same used in the other classifiers, with headers being modified from their implementation in QIIME (in QIIME the taxonomy and fasta sequences are in separate files; here the taxonomy information is part of the fasta header). These are the _same_ deprelicated, LCA-processed taxa used in QIIME-based classifiers. We reformatted our tidybug database data so that a single fasta contained the sequence and taxonomy information as follows:
```
## Take the dereplicated fasta file and split into 2 column:
zcat boldCOI.derep.fasta.gz | paste - - | sort -k1,1 > tmpfasta2merge.txt

## Take the metadata file and reformat to match amptk record headers:
zcat boldCOI.arth_derep.txt.gz | sed 's/^/>/g' | sed 's/;/,/g' | sed 's/__/:/g' | sed 's/;k:/;tax=k:/g' | sort -k1,1 > tmpamptkheaders.txt

## Join the two files by the common field with the sequence record ID and substitute spaces to create proper fasta format:
join -1 1 -2 1 tmpamptkheaders.txt tmpfasta2merge.txt -t $'\t' | sed 's/\t/;tax=/' | sed 's/\t/\n/' | sed 's/ /_/g' > tmp_forSINTAX.fasta.gz
```

We then used this file for classifying with SINTAX using [this VSEARCH release](https://github.com/torognes/vsearch/releases/tag/v2.13.4):
```
$VSEARCH \
--sintax CFMR_insect_mock4.fa \
--db tmp_forSINTAX.fasta.gz \
--tabbedout vsintax_out.txt \
--sintax_cutoff 0.7 \
--threads 12
```

**BOLD-API**: Querying records through BOLD directly is perhaps the most frequent means with which arthropod COI data are classified because it doesn't require any database construction. Several programs are available to access the BOLD API to query the representative sequences for potential BOLD reference matches through R (ex. see the "bold" package) and Python (ex. see "bold-retriever") scripts. We used the R `bold` package to perform this task in the [boldAPI_classification_mockSamples.R](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/boldAPI_classification_mockSamples.R) R script.

Each of QIIME `.qza` artifacts were then exported as `.txt` files. SINTAX and BOLD-API data and the QIIME classifier outputs are available in the [classify_comps directory of the repo](https://github.com/devonorourke/tidybug/tree/master/data/classify_comps/mock_comps).
> Note that each of these classifier files was further manipulated for use in the classifier comparisons such that there was one taxonomy string and a relative abundance of that classified sequence.

The resulting text files were used as inputs for the [figure6_classifier_MockComps.R](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/figure6_classifier_MockComps.R) R script to generate [Figure6](https://github.com/devonorourke/tidybug/blob/master/figures/figure6_classifierComps_wVals.png).

## Guano samples
We selected the DADA2-filtered guano data presented in the `sequence_filtering.md` script because the mock community analyses demonstrated denoising approaches having lower false positive error rates than the cluster-based approach; in addition the DADA2 approach yielded more unique sequences (13,407) than the other denoising program, Deblur (7,963), which gives us a greater sample size (unique sequences) to compare across classification platforms. We did not apply any additional filters described in the `sequence_filtering.md` document (i.e. this is the "basic" filter). The `dada2.arthseqs.qza` file served as input for the VSEARCH, BLAST, and Naive Bayes classifiers in QIIME2. The `.qza` file was exported into fasta format for evaluation with SINTAX and the BOLD API:
```
## export as fasta
qiime tools export --input-path dada2.arthseqs.qza --output-path dada2.arthASVs.fasta
```

Both the `.qza` and `.fasta` files are available in the [data/qiime directory of the Tidybug repo](https://github.com/devonorourke/tidybug/tree/master/data/qiime).

---

**BLAST**: As with the mock samples, we matched parameters between VSEARCH and BLAST.
```
qiime feature-classifier classify-consensus-blast \
--i-query dada2.arthseqs.qza \
--i-reference-reads boldCOI.derep.seqs.qza \
--i-reference-taxonomy boldCOI.derep.tax.qza \
--p-maxaccepts 1000 \
--p-perc-identity 0.97 \
--p-query-cov 0.89 \
--p-strand both \
--o-classification derep.blast.qza
```

**VSEARCH**: Parameters match BLAST. Not we did not invoke the "top-hits" parameter explored in the mock community samples.
```
## note file paths reflect our servers
qiime feature-classifier classify-consensus-vsearch \
--i-query dada2.arthseqs.qza \
--i-reference-reads boldCOI.derep.seqs.qza \
--i-reference-taxonomy boldCOI.derep.tax.qza \
--p-maxaccepts 1000 \
--p-perc-identity 0.97 \
--p-query-cov 0.89 \
--p-strand both \
--p-threads 24 \
--o-classification derep.vsearch.qza
```

**Naive Bayes**: Same classifier trained in mock data used here; same parameters also, using a 0.7 value for confidence (default).
```
qiime feature-classifier classify-sklearn \
  --i-reads dada2.arthseqs.qza --i-classifier nbClassifer_allboldDB.qza \
  --p-n-jobs 1 --p-reads-per-batch 2000 \
  --o-classification guano.nbayes_out.qza \
  --verbose
```

**SINTAX**
```
$VSEARCH \
--sintax dada2.arthASVs.fasta.gz \
--db tmp_forSINTAX.fasta.gz \
--tabbedout vsintax_out.txt \
--sintax_cutoff 0.7 \
--threads 12
```

**BOLD-API**: We used a tab delimited version of the DADA2-processed guano fasta `dada2.arthASVs.txt.gz` - available in the same [data/qiime directory of the Tidybug repo](https://github.com/devonorourke/tidybug/tree/master/data/qiime). We used this file to query the BOLD API in batches using the [boldAPI_classification.R](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/boldAPI_classification.R) R script.

The resulting `.qza` and `.txt` files containting taxonomic information were saved in the [/data/classify_comps/guano_comps section of the Repo](https://github.com/devonorourke/tidybug/tree/master/data/classify_comps/guano_comps). These files were processed with the [figure7_classifier_GuanoComparisons.R](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/figure7_classifier_GuanoComparisons.R) R script to generate [Figure 7](https://github.com/devonorourke/tidybug/blob/master/figures/figure7_classifierComps.png). The same data approach was used to generate [Figure S4](https://github.com/devonorourke/tidybug/blob/master/SupplementaryFiguresTables/figureS4_classifier_consequence.png) following the [figureS4_classifier_Consequence.R](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/figureS4_classifier_Consequence.R) script.
