# Overview
Sequence data was obtained and filtered as described in the [sequence_filtering.md](https://github.com/devonorourke/tidybug/blob/master/docs/sequence_filtering.md) document. As explained at the end of that document, the [sequence_filtering.R script](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/1_sequence_filtering.R) was applied to create a master file containing the necessary data structures for the relevant denoising comparisons and diversity estimates - this is the `all.filtmethods.df.csv.gz` file available in [this directory of the GitHub repo](https://github.com/devonorourke/tidybug/raw/master/data/text_tables). 

We imported the `all.filtmethods.df.csv.gz` file into a new R script to generate Figure 1
applying the  imported for all diversity estimate. We were interested in evaluating alpha and beta diversity estimates for mock community and guano data. The descriptions below briefly describe the R scripts used to achieve the specific plots and tables used in the manuscript.

## Alpha diversity
Scripts used to generate the associated figures and datasets for mock community analyses:
- The iNEXT package was used to generate the species accumulation curve presented for mock community samples. The [8_mock_iNEXTspecaccum.R](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/8_mock_iNEXTspecaccum.R) script was used to generate data for **Figure 4**
- Alpha diversity estimates for a range of Hill Numbers (q=0, q=1, and q=2) for mock community data were generated using the [10_AlphaDiversity_HillNumbers_MOCK.R](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/10_AlphaDiversity_HillNumbers_MOCK.R) script. These data were used to generate **Figure S4**.
- ANOVA and Dunn tests for mock communities were assessed using the [10_AlphaDiversity_HillNumbers_MOCK_dunnAndANOVAcalculations.R](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/10_HillNumbers_guano_AnovaAndDunn.R) script. The pairwise Dunn data are presented in **Tables S3-S5**, though text files are available for the [Dunn](https://github.com/devonorourke/tidybug/tree/master/data/text_tables/dunn_mock), [ANOVA](https://github.com/devonorourke/tidybug/tree/master/data/text_tables/anova_mock) and [Kruskal-Wallis](https://github.com/devonorourke/tidybug/tree/master/data/text_tables/kw_mock) outputs for each Hill value.

Scripts used to generate the associated figures and datasets for guano data analyses:
- Hill Number estimates of species diversity for guano samples were calculated using the [10_AlphaDiversity_HillNumbers_guano.R](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/10_AlphaDiversity_HillNumbers_guano.R) script, with resulting data used to generate **Figures 5 and S5**.
- ANOVA and Dunn tests for alpha diversity estimates of guano data were performed using the [10_AlphaDiversity_HillNumbers_GUANO_dunnAndANOVAcalculations.R](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/10_AlphaDiversity_HillNumbers_GUANO_dunnAndANOVAcalculations.R) script. Results of the pairwise Dunn tests are presented in **Tables S6-S8**. Text files for [Dunn](https://github.com/devonorourke/tidybug/tree/master/data/text_tables/dunn_guano), [Anova](https://github.com/devonorourke/tidybug/tree/master/data/text_tables/anova_guano), and [Kruskal-Wallis](https://github.com/devonorourke/tidybug/tree/master/data/text_tables/kw_guano) data are also available.

## Beta diversity
Comparisons of intra-sample diversity estimates were performed for mock and select guano data; all data were rarefied to a sampling depth fo 5000 sequences per sample.
- PERMANOVA were run with the Vegan function `adonis`
  - ADONIS tests for mock community samples were calculated using [14_adonisValues_MOCKdata.R](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/14_adonisValues_MOCKdata.R). ANOVA tables presented in **Tables S9-11**.
  - ADONIS tests for select guano samples were calculated using [14_adonisValues_FOXdata.R](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/14_adonisValues_FOXdata.R). ANOVA tables presented in **Tables S12-14**.
- Distance measures for mock data were generated using [13_betadist_mock.R](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/13_betadist_mock.R). These were not reported in this manuscript.

- We calculated dispersion estimates for distances for mock and select guano samples using Vegan function `betadisper`.
Disperion results for datasets were all non-significant for the main effects of Filtering Parameter and Denoising method and were not reported.
  - [16_MOCK-betadispersions.R](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/16_MOCK-betadispersions.R) was used to generate dispersion estimates for mock data.
  - [16_FOX_betadispersions.R](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/16_FOX_betadispersions.R) was used to generate dispersion estimates for guano data.
- Scree plots for guano data were used to identify an appropriate number of dimensions to run NMDS to generate ordination plots; see the [17_FOX_scree_and_NMDS.R](https://github.com/devonorourke/tidybug/blob/master/scripts/R_scripts/17_FOX_scree_and_NMDS.R) script for code used to generated the ordinations presented in **Figures S6-S8**.
