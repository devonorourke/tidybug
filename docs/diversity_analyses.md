# Overview
Sequence data filtered as described in the `sequence_filtering.md` document was imported for all diversity estimate. We were interested in evaluating alpha and beta diversity estimates for mock community and guano data. The descriptions below briefly describe the R scripts used to achieve the specific plots and tables used in the mansucript.

## Alpha diversity
Scripts used to generate the associated figures and datasets for mock community analyses:
- The iNEXT package was used to generate the species accumulation curve presented for mock community samples. The `8_mock_iNEXTspecaccum.R` script was used to generate data for **Figures 4** and **6**.
- Alpha diversity estimates for a range of Hill Numbers (q=0, q=1, and q=2) for mock community data were generated using the `10_AlphaDiversity_HillNumbers_MOCK.R` script. These data were used to generate **Figure 5**.
- ANOVA and Dunn tests for mock communities were assessed using the `10_AlphaDiversity_HillNumbers_MOCK_dunnAndANOVAcalculations.R` scripts. The outputs of this script generated **Tables S3-S8**.

Scripts used to generate the associated figures and datasets for guano data analyses:
- Hill Number estimates of species diversity for guano samples were calculated using the `10_AlphaDiversity_HillNumbers_guano.R` script, with resulting data used to generate **Figure 7**.
- ANOVA and Dunn tests for alpha diversity estimates of guano data were performed using the `10_HillNumbers_guano_AnovaAndDunn.R` script. Results of the ANOVA are presented in **Tables S9-S14**. Pairwise Dunn results are visualized in **Figures 8-10**.

## Beta diversity
Comparisons of intra-sample diversity estimates were performed for mock and select guano data; all data were rarefied to a sampling depth fo 5000 sequences per sample.
- PERMANOVA were run with the Vegan function `adonis`
  - ADONIS tests for mock community samples were calculated using `14_adonisValues_MOCKdata.R`. ANOVA tables presented in **Tables S16-18**.
  - ADONIS tests for select guano samples were calculated using `14_adonisValues_FOXdata.R`. ANOVA tables presented in **Tables S19-21**.
- Distance measures for mock data were generated using `13_betadist_mock.R` and plotted in **Figure 11**.

- We calculated dispersion estimates for distances for mock and select guano samples using Vegan function `betadisper`.
Disperion results for datasets were all non-significant for the main effects of Filtering Parameter and Denoising method and were not reported.
  - `16_MOCK-betadispersions.R` was used to generate dispersion estimates for mock data.
  - `16_FOX_betadispersions.R` was used to generate dispersion estimates for guano data.
- Scree plots for guano data were used to identify an appropriate number of dimensions to run NMDS to generate ordination plots; see `17_FOX_scree_and_NMDS.R` script resulted in **Figures 12-14**.
