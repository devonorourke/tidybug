Supplementary data tables and figures referenced in manuscript available within this page.

- Figure captions are embedded within each .png file.
- Table captions associated with each dataset are available below:

**TableS2**:
Alpha diversity estimates among four mock community samples. Number of sequences or sequence equivalents are provided for each combination of denoising program (Denoiser) and filtering regime (basic, standard, or extra) for each Hill Number (q = 0 is equivalent to observed richness; q = 1 is equivalent to Shannon’s Entropy; q = 2 is equivalent to Simpsons’s 1-D diversity).

**TableS3**:
Kruskal-Wallis statistic and Benjamini-Hochberg adjusted significance values for bat guano data. 

**TableS4**:
Dunn’s test for pairwise differences among denoising groups alpha diversity estimates for bat guano data.

**TableS5**:
Permutational Multivariate Analysis of Variance tests using guano samples collected from one location (Fox State Forest, Hillsborough NH) between April through October 2016. PERMANOVA run using three distance inputs: unweighted abundance metric (Dice-Sorensen) and weighted abundance metrics (Bray-Curtis and Morisita-Horn) on rarefied samples testing effects of denoising method (Method), filtering parameter (Filt), and date of sample collection (MonthStart).

**TableS6**:
Summary of per-library shared ASVs among denoising programs. Each element of the matrix represents the number of ASVs observed. Rows 2-4, 6-8, and 10-12 represent a unique combination of mock library sequenced (libA-libD) and whether or not the ASVs were Exact matches (100% identity) to expected mock sequence, Partial matches (>97% identity, not exact), or Miss (< 97% identity). Rows 5, 9, and 13 represent the number of ASVs shared across all four mock libraries for a given denoising/filtering pipeline. Columns 2-4, 6-8, and 10-12 represent the number of distinct ASVs observed for that mock sample, for a given denoising pipeline and filtering parameter. “Basic” represents default parameters for each filtering method; “Standard” requires a sample to have > 5000 reads, and an OTU to be present in > 1 sample; “Extra” includes “Standard” filters in addition to subtracting a fixed number of reads from all observations. CommonASVs represent the number of ASVs shared among all three of the denoising pipelines at a common filtering threshold. For example in the fifth colummn, CommonASVs-basic, there were 19 ASVs that were exact matches to the reference sequences shared across all three denoising pipelines for mock LibraryA when filtered with the Basic parameter, yet just 1 common ASV shared among the miss collection in the same library with the same filtering parameters.
