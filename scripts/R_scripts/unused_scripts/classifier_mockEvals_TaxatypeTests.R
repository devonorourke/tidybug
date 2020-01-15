library(tidyverse)
library(reshape2)

################################################################################
## note on updating this plot:
#0. may need to switch the data table used for "expected" as unclear what current ground truth is
#4. need to alter plot code so taxa Level are ordered properly (via levels=factor()...)
################################################################################


################################################################################
## data wrangling
################################################################################

## read in Data
readerfunction <- function(urlpath) {
  tmp <- read_delim(file = urlpath, delim = "\t", col_names = FALSE) %>%
    rename(., IMalias = X1, Taxon = X2) %>% 
    separate(., 
             col = Taxon, 
             into = c("kingdom", "phylum", "Class", "Order", "Family", "Genus", "Species"), 
             sep = ";") %>% 
    select(-kingdom, -phylum)
  tmp <- as.data.frame(apply(tmp, 2, function(y) (gsub(".__", "", y))))
  tmp <- as.data.frame(apply(tmp, 2, function(y) (gsub("^$|^ $", NA, y))))
  tmp <- as.data.frame(apply(tmp, 2, function(y) (gsub("^NA$", NA, y))))
  tmp[] <- lapply(tmp, as.character)
  tmp %>% gather(key = level, value = Taxon, -IMalias)
}

exp_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/expected_mockData/mock_expected_taxa.txt"
bl_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/blast/mock_blast_observedTaxa.txt"
vs_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/vsearch/mock_vsearch_observedTaxa.txt"
st_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/sintax/mock_sintax_observedTaxa.txt"
nb_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/nbayes/mock_nbayes_observedTaxa.txt"
bo_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/bold/boldAPI_mock_lca97.txt"

exp_df <- readerfunction(exp_url)
bl_df <- readerfunction(bl_url)
vs_df <- readerfunction(vs_url)
st_df <- readerfunction(st_url)
nb_df <- readerfunction(nb_url)
bo_df <- readerfunction(bo_url)

## evaluate the FalsePositive, FalseNegative, TruePositive, and TrueNegative conditions for each taxa, at each Level
## merge two datasets:
truthfunction <- function(obs_data, classifier) {
  selectTests <- c("TAR", "TDR", "F")
  inner_join(exp_df, obs_data, by = c("IMalias", "level"), suffix = c("_exp", "_obs")) %>% 
    mutate(truth = case_when(
      Taxon_exp == Taxon_obs  ~ "TP",
      is.na(Taxon_exp) == is.na(Taxon_obs)  ~ "TN",
      is.na(Taxon_exp) & !is.na(Taxon_obs)  ~ "FP",
      !is.na(Taxon_exp) & is.na(Taxon_obs)  ~ "FN")) %>% 
    select(-starts_with("Taxon")) %>%
    spread(key = level, value = truth) %>% 
    melt(., id.vars = "IMalias", variable.name = "TaxaLevel", value.name = "StateType") %>% 
    group_by(TaxaLevel, StateType) %>% 
    summarize(counts = n()) %>% 
    spread(StateType, counts, fill = 0) %>% 
    group_by(TaxaLevel) %>% 
    mutate(TAR = TP / (TP + FP)) %>% 
    mutate(TDR = TP / (TP + FN)) %>% 
    mutate(F = 2 * ((TAR*TDR)/(TAR+TDR))) %>% 
    melt(., id.vars = 'TaxaLevel', value.name = 'Score', variable.name = 'Metric') %>%
    filter(Metric %in% selectTests) %>% 
    mutate(Classifier = classifier)
}

bl_tmp_dat <- truthfunction(bl_df, "blast")
vs_tmp_dat <- truthfunction(vs_df, "vsearch")
st_tmp_dat <- truthfunction(st_df, "sintax")
nb_tmp_dat <- truthfunction(nb_df, "nbayes")
bo_tmp_dat <- truthfunction(bo_df, "bold")

## combine into single dataset for plotting:
plot_dat <- rbind(bl_tmp_dat, vs_tmp_dat, st_tmp_dat, nb_tmp_dat, bo_tmp_dat)

################################################################################
## plot
################################################################################

theme_devon <- function () { 
  theme_bw(base_size=12, base_family="Avenir") %+replace% 
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill="transparent", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA)
    )
}

## palette for plot
# red, orange, aqua, magenta, darkblue
pal5 <- c("#D55E00", '#E69F00', '#009E73', '#CC79A7', '#0072B2')

## set levels for plot 
plot_dat$TaxaLevel <- factor(plot_dat$TaxaLevel, 
                             levels = c("Class", "Order", "Family", "Genus", "Species"))

plot_dat$Classifier <- factor(plot_dat$Classifier,
                              levels = c("blast", "vsearch", "nbayes", "sintax", "bold"))

## and plot; save as 'cl_0_mockGroundTruth'; export at 800x400
ggplot(plot_dat, aes(x=TaxaLevel, y=Score, color=Classifier, group=Classifier)) +
  geom_point(size=1.5) +
  geom_line(alpha=0.5) +
  scale_color_manual(values = pal5) +
  theme_devon() +
  facet_grid(~ Metric) +
  theme(legend.position = "top") 


## export data as text files (one per Metric):
TAR_out <- plot_dat %>% 
  filter(Metric == "TAR") %>% 
  select(-Metric) %>% 
  mutate(Score = round(Score, 2)) %>% 
  spread(., key = TaxaLevel, value = Score)
write.table(TAR_out, 
            file = "~/Repos/tidybug/data/text_tables/classify_comps/mockTARdata.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE)

TDR_out <- plot_dat %>% 
  filter(Metric == "TDR") %>% 
  select(-Metric) %>% 
  mutate(Score = round(Score, 2)) %>% 
  spread(., key = TaxaLevel, value = Score)
write.table(TDR_out, 
            file = "~/Repos/tidybug/data/text_tables/classify_comps/mockTDRdata.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE)

F_out <- plot_dat %>% 
  filter(Metric == "F") %>% 
  select(-Metric) %>% 
  mutate(Score = round(Score, 2)) %>% 
  spread(., key = TaxaLevel, value = Score)
write.table(F_out, 
            file = "~/Repos/tidybug/data/text_tables/classify_comps/mockFdata.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE)
