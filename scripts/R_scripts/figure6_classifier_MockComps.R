library(tidyverse)
library(reshape2)

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
  tmp <- as.data.frame(apply(tmp, 2, function(y) (gsub("Unassigned", NA, y))))
  tmp[] <- lapply(tmp, as.character)
  tmp %>% gather(key = level, value = Taxon, -IMalias)
}

exp_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/expected_mockData/mock_expected_taxa.txt"

bo95_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/bold/boldAPI_mock_lca95.txt"
bo97_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/bold/boldAPI_mock_lca97.txt"
bo99_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/bold/boldAPI_mock_lca99.txt"

nb30_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/nbayes/nbayes_conf30.txt"
nb50_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/nbayes/nbayes_conf50.txt"
nb70_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/nbayes/nbayes_conf70.txt"
nb80_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/nbayes/nbayes_conf80.txt"
nb90_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/nbayes/nbayes_conf90.txt"

sn30_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/sintax/sin_c30.txt"
sn50_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/sintax/sin_c50.txt"
sn70_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/sintax/sin_c70.txt"
sn80_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/sintax/sin_c80.txt"
sn90_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/sintax/sin_c90.txt"

bl95_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/bold/bl95.txt"
bl97_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/bold/bl97.txt"
bl99_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/bold/bl99.txt"

vs95_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/vsearch/vs95.txt"
vs97_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/vsearch/vs97.txt"
vs99_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/vsearch/vs99.txt"
vs95top_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/vsearch/vs95top.txt"
vs97top_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/vsearch/vs97top.txt"
vs99top_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/vsearch/vs99top.txt"

exp_df <- readerfunction(exp_url)

nb30_df <- readerfunction(nb30_url)
nb50_df <- readerfunction(nb50_url)
nb70_df <- readerfunction(nb70_url)
nb80_df <- readerfunction(nb80_url)
nb90_df <- readerfunction(nb90_url)

sn30_df <- readerfunction(sn30_url)
sn50_df <- readerfunction(sn50_url)
sn70_df <- readerfunction(sn70_url)
sn80_df <- readerfunction(sn80_url)
sn90_df <- readerfunction(sn90_url)

bl95_df <- readerfunction(bl95_url)
bl97_df <- readerfunction(bl97_url)
bl99_df <- readerfunction(bl99_url)

bo95_df <- readerfunction(bo95_url)
bo97_df <- readerfunction(bo97_url)
bo99_df <- readerfunction(bo99_url)

vs95_df <- readerfunction(vs95_url)
vs97_df <- readerfunction(vs97_url)
vs99_df <- readerfunction(vs99_url)
vs95top_df <- readerfunction(vs95top_url)
vs97top_df <- readerfunction(vs97top_url)
vs99top_df <- readerfunction(vs99top_url)

truthfunction <- function(obs_data, classifier, parameter) {
  selectTests <- c("TAR", "TDR", "F")
  tmp1 <- inner_join(exp_df, obs_data, by = c("IMalias", "level"), suffix = c("_exp", "_obs")) %>% 
    mutate(truth = case_when(
      Taxon_exp == Taxon_obs  ~ "TP",
      is.na(Taxon_exp) == is.na(Taxon_obs)  ~ "TN",
      is.na(Taxon_exp) & !is.na(Taxon_obs)  ~ "FP",
      !is.na(Taxon_exp) & is.na(Taxon_obs)  ~ "FN")) %>% 
    select(-starts_with("Taxon")) %>%
    group_by(level, truth) %>% 
    summarize(counts = n()) %>% 
    spread(truth, counts, fill = 0)
  if (is.null(tmp1$FP)) {
    tmp2 <- data.frame(tmp1, FP = rep(0, 5))
  } else {
    tmp2 <- tmp1
  }
  tmp2 %>% 
    group_by(level) %>% 
    mutate(TAR = TP / (TP + FP)) %>% 
    mutate(TDR = TP / (TP + FN)) %>% 
    mutate(F = 2 * ((TAR*TDR)/(TAR+TDR))) %>% 
    melt(., id.vars = 'level', value.name = 'Score', variable.name = 'Metric') %>%
    filter(Metric %in% selectTests) %>% 
    mutate(Classifier = classifier) %>% 
    mutate(Parameter = parameter)
}

bo95_tmp_dat <- truthfunction(bo95_df, "boldAPI+LCA", "95")
bo97_tmp_dat <- truthfunction(bo97_df, "boldAPI+LCA", "97")
bo99_tmp_dat <- truthfunction(bo99_df, "boldAPI+LCA", "99")

sn30_tmp_dat <- truthfunction(sn30_df, "SINTAX (vsearch)", "30")
sn50_tmp_dat <- truthfunction(sn50_df, "SINTAX (vsearch)", "50")
sn70_tmp_dat <- truthfunction(sn70_df, "SINTAX (vsearch)", "70")
sn80_tmp_dat <- truthfunction(sn80_df, "SINTAX (vsearch)", "80")
sn90_tmp_dat <- truthfunction(sn90_df, "SINTAX (vsearch)", "90")

nb30_tmp_dat <- truthfunction(nb30_df, "NaiveBayes (q2)", "30")
nb50_tmp_dat <- truthfunction(nb50_df, "NaiveBayes (q2)", "50")
nb70_tmp_dat <- truthfunction(nb70_df, "NaiveBayes (q2)", "70")
nb80_tmp_dat <- truthfunction(nb80_df, "NaiveBayes (q2)", "80")
nb90_tmp_dat <- truthfunction(nb90_df, "NaiveBayes (q2)", "90")

bl95_tmp_dat <- truthfunction(bl95_df, "BLAST+LCA (q2)", "95")
bl97_tmp_dat <- truthfunction(bl97_df, "BLAST+LCA (q2)", "97")
bl99_tmp_dat <- truthfunction(bl99_df, "BLAST+LCA (q2)", "99")

vs95_tmp_dat <- truthfunction(vs95_df, "VSEARCH+LCA (q2)", "95")
vs97_tmp_dat <- truthfunction(vs97_df, "VSEARCH+LCA (q2)", "97")
vs99_tmp_dat <- truthfunction(vs99_df, "VSEARCH+LCA (q2)", "99")
vs95top_tmp_dat <- truthfunction(vs95top_df, "VSEARCH+LCA+top_hit (q2)", "95")
vs97top_tmp_dat <- truthfunction(vs97top_df, "VSEARCH+LCA+top_hit (q2)", "97")
vs99top_tmp_dat <- truthfunction(vs99top_df, "VSEARCH+LCA+top_hit (q2)", "99")


plot_dat <- rbind(bo95_tmp_dat, bo97_tmp_dat, bo99_tmp_dat,
                  sn30_tmp_dat, sn50_tmp_dat, sn70_tmp_dat, sn80_tmp_dat, sn90_tmp_dat,
                  nb30_tmp_dat, nb50_tmp_dat, nb70_tmp_dat, nb80_tmp_dat, nb90_tmp_dat,
                  bl95_tmp_dat, bl97_tmp_dat, bl99_tmp_dat,
                  vs95_tmp_dat, vs97_tmp_dat, vs99_tmp_dat,
                  vs95top_tmp_dat, vs97top_tmp_dat, vs99top_tmp_dat)


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
#pal5 <- c("#D55E00", '#E69F00', 'bur' '#009E73', '#CC79A7', '#0072B2')
pal6 <- c("#D55E00", 'orange', 'darkorange', '#009E73', '#CC79A7', '#0072B2')

## set levels for plot 
plot_dat$level <- factor(plot_dat$level, 
                         levels = c("Class", "Order", "Family", "Genus", "Species"))

plot_dat$Classifier <- factor(plot_dat$Classifier,
                              levels = c("BLAST+LCA (q2)", "VSEARCH+LCA (q2)", "VSEARCH+LCA+top_hit (q2)","NaiveBayes (q2)", "SINTAX (vsearch)", "boldAPI+LCA"))

## and plot; save as 'figure6_classifierComps_wVals.png'; will need to modify in Illustrator to fix internal integers that overlap
pdodge_val = 0.75
ymin_val = 0.3

ggplot(plot_dat %>% filter(Metric == "TDR"), 
       aes(x=level, y=Score, fill=Classifier, group=Classifier, label=Parameter)) +
  geom_linerange(aes(ymin=ymin_val, ymax=Score),
                 color="gray50",
                 position = position_dodge(width = pdodge_val + 0.5),
                 size = 1.5) +
  geom_rect(aes(xmin = 1, ymin = Score - 0.03, xmax = 2, ymax = Score + 0.03, colour=Classifier),
            position = position_dodge(pdodge_val + 0.5), fill="white") +
  facet_wrap(~ level, scales = "free_x") +
  scale_color_manual(values = pal6) +
  geom_text(size = 5, position = position_dodge(pdodge_val + 0.5), color="black") +
  theme_devon() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size=11),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        strip.text = element_text(size=12),
        legend.position="top", legend.text = element_text(size=13),legend.title = element_text(size=12)) +
  labs(y = "TDR score\n", x="") +
  scale_y_continuous(limits = c(ymin_val,1.1))
