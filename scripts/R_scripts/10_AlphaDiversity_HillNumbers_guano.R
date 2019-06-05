library(tidyverse)
library(reshape2)
library(vegan)
library(phyloseq)

theme_devon <- function () { 
  theme_bw(base_size=12, base_family="Avenir") %+replace% 
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill="transparent", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA)
    )
}

## import data and select mock samples:
df <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/all.filtmethods.df.csv.gz")
## filter out mock samples:
df <- df %>% filter(SampleType != "mock") %>% select(-StudyID, -Alias)
df$Labeler <- paste(df$Method, df$Filt, df$Library, sep="-")

## function to calculate Hill Numbers per Method + Filt (grouping all guano data among all libraries)
hill.function.wrare <- function(data, filter_exp, filter_exp2, filter_exp3) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  filter_exp_enq3 <- enquo(filter_exp3)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2) %>% filter(!!filter_exp_enq3)
  Method <- tmp.df %>% distinct(Method)
  Filt <- tmp.df %>% distinct(Filt)
  tmp.mat <- dcast(tmp.df, HashID ~ SeqID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$HashID
  tmp.mat$HashID <- NULL
  tmp_otutbl <- otu_table(tmp.mat, taxa_are_rows = TRUE)
  tmp_rphy <- rarefy_even_depth(tmp_otutbl, sample.size = 5000, replace = FALSE, rngseed = 123)
  tmp_raremat = as(otu_table(tmp_rphy), "matrix")
  tmp.hill <- renyi(t(tmp_raremat), scales = c(0,1,2), hill=TRUE)
  tmp.hill <- data.frame(tmp.hill, Method, Filt, row.names = NULL)
  colnames(tmp.hill)[1:3] <- c("q=0", "q=1", "q=2")
  tmp_out <- gather(tmp.hill, key="Hill_qType", value = "Hill_value", c("q=0", "q=1", "q=2"))
}

rarefied.dada2.basic <- hill.function.wrare(df, Method=="dada2", Filt=="basic", SampleType=="sample")
rarefied.dada2.standard <- hill.function.wrare(df, Method=="dada2", Filt=="standard", SampleType=="sample")
rarefied.dada2.extra <- hill.function.wrare(df, Method=="dada2", Filt=="extra", SampleType=="sample")
rarefied.deblur.basic <- hill.function.wrare(df, Method=="deblur", Filt=="basic", SampleType=="sample")
rarefied.deblur.standard <- hill.function.wrare(df, Method=="deblur", Filt=="standard", SampleType=="sample")
rarefied.deblur.extra <- hill.function.wrare(df, Method=="deblur", Filt=="extra", SampleType=="sample")
rarefied.vsearch.basic <- hill.function.wrare(df, Method=="vsearch", Filt=="basic", SampleType=="sample")
rarefied.vsearch.standard <- hill.function.wrare(df, Method=="vsearch", Filt=="standard", SampleType=="sample")
rarefied.vsearch.extra <- hill.function.wrare(df, Method=="vsearch", Filt=="extra", SampleType=="sample")

## merge into single dataframe
all.guano.hill.rare <- rbind(rarefied.dada2.basic, rarefied.dada2.standard, rarefied.dada2.extra, 
                             rarefied.deblur.basic, rarefied.deblur.standard, rarefied.deblur.extra, 
                             rarefied.vsearch.basic, rarefied.vsearch.standard, rarefied.vsearch.extra)

rm(list=ls(pattern = "rarefied*"))

## set the levels
all.guano.hill.rare$Filt <- factor(all.guano.hill.rare$Filt, levels = c("basic", "standard", "extra"))

## palette:
pal3 <- c('#9f9244', '#6c42b8', '#628a47')

## plot; save as 5_figure_guano_HillDiversities_perPipeline ; export at 750x750
ggplot(data=all.guano.hill.rare, aes(x=Method, y=Hill_value, fill=Method)) +
  geom_jitter(alpha = 0.05) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = pal3) +
  facet_grid(Filt ~ Hill_qType) +
  scale_y_continuous(trans = "log2") +
  labs(x="", y="sequence variant equivalents", color="", shape="") +
  theme_devon() + 
  theme(legend.position = "none",
        strip.text = element_text(size = 12))


################################################################################
## same process as above, but we do not rarefy...
################################################################################

## function to calculate Hill Numbers per Method + Filt (grouping all guano data among all libraries)
hill.function.norare <- function(data, filter_exp, filter_exp2, filter_exp3) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  filter_exp_enq3 <- enquo(filter_exp3)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2) %>% filter(!!filter_exp_enq3)
  Method <- tmp.df %>% distinct(Method)
  Filt <- tmp.df %>% distinct(Filt)
  tmp.mat <- dcast(tmp.df, HashID ~ SeqID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$HashID
  tmp.mat$HashID <- NULL
  tmp.hill <- renyi(t(tmp.mat), scales = c(0,1,2), hill=TRUE)
  tmp.hill <- data.frame(tmp.hill, Method, Filt, row.names = NULL)
  colnames(tmp.hill)[1:3] <- c("q=0", "q=1", "q=2")
  tmp_out <- gather(tmp.hill, key="Hill_qType", value = "Hill_value", c("q=0", "q=1", "q=2"))
}

unrarefied.dada2.basic <- hill.function.norare(df, Method=="dada2", Filt=="basic", SampleType=="sample")
unrarefied.dada2.standard <- hill.function.norare(df, Method=="dada2", Filt=="standard", SampleType=="sample")
unrarefied.dada2.extra <- hill.function.norare(df, Method=="dada2", Filt=="extra", SampleType=="sample")
unrarefied.deblur.basic <- hill.function.norare(df, Method=="deblur", Filt=="basic", SampleType=="sample")
unrarefied.deblur.standard <- hill.function.norare(df, Method=="deblur", Filt=="standard", SampleType=="sample")
unrarefied.deblur.extra <- hill.function.norare(df, Method=="deblur", Filt=="extra", SampleType=="sample")
unrarefied.vsearch.basic <- hill.function.norare(df, Method=="vsearch", Filt=="basic", SampleType=="sample")
unrarefied.vsearch.standard <- hill.function.norare(df, Method=="vsearch", Filt=="standard", SampleType=="sample")
unrarefied.vsearch.extra <- hill.function.norare(df, Method=="vsearch", Filt=="extra", SampleType=="sample")

## merge into single dataframe
all.guano.hill.norare <- rbind(unrarefied.dada2.basic, unrarefied.dada2.standard, unrarefied.dada2.extra, 
                             unrarefied.deblur.basic, unrarefied.deblur.standard, unrarefied.deblur.extra, 
                             unrarefied.vsearch.basic, unrarefied.vsearch.standard, unrarefied.vsearch.extra)

rm(list=ls(pattern = "unrarefied*"))


## set the levels
all.guano.hill.norare$Filt <- factor(all.guano.hill.norare$Filt, levels = c("basic", "standard", "extra"))

## plot; save as 5_figure_guano_HillDiversities_perPipeline_noRare ; export at 750x750
ggplot(data=all.guano.hill.norare, aes(x=Method, y=Hill_value, fill=Method)) +
  geom_jitter(alpha = 0.05) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = pal3) +
  facet_grid(Filt ~ Hill_qType) +
  scale_y_continuous(trans = "log2") +
  labs(x="", y="sequence variant equivalents", color="", shape="") +
  theme_devon() + 
  theme(legend.position = "none",
        strip.text = element_text(size = 12))
