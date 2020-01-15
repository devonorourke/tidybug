library(tidyverse)
library(reshape2)
library(vegan)
library(phyloseq)
library(FSA)

theme_devon <- function () { 
  theme_bw(base_size=12, base_family="Avenir") %+replace% 
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill="transparent", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA)
    )
}

## import data:
df <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/all.filtmethods.df.csv.gz")
df <- df %>% filter(SampleType == "mock")
df$Labeler <- paste(df$Method, df$Filt, df$Library, sep="-")

################################################################################
## data generation and stats behind rarefied mock data
################################################################################

## function to calculate Hill Numbers per Method + Filt (grouping all guano data among all libraries)
hill.function.wrare <- function(data, filter_exp, filter_exp2) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
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

rarefied.dada2.basic <- hill.function.wrare(df, Method=="dada2", Filt=="basic")
rarefied.dada2.standard <- hill.function.wrare(df, Method=="dada2", Filt=="standard")
rarefied.dada2.extra <- hill.function.wrare(df, Method=="dada2", Filt=="extra")
rarefied.deblur.basic <- hill.function.wrare(df, Method=="deblur", Filt=="basic")
rarefied.deblur.standard <- hill.function.wrare(df, Method=="deblur", Filt=="standard")
rarefied.deblur.extra <- hill.function.wrare(df, Method=="deblur", Filt=="extra")
rarefied.vsearch.basic <- hill.function.wrare(df, Method=="vsearch", Filt=="basic")
rarefied.vsearch.standard <- hill.function.wrare(df, Method=="vsearch", Filt=="standard")
rarefied.vsearch.extra <- hill.function.wrare(df, Method=="vsearch", Filt=="extra")

## merge into single dataframe
all.guano.hill.rare <- rbind(rarefied.dada2.basic, rarefied.dada2.standard, rarefied.dada2.extra, 
                             rarefied.deblur.basic, rarefied.deblur.standard, rarefied.deblur.extra, 
                             rarefied.vsearch.basic, rarefied.vsearch.standard, rarefied.vsearch.extra)

rm(list=ls(pattern = "rarefied*"))
all.guano.hill.rare$Method <- as.factor(all.guano.hill.rare$Method)

## split df by Hill q-value:
q0.hill.df.rare <- all.guano.hill.rare %>% filter(Hill_qType == "q=0")
q1.hill.df.rare <- all.guano.hill.rare %>% filter(Hill_qType == "q=1")
q2.hill.df.rare <- all.guano.hill.rare %>% filter(Hill_qType == "q=2")


######## run the ANOVAs on unrarefied data
q0.anova <- aov(Hill_value ~ Method * Filt, data=q0.hill.df.rare)
capture.output(summary(q0.anova), file = "~/Repos/tidybug/data/text_tables/anova_mock/mock_rare.anova.q0.txt")

q1.anova <- aov(Hill_value ~ Method * Filt, data=q1.hill.df.rare)
capture.output(summary(q1.anova), file = "~/Repos/tidybug/data/text_tables/anova_mock/mock_rare.anova.q1.txt")

q2.anova <- aov(Hill_value ~ Method * Filt, data=q2.hill.df.rare)
capture.output(summary(q2.anova), file = "~/Repos/tidybug/data/text_tables/anova_mock/mock_rare.anova.q2.txt")


## Kruskal-Wallis:
kw_df <- all.guano.hill.rare
kw_df$Labeler <- paste(kw_df$Method, kw_df$Filt, sep="-")
kw_df$Labeler <- as.factor(kw_df$Labeler)
kw_df_q0 <- kw_df %>% filter(Hill_qType == "q=0")
kw_df_q1 <- kw_df %>% filter(Hill_qType == "q=1")
kw_df_q2 <- kw_df %>% filter(Hill_qType == "q=2")

capture.output(kruskal.test(Hill_value ~ Labeler, data = kw_df_q0),
               file = "~/Repos/tidybug/data/text_tables/kw_mock/mock_rare_kw.q0.txt")
capture.output(kruskal.test(Hill_value ~ Labeler, data = kw_df_q1),
               file = "~/Repos/tidybug/data/text_tables/kw_mock/mock_rare_kw.q1.txt")
capture.output(kruskal.test(Hill_value ~ Labeler, data = kw_df_q2),
               file = "~/Repos/tidybug/data/text_tables/kw_mock/mock_rare_kw.q2.txt")


##run Dunn's test for pairwise differences in diversity estimates
## add Labeler to data frame
dunn_df <- all.guano.hill.rare
dunn_df$Labeler <- paste(dunn_df$Method, dunn_df$Filt, sep="-")
dunn_df$Labeler <- as.factor(dunn_df$Labeler)
dunn_df_q0 <- filter(dunn_df, Hill_qType=="q=0")
dunn_df_q1 <- filter(dunn_df, Hill_qType=="q=1")
dunn_df_q2 <- filter(dunn_df, Hill_qType=="q=2")

dunnout_q0 <- dunnTest(Hill_value ~ Labeler, data=dunn_df_q0, method="bh")
dunnout_q0 <- dunnout_q0$res %>% arrange(P.adj)
write.csv(dunnout_q0, file="~/Repos/tidybug/data/text_tables/dunn_mock/mock_rare_dunn_q0.csv", quote = FALSE, row.names = FALSE)

dunnout_q1 <- dunnTest(Hill_value ~ Labeler, data=dunn_df_q1, method="bh")
dunnout_q1 <- dunnout_q1$res %>% arrange(P.adj)
write.csv(dunnout_q1, file="~/Repos/tidybug/data/text_tables/dunn_mock/mock_rare_dunn_q1.csv", quote = FALSE, row.names = FALSE)

dunnout_q2 <- dunnTest(Hill_value ~ Labeler, data=dunn_df_q2, method="bh")
dunnout_q2 <- dunnout_q2$res %>% arrange(P.adj)
write.csv(dunnout_q2, file="~/Repos/tidybug/data/text_tables/dunn_mock/mock_rare_dunn_q2.csv", quote = FALSE, row.names = FALSE)


################################################################################
## data generation and stats behind rarefied mock data
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


## split df by Hill q-value:
q0.hill.df.norare <- all.guano.hill.norare %>% filter(Hill_qType == "q=0")
q1.hill.df.norare <- all.guano.hill.norare %>% filter(Hill_qType == "q=1")
q2.hill.df.norare <- all.guano.hill.norare %>% filter(Hill_qType == "q=2")


######## run the ANOVAs on unrarefied data
q0.anova_norare <- aov(Hill_value ~ Method * Filt, data=q0.hill.df.norare)
capture.output(summary(q0.anova_norare), file = "~/Repos/tidybug/data/text_tables/anova_mock/mock_norare.anova.q0.txt")

q1.anova_norare <- aov(Hill_value ~ Method * Filt, data=q1.hill.df.norare)
capture.output(summary(q1.anova_norare), file = "~/Repos/tidybug/data/text_tables/anova_mock/mock_norare.anova.q1.txt")

q2.anova_norare <- aov(Hill_value ~ Method * Filt, data=q2.hill.df.norare)
capture.output(summary(q2.anova_norare), file = "~/Repos/tidybug/data/text_tables/anova_mock/mock_norare.anova.q2.txt")

## Kruskal-Wallis:
kw_df_norare <- all.guano.hill.norare
kw_df_norare$Labeler <- paste(kw_df_norare$Method, kw_df_norare$Filt, sep="-")
kw_df_norare$Labeler <- as.factor(kw_df_norare$Labeler)
kw_df_norare_q0 <- kw_df_norare %>% filter(Hill_qType == "q=0")
kw_df_norare_q1 <- kw_df_norare %>% filter(Hill_qType == "q=1")
kw_df_norare_q2 <- kw_df_norare %>% filter(Hill_qType == "q=2")

capture.output(kruskal.test(Hill_value ~ Labeler, data = kw_df_norare_q0), 
               file="~/Repos/tidybug/data/text_tables/kw_mock/mock_norare_kw_q0.txt", quote = FALSE, row.names = FALSE)
capture.output(kruskal.test(Hill_value ~ Labeler, data = kw_df_norare_q1),
               file="~/Repos/tidybug/data/text_tables/kw_mock/mock_norare_kw_q1.txt", quote = FALSE, row.names = FALSE)
capture.output(kruskal.test(Hill_value ~ Labeler, data = kw_df_norare_q2),
               file="~/Repos/tidybug/data/text_tables/kw_mock/mock_norare_kw_q2.txt", quote = FALSE, row.names = FALSE)

##run Dunn's test for pairwise differences in diversity estimates
## add Labeler to data frame
dunn_df_norare <- all.guano.hill.rare
dunn_df_norare$Labeler <- paste(dunn_df_norare$Method, dunn_df_norare$Filt, sep="-")
dunn_df_norare$Labeler <- as.factor(dunn_df_norare$Labeler)
dunn_df_norare_q0 <- filter(dunn_df_norare, Hill_qType=="q=0")
dunn_df_norare_q1 <- filter(dunn_df_norare, Hill_qType=="q=1")
dunn_df_norare_q2 <- filter(dunn_df_norare, Hill_qType=="q=2")

dunnout_norare_q0 <- dunnTest(Hill_value ~ Labeler, data=dunn_df_norare_q0, method="bh")
dunnout_norare_q0 <- dunnout_norare_q0$res %>% arrange(P.adj)
write.csv(dunnout_norare_q0, file="~/Repos/tidybug/data/text_tables/dunn_mock/mock_norare_dunn_q0.csv", quote = FALSE, row.names = FALSE)

dunnout_norare_q1 <- dunnTest(Hill_value ~ Labeler, data=dunn_df_norare_q1, method="bh")
dunnout_norare_q1 <- dunnout_norare_q1$res %>% arrange(P.adj)
write.csv(dunnout_norare_q1, file="~/Repos/tidybug/data/text_tables/dunn_mock/mock_norare_dunn_q1.csv", quote = FALSE, row.names = FALSE)

dunnout_norare_q2 <- dunnTest(Hill_value ~ Labeler, data=dunn_df_norare_q2, method="bh")
dunnout_norare_q2 <- dunnout_norare_q2$res %>% arrange(P.adj)
write.csv(dunnout_norare_q2, file="~/Repos/tidybug/data/text_tables/dunn_mock/mock_norare_dunn_q2.csv", quote = FALSE, row.names = FALSE)
