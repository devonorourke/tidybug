library(tidyverse)
library(reshape2)
library(vegan)
library(cowplot)

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

## Hull Numbers function applied to calculate diversity measures: observed OTUs, simpson, and shannon
hillnumber.function <- function(data, filter_exp, filter_exp2) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  Method <- tmp.df %>% distinct(Method)
  Filt <- tmp.df %>% distinct(Filt)
  tmp.mat <- dcast(tmp.df, SeqID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SeqID
  SeqID <- row.names(tmp.mat)
  tmp.mat$SeqID <- NULL
  tmp.hill <- renyi(tmp.mat, scales = c(0,1,2), hill=TRUE)
  tmp <- data.frame(tmp.hill, SeqID, Method, Filt, row.names = NULL)
  colnames(tmp)[1:3] <- c("q=0", "q=1", "q=2")
  gather(tmp, key="Hill_qType", value = "Hill_value", c('q=0', 'q=1', 'q=2'))
}

## calculate Hill Numbers for unrarefied data
dada2.basic <- hillnumber.function(df, Method=="dada2", Filt=="basic")
dada2.standard <- hillnumber.function(df, Method=="dada2", Filt=="standard")
dada2.extra <- hillnumber.function(df, Method=="dada2", Filt=="extra")
deblur.basic <- hillnumber.function(df, Method=="deblur", Filt=="basic")
deblur.standard <- hillnumber.function(df, Method=="deblur", Filt=="standard")
deblur.extra <- hillnumber.function(df, Method=="deblur", Filt=="extra")
vsearch.basic <- hillnumber.function(df, Method=="vsearch", Filt=="basic")
vsearch.standard <- hillnumber.function(df, Method=="vsearch", Filt=="standard")
vsearch.extra <- hillnumber.function(df, Method=="vsearch", Filt=="extra")

## merge into single dataframe
all.df.hill <- rbind(dada2.basic, dada2.standard, dada2.extra, deblur.basic, deblur.standard, deblur.extra, vsearch.basic, vsearch.standard, vsearch.extra)
rm(dada2.basic, dada2.standard, dada2.extra, deblur.basic, deblur.standard, deblur.extra, vsearch.basic, vsearch.standard, vsearch.extra)

## set the levels
all.df.hill$Filt <- factor(all.df.hill$Filt, levels = c("basic", "standard", "extra"))
## add labeler to make facet labels easier to understand
all.df.hill$Labeler <- paste(all.df.hill$Method, all.df.hill$Hill_qType, sep="  ")

## split df by Hill q-value:
q0.hill.df <- all.df.hill %>% filter(Hill_qType == "q=0")
q1.hill.df <- all.df.hill %>% filter(Hill_qType == "q=1")
q2.hill.df <- all.df.hill %>% filter(Hill_qType == "q=2")


######## run the ANOVAs on unrarefied data
q0.anova <- aov(Hill_value ~ Method * Filt, data=q0.hill.df)
summary(q0.anova)

q1.anova <- aov(Hill_value ~ Method * Filt, data=q1.hill.df)
summary(q1.anova)

q2.anova <- aov(Hill_value ~ Method * Filt, data=q2.hill.df)
summary(q2.anova)

## create similar dataframe with rarefied data:
rrarewithdrop <- function(x, sample) {
  rrarefy(x[rowSums(x) >= sample, , drop = FALSE], sample)
}

rare.function <- function(data, filter_exp, filter_exp2) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  meta.df <- tmp.df %>% distinct(SeqID, Library, Method, Filt, SampleType)
  tmp.mat <- dcast(tmp.df, SeqID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SeqID
  SeqID <- row.names(tmp.mat)
  tmp.mat$SeqID <- NULL
  tmp.mat2 <- data.frame(rrarewithdrop(tmp.mat, 5000))
  tmp.mat2$SeqID <- row.names(tmp.mat2)
  tmp.df2 <- melt(data = tmp.mat2, id.vars = "SeqID",
                  variable.name = "HashID", 
                  value.name = "Reads") %>% 
    filter(Reads > 0) %>%
    mutate(HashID=sub('.', '', HashID))
  tmp.df2 <- merge(tmp.df2, meta.df)
}

tmp.dada2.basic <- rare.function(df, Method=="dada2", Filt=="basic")
tmp.dada2.standard <- rare.function(df, Method=="dada2", Filt=="standard")
tmp.dada2.extra <- rare.function(df, Method=="dada2", Filt=="extra")
tmp.deblur.basic <- rare.function(df, Method=="deblur", Filt=="basic")
tmp.deblur.standard <- rare.function(df, Method=="deblur", Filt=="standard")
tmp.deblur.extra <- rare.function(df, Method=="deblur", Filt=="extra")
tmp.vsearch.basic <- rare.function(df, Method=="vsearch", Filt=="basic")
tmp.vsearch.standard <- rare.function(df, Method=="vsearch", Filt=="standard")
tmp.vsearch.extra <- rare.function(df, Method=="vsearch", Filt=="extra")
rare.df <- rbind(tmp.dada2.basic, tmp.dada2.standard, tmp.dada2.extra,
                 tmp.deblur.basic, tmp.deblur.standard, tmp.deblur.extra,
                 tmp.vsearch.basic, tmp.vsearch.standard, tmp.vsearch.extra)
rm(list=ls(pattern="tmp.*"))

## calculate Hill Numbers for rarefied data
r.dada2.basic <- hillnumber.function(rare.df, Method=="dada2", Filt=="basic")
r.dada2.standard <- hillnumber.function(rare.df, Method=="dada2", Filt=="standard")
r.dada2.extra <- hillnumber.function(rare.df, Method=="dada2", Filt=="extra")
r.deblur.basic <- hillnumber.function(rare.df, Method=="deblur", Filt=="basic")
r.deblur.standard <- hillnumber.function(rare.df, Method=="deblur", Filt=="standard")
r.deblur.extra <- hillnumber.function(rare.df, Method=="deblur", Filt=="extra")
r.vsearch.basic <- hillnumber.function(rare.df, Method=="vsearch", Filt=="basic")
r.vsearch.standard <- hillnumber.function(rare.df, Method=="vsearch", Filt=="standard")
r.vsearch.extra <- hillnumber.function(rare.df, Method=="vsearch", Filt=="extra")

## merge into single dataframe
all.r.df.hill <- rbind(r.dada2.basic, r.dada2.standard, r.dada2.extra, r.deblur.basic, r.deblur.standard, r.deblur.extra, r.vsearch.basic, r.vsearch.standard, r.vsearch.extra)
rm(r.dada2.basic, r.dada2.standard, r.dada2.extra, r.deblur.basic, r.deblur.standard, r.deblur.extra, r.vsearch.basic, r.vsearch.standard, r.vsearch.extra)

## split df by Hill q-value:
r.q0.hill.df <- all.r.df.hill %>% filter(Hill_qType == "q=0")
r.q1.hill.df <- all.r.df.hill %>% filter(Hill_qType == "q=1")
r.q2.hill.df <- all.r.df.hill %>% filter(Hill_qType == "q=2")

### run the ANOVAs on rarefied data
### going to rarefy data to 5000 reads
r.q0.anova <- aov(Hill_value ~ Method * Filt, data=r.q0.hill.df)
summary(r.q0.anova)

r.q1.anova <- aov(Hill_value ~ Method * Filt, data=r.q1.hill.df)
summary(r.q1.anova)

r.q2.anova <- aov(Hill_value ~ Method * Filt, data=r.q2.hill.df)
summary(r.q2.anova)


## Kruskal-Wallis:
kw_df <- all.df.hill
kw_df$Labeler <- paste(kw_df$Method, kw_df$Filt, sep="-")
kw_df$Labeler <- as.factor(kw_df$Labeler)
kw_df_q0 <- kw_df %>% filter(Hill_qType == "q=0")
kw_df_q1 <- kw_df %>% filter(Hill_qType == "q=1")
kw_df_q2 <- kw_df %>% filter(Hill_qType == "q=2")

kruskal.test(Hill_value ~ Labeler, data = kw_df_q0) ## significant difference: chi-squared = 32.381, df = 8, p-value = 7.959e-05
kruskal.test(Hill_value ~ Labeler, data = kw_df_q1) ## significant difference: chi-squared = 33.29, df = 8, p-value = 5.46e-05
kruskal.test(Hill_value ~ Labeler, data = kw_df_q2) ## significant difference: chi-squared = 35, df = 8, p-value = 2.674e-05


##run Dunn's test for pairwise differences in diversity estimates
## add Labeler to data frame
dunn_df <- all.df.hill
dunn_df$Labeler <- paste(dunn_df$Method, dunn_df$Filt, sep="-")
dunn_df_q0 <- filter(dunn_df, Hill_qType=="q=0")
dunn_df_q1 <- filter(dunn_df, Hill_qType=="q=1")
dunn_df_q2 <- filter(dunn_df, Hill_qType=="q=2")

dunnout_q0 <- dunnTest(Hill_value ~ Labeler, data=dunn_df_q0, method="bh")
dunnout_q0 <- dunnout_q0$res %>% arrange(P.adj)
write.csv(dunnout_q0, file="~/Repos/tidybug/data/text_tables/mock_dunn_q0.csv", quote = FALSE, row.names = FALSE)

dunnout_q1 <- dunnTest(Hill_value ~ Labeler, data=dunn_df_q1, method="bh")
dunnout_q1 <- dunnout_q1$res %>% arrange(P.adj)
write.csv(dunnout_q1, file="~/Repos/tidybug/data/text_tables/mock_dunn_q1.csv", quote = FALSE, row.names = FALSE)

dunnout_q2 <- dunnTest(Hill_value ~ Labeler, data=dunn_df_q2, method="bh")
dunnout_q2 <- dunnout_q2$res %>% arrange(P.adj)
write.csv(dunnout_q2, file="~/Repos/tidybug/data/text_tables/mock_dunn_q2.csv", quote = FALSE, row.names = FALSE)
