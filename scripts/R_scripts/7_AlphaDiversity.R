# Plots used in paper

library(tidyverse)
library(reshape2)
library(vegan)


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
df <- df %>% filter(SampleType != "mock")

## alpha function applied to calculate diversity measures: observed OTUs, simpson, and shannon
alpha.function <- function(data, filter_exp, filter_exp2) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  tmp1 <- tmp.df %>% distinct(Method)
  tmp2 <- tmp.df %>% distinct(Filt)
  Label <- paste(tmp1, tmp2, sep="-")
  tmp.mat <- dcast(tmp.df, SeqID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SeqID
  tmp.mat$SeqID <- NULL
  SeqID <- row.names(tmp.mat)
  Simpson <- diversity(tmp.mat, index="simpson")
  Shannon <- diversity(tmp.mat, index="shannon")
  OTUs <- specnumber(tmp.mat)
  data.frame(SeqID, OTUs, Simpson, Shannon, row.names = NULL)
}

## calculate alphas
dada2.basic <- alpha.function(df, Method=="dada2", Filt=="basic")
dada2.standard <- alpha.function(df, Method=="dada2", Filt=="standard")
dada2.extra <- alpha.function(df, Method=="dada2", Filt=="extra")
deblur.basic <- alpha.function(df, Method=="deblur", Filt=="basic")
deblur.standard <- alpha.function(df, Method=="deblur", Filt=="standard")
deblur.extra <- alpha.function(df, Method=="deblur", Filt=="extra")
vsearch.basic <- alpha.function(df, Method=="vsearch", Filt=="basic")
vsearch.standard <- alpha.function(df, Method=="vsearch", Filt=="standard")
vsearch.extra <- alpha.function(df, Method=="vsearch", Filt=="extra")

## merge into single dataframe
all <- rbind(dada2.basic, dada2.standard, dada2.extra, deblur.basic, deblur.standard, deblur.extra, vsearch.basic, vsearch.standard, vsearch.extra)
rm(dada2.basic, dada2.standard, dada2.extra, deblur.basic, deblur.standard, deblur.extra, vsearch.basic, vsearch.standard, vsearch.extra)
rm(df)

## split the $Label
all.simpson <- separate(all.simpson, col=Label, into=c("Method", "Filt", "AlphaTest"), sep = "-")

## set the levels
all.simpson$Filt <- factor(all.simpson$Filt,levels = c("basic", "standard", "extra"))

## set the palette to match figure 6
pipepal3 <- c('#9f9244', '#6c42b8', '#628a47', '#a9d190')

## the full plot is really noisy if you plot by Method
ggplot(data = all.simpson, aes(x = Method, y = AlphaVal, color=Method)) +
  #geom_hline(yintercept = 50, linetype="dotted", color="firebrick", size=1) +
  geom_jitter(alpha=0.55, width = 0.25) +
  facet_grid(Filt ~ .) +
  scale_color_manual(values=pipepal3) +
  labs(title="", x="", y="Simpson 1-D", color="") +
  theme_devon() +
  theme(legend.position="top", legend.text = element_text(size = 12), axis.text.x = element_text(angle=22.5, hjust = 1))



## color scheme below:
#Vse.arch: dark=#9f9244; light=#ebdb8e
#  Da.da2: dark=#6c42b8 ; light=#c8b2e8
#  De.blur: dark=#628a47; light=#a9d190

######## Stat test with Kruskal-Wallis
tmp <- all.simpson
tmp$Label <- paste(all.simpson$Method, all.simpson$Filt, sep="-")
tmp$Label <- as.factor(tmp$Label)
kruskal.test(AlphaVal ~ Label, data=tmp)


## Try plotting by the $Label, and showing boxplot
ggplot(data = tmp, aes(x = Label, y = AlphaVal, color=Method)) +
  geom_boxplot() +
  geom_jitter(alpha=0.15, width = 0.25) +
  #facet_grid(Filt ~ .) +
  scale_color_manual(values=pipepal3) +
  labs(title="", x="", y="Simpson 1-D", color="") +
  theme_devon() +
  theme(legend.position="top", legend.text = element_text(size = 12), axis.text.x = element_text(angle=22.5, hjust = 1))



### rerun, but with alpha function rarefying data first!