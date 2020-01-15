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
df <- df %>% filter(SampleType != "mock")

## Hull Numbers function applied to calculate diversity measures: observed OTUs, simpson, and shannon
hillcurve.function <- function(data, filter_exp, filter_exp2) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  depth <- tmp.df %>% group_by(SeqID) %>% summarise(TotalReads = sum(Reads))
  Method <- tmp.df %>% distinct(Method)
  Filt <- tmp.df %>% distinct(Filt)
  tmp.mat <- dcast(tmp.df, SeqID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SeqID
  SeqID <- row.names(tmp.mat)
  tmp.mat$SeqID <- NULL
  tmp.hill <- renyi(tmp.mat, scales = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2), hill=TRUE)
  tmp <- data.frame(tmp.hill, SeqID, Method, Filt, row.names = NULL)
  tmp <- merge(tmp, depth)
}

## calculate Hill Numbers for unrarefied data
dada2.basic <- hillcurve.function(df, Method=="dada2", Filt=="basic")
dada2.standard <- hillcurve.function(df, Method=="dada2", Filt=="standard")
dada2.extra <- hillcurve.function(df, Method=="dada2", Filt=="extra")
deblur.basic <- hillcurve.function(df, Method=="deblur", Filt=="basic")
deblur.standard <- hillcurve.function(df, Method=="deblur", Filt=="standard")
deblur.extra <- hillcurve.function(df, Method=="deblur", Filt=="extra")
vsearch.basic <- hillcurve.function(df, Method=="vsearch", Filt=="basic")
vsearch.standard <- hillcurve.function(df, Method=="vsearch", Filt=="standard")
vsearch.extra <- hillcurve.function(df, Method=="vsearch", Filt=="extra")

rm(df)

## merge into single dataframe
all.df.hill <- rbind(dada2.basic, dada2.standard, dada2.extra, deblur.basic, deblur.standard, deblur.extra, vsearch.basic, vsearch.standard, vsearch.extra)
rm(dada2.basic, dada2.standard, dada2.extra, deblur.basic, deblur.standard, deblur.extra, vsearch.basic, vsearch.standard, vsearch.extra)

## gather q values into single column
colnames(all.df.hill)[2:10] <- c('0', '0.25', '0.5', '0.75', '1', '1.25', '1.5', '1.75', '2')
all.df.hill <- gather(all.df.hill, key="Hill_q", value="Hill_val", c('0', '0.25', '0.5', '0.75', '1', '1.25', '1.5', '1.75', '2'))

all.df.hill$logReads <- log(all.df.hill$TotalReads)

## reorder levels for plot
all.df.hill$Filt <- factor(all.df.hill$Filt, levels=c("basic", "standard", "extra"))

## plot; save as 11_HillCurve_guanodata_withReadDepths; export at 1000x1000
ggplot(all.df.hill, aes(x=Hill_q, y=Hill_val, color=logReads)) + 
  geom_jitter(data=all.df.hill %>% filter(TotalReads >= 5000), alpha=0.75) + 
  geom_jitter(data=all.df.hill %>% filter(TotalReads < 5000), alpha=0.5) + 
  scale_color_viridis_c(option = "plasma") +
  facet_grid(Filt ~ Method, scales = "free_y") +
  theme_devon() +
  theme(panel.background = element_rect(fill = "gray70")) +
  labs(title="", x="Hill number", y="Diversity", color="natural log \nReads per sample")

