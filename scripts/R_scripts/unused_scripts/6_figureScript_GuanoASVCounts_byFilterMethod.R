# Plots used in paper

library(tidyverse)
library(scales)
library(vegan)
library(reshape2)

## plotting setup:
## generate 3 color palette to distinguish between filtering pipelines:
pal3 <- c('#9f9244', '#6c42b8', '#628a47', '#a9d190')
# theme function for custom plot style:
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

## create similar dataframe with rarefied data:
rrarewithdrop <- 
  function(x, sample) 
  {
    rrarefy(x[rowSums(x) >= sample, , drop = FALSE], sample)
  }

rare.function <- function(data, filter_exp, filter_exp2) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  meta.df <- tmp.df %>% distinct(SeqID, Library, Method, Filt, SampleType, Labeler)
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

## generate rarefied dataset:
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


## combine both datasets and add Rarefy label
df$RarefyType <- "unrarefied"
rare.df$RarefyType <- "rarefied"
df <- rbind(df, rare.df)
rm(rare.df)

## make separate plots for rarefied and unrarefied data
## generate summaries for number of ASVs per sample, per pipeline and filtering method and rarefying type
unrare.df.ASVcounts <- df %>% filter(RarefyType=="unrarefied") %>%group_by(Method, Filt, SeqID, Library) %>% summarise(ASVcounts=n())
rare.df.ASVcounts <- df %>% filter(RarefyType=="rarefied") %>%group_by(Method, Filt, SeqID, Library) %>% summarise(ASVcounts=n())

## reset the levels for plot
unrare.df.ASVcounts$Filt <- factor(unrare.df.ASVcounts$Filt,levels = c("basic", "standard", "extra"))
unrare.df.ASVcounts$Method <- factor(unrare.df.ASVcounts$Method,levels = c("dada2", "deblur", "vsearch"))
rare.df.ASVcounts$Filt <- factor(rare.df.ASVcounts$Filt,levels = c("basic", "standard", "extra"))
rare.df.ASVcounts$Method <- factor(rare.df.ASVcounts$Method,levels = c("dada2", "deblur", "vsearch"))


## plot unrarefied first; save as 6_GuanoASVCounts_byFilterMethod_Unrarefied; export at 900x900
ggplot(data = unrare.df.ASVcounts, aes(x = Filt, y = ASVcounts, color=Filt)) +
  geom_jitter(alpha=0.2, width = 0.25) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(Method ~ Library) +
  ylim(0,250) +
  scale_color_manual(values=pal3) +
  labs(title="", x="", y="sequence variants per sample", color="",
       caption = "60 outliers with > 250 sequence variants per sample not shown \n All outliers associated with Vsearch-filtered samples") +
  theme_devon() +
  theme(legend.position="top", legend.text = element_text(size = 12), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank())

## plot rarefied data; save as 6_GuanoASVCounts_byFilterMethod_Rarefied; export at 900x900
ggplot(data = rare.df.ASVcounts, aes(x = Filt, y = ASVcounts, color=Filt)) +
  geom_jitter(alpha=0.2, width = 0.25) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(Method ~ Library) +
  ylim(0,250) +
  scale_color_manual(values=pal3) +
  labs(title="", x="", y="sequence variants per sample", color="") +
  theme_devon() +
  theme(legend.position="top", legend.text = element_text(size = 12), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank())