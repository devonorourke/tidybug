library(tidyverse)
library(reshape2)
library(vegan)
#library(cowplot)

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
meta <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/metadata/meta_nomock.csv/")
meta <- meta %>% select(-StudyID, -SampleType, -Library)
df <- merge(df, meta)
rm(meta)
fox <- df %>% filter(Site=="FOX")
rm(df)
## only selected months
SelectMonths <- c("April", "May", "September", "October")
fox <- fox %>% filter(MonthStart %in% SelectMonths)


##function applied to calculate indices for beta diversity
betadiv.function <- function(data, filter_exp, filter_exp2, arg3) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  tmp.mat <- dcast(tmp.df, SeqID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SeqID
  SeqID <- row.names(tmp.mat)
  tmp.mat$SeqID <- NULL
  vegdist(tmp.mat, arg3)
}

## calculate beta distance - index defined in output name
## for bray-curtis
bray.dada2.basic <- betadiv.function(fox, Method=="dada2", Filt=="basic", "bray")
bray.dada2.standard <- betadiv.function(fox, Method=="dada2", Filt=="standard", "bray")
bray.dada2.extra <- betadiv.function(fox, Method=="dada2", Filt=="extra", "bray")
bray.deblur.basic <- betadiv.function(fox, Method=="deblur", Filt=="basic", "bray")
bray.deblur.standard <- betadiv.function(fox, Method=="deblur", Filt=="standard", "bray")
bray.deblur.extra <- betadiv.function(fox, Method=="deblur", Filt=="extra", "bray")
bray.vsearch.basic <- betadiv.function(fox, Method=="vsearch", Filt=="basic", "bray")
bray.vsearch.standard <- betadiv.function(fox, Method=="vsearch", Filt=="standard", "bray")
bray.vsearch.extra <- betadiv.function(fox, Method=="vsearch", Filt=="extra", "bray")

## for dice-sorensen
## modifying function to set quantiative index calculation to FALSE
betadivBINARY.function <- function(data, filter_exp, filter_exp2, arg3) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  tmp.mat <- dcast(tmp.df, SeqID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SeqID
  SeqID <- row.names(tmp.mat)
  tmp.mat$SeqID <- NULL
  vegdist(tmp.mat, arg3, binary=TRUE)
}

dice.dada2.basic <- betadivBINARY.function(fox, Method=="dada2", Filt=="basic", "bray")
dice.dada2.standard <- betadivBINARY.function(fox, Method=="dada2", Filt=="standard", "bray")
dice.dada2.extra <- betadivBINARY.function(fox, Method=="dada2", Filt=="extra", "bray")
dice.deblur.basic <- betadivBINARY.function(fox, Method=="deblur", Filt=="basic", "bray")
dice.deblur.standard <- betadivBINARY.function(fox, Method=="deblur", Filt=="standard", "bray")
dice.deblur.extra <- betadivBINARY.function(fox, Method=="deblur", Filt=="extra", "bray")
dice.vsearch.basic <- betadivBINARY.function(fox, Method=="vsearch", Filt=="basic", "bray")
dice.vsearch.standard <- betadivBINARY.function(fox, Method=="vsearch", Filt=="standard", "bray")
dice.vsearch.extra <- betadivBINARY.function(fox, Method=="vsearch", Filt=="extra", "bray")

## function to perform adonis test for each distance matrix
## generates needed metadata internally
adonis.function <- function(data, filter_exp, filter_exp2, arg3) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  tmp.meta <- 
  tmp.mat <- dcast(tmp.df, SeqID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SeqID
  SeqID <- row.names(tmp.mat)
  tmp.mat$SeqID <- NULL
  vegdist(tmp.mat, arg3)
}

## could add in the following function within the above function to generate adonis output on the fly...
tmp.df <- fox %>% filter(Method=="deblur" & Filt=="extra")
tmp.meta <- tmp.df %>% distinct(SeqID, MonthStart)
adonis(deblur.extra ~ MonthStart, data=tmp.meta)

tmp.df <- fox %>% filter(Method=="deblur" & Filt=="basic")
tmp.meta <- tmp.df %>% distinct(SeqID, MonthStart)
adonis(deblur.basic ~ MonthStart, data=tmp.meta)


## What about plotting dispersion?
tmp.df <- fox %>% filter(Method=="deblur" & Filt=="extra")
tmp.meta <- tmp.df %>% distinct(SeqID, MonthStart)
groups <- as.character(tmp.meta$MonthStart)
tmp.mat <- dcast(tmp.df, SeqID ~ HashID, value.var = "Reads", fill=0)
row.names(tmp.mat) <- tmp.mat$SeqID
SeqID <- row.names(tmp.mat)
tmp.mat$SeqID <- NULL
braydisp.deblur.extra <- betadisper(bray.deblur.extra, groups, type = c("median"))
anova(braydisp.deblur.extra)
## basic plot
plot(braydisp.deblur.extra)
## no labels, but added ellipses
plot(braydisp.deblur.extra, label = F, hull = F, ellipse = T)
## boxplot
boxplot(braydisp.deblur.extra)
## significane tests
TukeyHSD(braydisp.deblur.extra)


## What about plotting dispersion?
dicedisp.deblur.extra <- betadisper(dice.deblur.extra, groups, type = c("median"))
anova(dicedisp.deblur.extra)
## basic plot
plot(dicedisp.deblur.extra)
## no labels, but added ellipses
plot(dicedisp.deblur.extra, label = F, hull = F, ellipse = T)
## boxplot
boxplot(dicedisp.deblur.extra)
## significane tests
TukeyHSD(dicedisp.deblur.extra)


#### enhanced plotting
## gather the "scores" -- PCoA points...
dicedisp.deblur.extra_sco <- as.data.frame(scores(dicedisp.deblur.extra, display = "sites"))
dicedisp.deblur.extra_sco$SeqID <- row.names(dicedisp.deblur.extra_sco)
row.names(dicedisp.deblur.extra_sco) <- NULL
colnames(dicedisp.deblur.extra_sco) <- c("pc1end", "pc2end", "SeqID")
dicedisp.deblur.extra_sco$Type <- "sample"
dicedisp.deblur.extra_sco <- merge(dicedisp.deblur.extra_sco, tmp.meta)
dicedisp.deblur.extra_sco
## gather the centroids
dicedisp.deblur.extra_centroids <- as.data.frame(scores(dicedisp.deblur.extra, display="centroids"))
dicedisp.deblur.extra_centroids$SeqID <- row.names(dicedisp.deblur.extra_centroids)
row.names(dicedisp.deblur.extra_centroids) <- NULL
colnames(dicedisp.deblur.extra_centroids) <- c("pc1start", "pc2start", "MonthStart")
dicedisp.deblur.extra_centroids

## combine data as single data.frame
dicedisp.deblur.extra.dfplot <- merge(dicedisp.deblur.extra_sco, dicedisp.deblur.extra_centroids)
dicedisp.deblur.extra.dfplot

## plotting:
## define color palette:

## plot
ggplot(dicedisp.deblur.extra.dfplot, aes(x=pc1end, y=pc2end, shape=MonthStart, color=MonthStart)) + 
  geom_point(data=dicedisp.deblur.extra.dfplot, aes(x=pc1end, y=pc2end, color=MonthStart, shape=MonthStart)) +
  geom_segment(data=dicedisp.deblur.extra.dfplot, aes(x=pc1end, y=pc2end, xend=pc1start, yend=pc2start), alpha=0.5) +
  scale_color_manual() + 
  theme_devon()

