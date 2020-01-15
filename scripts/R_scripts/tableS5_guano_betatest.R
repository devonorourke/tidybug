library(tidyverse)
library(reshape2)
library(vegan)

## import data:
df <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/all.filtmethods.df.csv.gz")
meta <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/metadata/meta_nomock.csv") %>% 
  select(-StudyID, -SampleType, -Library)
df <- merge(df, meta)
SelectMonths <- c("April", "May", "September", "October")
fox <- df %>% 
  filter(Site=="FOX", 
         SampleType=="sample",
         MonthStart %in% SelectMonths) %>% 
  mutate(Labeler=paste(Method, Filt, MonthStart, sep="-"),
         bigID=paste(SeqID, Method, Filt, MonthStart, sep = "-"))

#rm(df, meta, SelectMonths)

## rarefy data:
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

tmp.dada2.basic <- rare.function(fox, Method=="dada2", Filt=="basic")
tmp.dada2.standard <- rare.function(fox, Method=="dada2", Filt=="standard")
tmp.dada2.extra <- rare.function(fox, Method=="dada2", Filt=="extra")
tmp.deblur.basic <- rare.function(fox, Method=="deblur", Filt=="basic")
tmp.deblur.standard <- rare.function(fox, Method=="deblur", Filt=="standard")
tmp.deblur.extra <- rare.function(fox, Method=="deblur", Filt=="extra")
tmp.vsearch.basic <- rare.function(fox, Method=="vsearch", Filt=="basic")
tmp.vsearch.standard <- rare.function(fox, Method=="vsearch", Filt=="standard")
tmp.vsearch.extra <- rare.function(fox, Method=="vsearch", Filt=="extra")
r.fox.df <- rbind(tmp.dada2.basic, tmp.dada2.standard, tmp.dada2.extra,
                   tmp.deblur.basic, tmp.deblur.standard, tmp.deblur.extra,
                   tmp.vsearch.basic, tmp.vsearch.standard, tmp.vsearch.extra)
rm(list=ls(pattern="tmp.*"))

## run Adonis
## ADONIS function for RAREFIED data from FOX site for select MONTHS
rrarewithdrop <- function(x, sample) {
  rrarefy(x[rowSums(x) >= sample, , drop = FALSE], sample)
}

adonis.function <- function(betatest, binaryval) {
  tmp.meta <- fox %>% distinct(bigID, MonthStart, Filt, Method)
  tmp.mat <- dcast(fox, bigID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$bigID
  tmp.mat$bigID <- NULL
  tmp.mat2 <- data.frame(rrarewithdrop(tmp.mat, 5000))
  rm(tmp.mat)
  namefilter <- row.names(tmp.mat2)
  tmp.meta <- tmp.meta %>% filter(bigID %in% namefilter)  ## redo meta data to drop any samples that were discarded by rarefying
  tmp.betadiv <- vegdist(tmp.mat2, betatest, binary = binaryval)
  tmp.adonis <- adonis(tmp.betadiv ~ Method * Filt * MonthStart, data=tmp.meta) ## including 3-way interaction; if all interactions are non-significant, can test dispersions for each factor separately
  tmp.adonis <- as.data.frame(tmp.adonis$aov.tab)
  tmp.out <- data.frame(tmp.adonis, betatest, binaryval)
  tmp.out$TestGroup <- row.names(tmp.out)
  tmp.out
}

r.adonis.dice <- adonis.function("bray", TRUE)
r.adonis.dice$betatest <- gsub("bray", "dice", r.adonis.dice$betatest)
r.adonis.bray <- adonis.function("bray", FALSE)
r.adonis.mori <- adonis.function("morisita", FALSE)

r.all_adonis <- rbind(r.adonis.dice, r.adonis.bray, r.adonis.mori)
row.names(r.all_adonis) <- NULL
r.all_adonis <- r.all_adonis %>% 
  select(TestGroup, betatest, Df, SumsOfSqs, MeanSqs, F.Model, R2, Pr..F.)

write.csv(r.all_adonis, file="/Users/devonorourke/Documents/nau_projects/guano/mole_ecol_methods_paper/fox_adonis-rarefied.out.csv", 
          quote=FALSE, row.names = FALSE)
