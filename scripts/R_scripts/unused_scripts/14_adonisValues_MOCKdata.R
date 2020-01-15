library(tidyverse)
library(reshape2)
library(vegan)


## import data:
df <- read_csv("~/Repos/tidybug/data/text_tables/all.filtmethods.df.csv.gz")
meta <- read_csv("~/Repos/tidybug/data/metadata/meta_nomock.csv")
mock <- df %>% filter(SampleType == "mock")
mock$Labeler <- paste(mock$Method, mock$Filt, sep="-")
mock$bigID <- paste(mock$SeqID, mock$Method, mock$Filt, sep="-")
mock.meta <- mock %>% group_by(Method, Filt) %>% distinct(bigID)

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

tmp.dada2.basic <- rare.function(mock, Method=="dada2", Filt=="basic")
tmp.dada2.standard <- rare.function(mock, Method=="dada2", Filt=="standard")
tmp.dada2.extra <- rare.function(mock, Method=="dada2", Filt=="extra")
tmp.deblur.basic <- rare.function(mock, Method=="deblur", Filt=="basic")
tmp.deblur.standard <- rare.function(mock, Method=="deblur", Filt=="standard")
tmp.deblur.extra <- rare.function(mock, Method=="deblur", Filt=="extra")
tmp.vsearch.basic <- rare.function(mock, Method=="vsearch", Filt=="basic")
tmp.vsearch.standard <- rare.function(mock, Method=="vsearch", Filt=="standard")
tmp.vsearch.extra <- rare.function(mock, Method=="vsearch", Filt=="extra")
r.mock.df <- rbind(tmp.dada2.basic, tmp.dada2.standard, tmp.dada2.extra,
                   tmp.deblur.basic, tmp.deblur.standard, tmp.deblur.extra,
                   tmp.vsearch.basic, tmp.vsearch.standard, tmp.vsearch.extra)
rm(list=ls(pattern="tmp.*"))

## run Adonis
adonis.function <- function(data, betatest, binaryval) {
  tmp.meta <- data %>% group_by(Method, Filt) %>% distinct(bigID)
  tmp.mat <- dcast(data, bigID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$bigID
  bigID <- row.names(tmp.mat)
  tmp.mat$bigID <- NULL
  tmp.dist <- vegdist(tmp.mat, betatest, binary = binaryval)
  tmp.adonis <- adonis(tmp.dist ~ Method * Filt, data=tmp.meta)
  tmp.adonis <- as.data.frame(tmp.adonis$aov.tab)
  tmp.out <- data.frame(tmp.adonis, betatest)
  tmp.out$TestGroup <- row.names(tmp.out)
  tmp.out
}


r.mock.df$Labeler <- paste(r.mock.df$Method, r.mock.df$Filt, sep="-")
r.mock.df$bigID <- paste(r.mock.df$SeqID, r.mock.df$Method, r.mock.df$Filt, sep="-")
r.mock.df.meta <- r.mock.df %>% group_by(Method, Filt) %>% distinct(bigID)


r.adonis.dice <- adonis.function(r.mock.df, "bray", TRUE)
r.adonis.dice$betatest <- gsub("bray", "dice", r.adonis.dice$betatest)
r.adonis.bray <- adonis.function(mock, "bray", FALSE)
r.adonis.mori <- adonis.function(mock, "morisita", FALSE)

r.all_adonis <- rbind(r.adonis.dice, r.adonis.bray, r.adonis.mori)
row.names(r.all_adonis) <- NULL

write.csv(r.all_adonis, file="~/Repos/tidybug/data/text_tables/mock_adonis-rarefied.out.csv", 
          quote=FALSE, row.names = FALSE)


#### run the same analyses on unrarefied data? 
## demonstrates how with Dice measure, no main effects significant because of unequal sampling depths
adonis.dice <- adonis.function(mock, "bray", TRUE)
adonis.dice$betatest <- gsub("bray", "dice", adonis.dice$betatest)
adonis.bray <- adonis.function(mock, "bray", FALSE)
adonis.mori <- adonis.function(mock, "morisita", FALSE)

all_adonis <- rbind(adonis.dice, adonis.bray, adonis.mori)
row.names(all_adonis) <- NULL

write.csv(all_adonis, file="~/Repos/tidybug/data/text_tables/mock_adonis.out.csv", 
          quote=FALSE, row.names = FALSE)