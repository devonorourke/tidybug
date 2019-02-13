library(tidyverse)
library(reshape2)
library(vegan)


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
rm(SelectMonths)


adonis.function <- function(data, filter_exp, filter_exp2, betatest, binaryval) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  Method <- tmp.df %>% distinct(Method)
  Filt <- tmp.df %>% distinct(Filt)
  BetaType <- binaryval
  tmp.meta <- tmp.df %>% distinct(SeqID, MonthStart, WOY)
  tmp.mat <- dcast(tmp.df, SeqID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SeqID
  SeqID <- row.names(tmp.mat)
  tmp.mat$SeqID <- NULL
  tmp.dist <- vegdist(tmp.mat, betatest, binary = binaryval)
  tmp.adonis <- adonis(tmp.dist ~ MonthStart, data=tmp.meta)
  tmp.adonis <- as.data.frame(tmp.adonis$aov.tab)
  tmp.out <- data.frame(tmp.adonis, Method, Filt, BetaType)
  tmp.out$TestGroup <- row.names(tmp.out)
  tmp.out
}


adonis.bray.dada2.basic <- adonis.function(fox, Method=="dada2", Filt=="basic", "bray", TRUE)
adonis.dice.dada2.basic <- adonis.function(fox, Method=="dada2", Filt=="basic", "bray", FALSE)
adonis.bray.dada2.standard <- adonis.function(fox, Method=="dada2", Filt=="standard", "bray", TRUE)
adonis.dice.dada2.standard <- adonis.function(fox, Method=="dada2", Filt=="standard", "bray", FALSE)
adonis.bray.dada2.extra <- adonis.function(fox, Method=="dada2", Filt=="extra", "bray", TRUE)
adonis.dice.dada2.extra <- adonis.function(fox, Method=="dada2", Filt=="extra", "bray", FALSE)

adonis.bray.deblur.basic <- adonis.function(fox, Method=="deblur", Filt=="basic", "bray", TRUE)
adonis.dice.deblur.basic <- adonis.function(fox, Method=="deblur", Filt=="basic", "bray", FALSE)
adonis.bray.deblur.standard <- adonis.function(fox, Method=="deblur", Filt=="standard", "bray", TRUE)
adonis.dice.deblur.standard <- adonis.function(fox, Method=="deblur", Filt=="standard", "bray", FALSE)
adonis.bray.deblur.extra <- adonis.function(fox, Method=="deblur", Filt=="extra", "bray", TRUE)
adonis.dice.deblur.extra <- adonis.function(fox, Method=="deblur", Filt=="extra", "bray", FALSE)

adonis.bray.vsearch.basic <- adonis.function(fox, Method=="vsearch", Filt=="basic", "bray", TRUE)
adonis.dice.vsearch.basic <- adonis.function(fox, Method=="vsearch", Filt=="basic", "bray", FALSE)
adonis.bray.vsearch.standard <- adonis.function(fox, Method=="vsearch", Filt=="standard", "bray", TRUE)
adonis.dice.vsearch.standard <- adonis.function(fox, Method=="vsearch", Filt=="standard", "bray", FALSE)
adonis.bray.vsearch.extra <- adonis.function(fox, Method=="vsearch", Filt=="extra", "bray", TRUE)
adonis.dice.vsearch.extra <- adonis.function(fox, Method=="vsearch", Filt=="extra", "bray", FALSE)

all.adonis.df <- rbind(adonis.dice.vsearch.extra, adonis.dice.vsearch.standard, adonis.dice.vsearch.basic,
                       adonis.bray.vsearch.extra, adonis.bray.vsearch.standard, adonis.bray.vsearch.basic,
                       adonis.dice.dada2.extra, adonis.dice.dada2.standard, adonis.dice.dada2.basic,
                       adonis.bray.dada2.extra, adonis.bray.dada2.standard, adonis.bray.dada2.basic,
                       adonis.dice.deblur.extra, adonis.dice.deblur.standard, adonis.dice.deblur.basic,
                       adonis.bray.deblur.extra, adonis.bray.deblur.standard, adonis.bray.deblur.basic)
row.names(all.adonis.df) <- NULL
all.adonis.df$BetaType[which(all.adonis.df$BetaType==TRUE)] <- "bray"
all.adonis.df$BetaType[which(all.adonis.df$BetaType==FALSE)] <- "dice"

rm(adonis.dice.vsearch.extra, adonis.dice.vsearch.standard, adonis.dice.vsearch.basic,
   adonis.bray.vsearch.extra, adonis.bray.vsearch.standard, adonis.bray.vsearch.basic,
   adonis.dice.dada2.extra, adonis.dice.dada2.standard, adonis.dice.dada2.basic,
   adonis.bray.dada2.extra, adonis.bray.dada2.standard, adonis.bray.dada2.basic,
   adonis.dice.deblur.extra, adonis.dice.deblur.standard, adonis.dice.deblur.basic,
   adonis.bray.deblur.extra, adonis.bray.deblur.standard, adonis.bray.deblur.basic)

write.csv(all.adonis.df, file="~/Repos/tidybug/data/text_tables/adonis.out.csv", quote=FALSE, row.names = FALSE)
