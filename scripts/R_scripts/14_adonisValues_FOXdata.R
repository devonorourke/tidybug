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
fox$Labeler <- paste(fox$Method, fox$Filt, fox$MonthStart, sep="-")
fox$bigID <- paste(fox$SeqID, fox$Method, fox$Filt, fox$MonthStart, sep="-")
str(tmp.meta)

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
r.all_adonis$binaryval <- NULL

#write.csv(r.all_adonis, file="~/Repos/tidybug/data/text_tables/adonis_FOX.out.csv", quote=FALSE, row.names = FALSE)