library(tidyverse)
library(reshape2)
library(vegan)

## import data:
df <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/all.filtmethods.df.csv.gz")
mock <- df %>% filter(SampleType == "mock")

## Hull Numbers function applied to calculate diversity measures: observed OTUs, simpson, and shannon
## unrarefied data only...
hillnumber.function <- function(data, filter_exp, filter_exp2) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.mock <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  meta.mock <- tmp.mock %>% distinct(SeqID, Library, SampleType)
  Method <- tmp.mock %>% distinct(Method)
  Filt <- tmp.mock %>% distinct(Filt)
  tmp.mat <- dcast(tmp.mock, SeqID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SeqID
  SeqID <- row.names(tmp.mat)
  tmp.mat$SeqID <- NULL
  tmp.hill <- renyi(tmp.mat, scales = c(0,1,2), hill=TRUE)
  tmp <- data.frame(tmp.hill, SeqID, Method, Filt, row.names = NULL)
  colnames(tmp)[1:3] <- c("q=0", "q=1", "q=2")
  tmp <- gather(tmp, key="Hill_qType", value = "Hill_value", c('q=0', 'q=1', 'q=2'))
  tmp <- merge(tmp, meta.mock)
}

## calculate Hill Numbers for unrarefied data
dada2.basic <- hillnumber.function(mock, Method=="dada2", Filt=="basic")
dada2.standard <- hillnumber.function(mock, Method=="dada2", Filt=="standard")
dada2.extra <- hillnumber.function(mock, Method=="dada2", Filt=="extra")
deblur.basic <- hillnumber.function(mock, Method=="deblur", Filt=="basic")
deblur.standard <- hillnumber.function(mock, Method=="deblur", Filt=="standard")
deblur.extra <- hillnumber.function(mock, Method=="deblur", Filt=="extra")
vsearch.basic <- hillnumber.function(mock, Method=="vsearch", Filt=="basic")
vsearch.standard <- hillnumber.function(mock, Method=="vsearch", Filt=="standard")
vsearch.extra <- hillnumber.function(mock, Method=="vsearch", Filt=="extra")

## merge into single dataframe
all.mock.hill <- rbind(dada2.basic, dada2.standard, dada2.extra, deblur.basic, deblur.standard, deblur.extra, vsearch.basic, vsearch.standard, vsearch.extra)
rm(dada2.basic, dada2.standard, dada2.extra, deblur.basic, deblur.standard, deblur.extra, vsearch.basic, vsearch.standard, vsearch.extra)
rm(df, mock)

## add labeler to make facet labels easier to understand
## add $Labeler and reformat from character to factor
all.mock.hill$Labeler <- paste(all.mock.hill$Method, all.mock.hill$Filt, sep="-")
all.mock.hill$Labeler <- as.factor(all.mock.hill$Labeler)

## split mock by Hill q-value:
q0.hill.mock <- all.mock.hill %>% filter(Hill_qType == "q=0")
q1.hill.mock <- all.mock.hill %>% filter(Hill_qType == "q=1")
q2.hill.mock <- all.mock.hill %>% filter(Hill_qType == "q=2")

## 3way anova, one per Hill Number
q0.anova <- aov(Hill_value ~ Method * Filt, data=q0.hill.mock)
summary(q0.anova)
  ## Method and Filt and interaction Method*Filt all higly significant
q1.anova <- aov(Hill_value ~ Method * Filt, data=q1.hill.mock)
summary(q1.anova)
  ## same as q=0; both factors and interaction significant
q2.anova <- aov(Hill_value ~  Method * Filt, data=q2.hill.mock)
summary(q2.anova)
  ## Only Method significant - not Filt or interaction term...
  ## when evaluating dominant sequences, the filtering Method is what matters


## dunn test and resulting pairwise correlation heatmap:
## big messy function to generate the data.frame for plotting a correlation matrix of the pairwise values from Dunn test
corrplot.function <- function(data, filter_exp) {
  tmp.dunn = dunnTest(Hill_value ~ Labeler, data=data, method="bh")
  tmp.dunn.mock <- tmp.dunn$res
  tmp.dunn.mock <- separate(data=tmp.dunn.mock, col = "Comparison", into=c("Group1", "Group2"), sep=" - ")
  tmp.dunn.mock$P.adj <- round(tmp.dunn.mock$P.adj,3)
  tmp.dunn.mat <- dcast(tmp.dunn.mock, Group1 ~ Group2, value.var = "P.adj", drop = FALSE)
  row.names(tmp.dunn.mat) <- tmp.dunn.mat$Group1
  tmp.dunn.mat$Group1 <- NA
  tmp.dunn.mat <- as.matrix(tmp.dunn.mat)
  botrow <- rep(NA, length(colnames(tmp.dunn.mat)))
  tmp.dunn.mat <- rbind(tmp.dunn.mat, botrow)
  lengthcols <- length(colnames(tmp.dunn.mat))
  selectrow <- colnames(tmp.dunn.mat)[lengthcols]
  rownames(tmp.dunn.mat)[lengthcols] <- selectrow
  rownames(tmp.dunn.mat)[1] -> selectcol
  colnames(tmp.dunn.mat)[1] <- selectcol
  colnames(tmp.dunn.mat) <- rownames(tmp.dunn.mat)
  tmp.dunn.mat[lower.tri(tmp.dunn.mat)] = t(tmp.dunn.mat)[lower.tri(tmp.dunn.mat)]
  tmp.dunn.mat[is.na(tmp.dunn.mat)] <- 1
  tmp2 <- tmp.dunn.mat + 10
  tmp3 <- tril(tmp2)
  tmp3 <- tmp3 - 10
  tmp4 <- as.matrix(tmp3)
  melt(tmp4, varnames = c("Group1", "Group2"), value.name = "P.adj") %>% filter(P.adj >= 0)
}


q0.dunn = dunnTest(Hill_value ~ Labeler, data=q0.hill.mock, method="bh")
q0.dunn.mock <- q0.dunn$res
q0.dunn.mock <- separate(data=q0.dunn.mock, col = "Comparison", into=c("Group1", "Group2"), sep=" - ")
q0.dunn.mock$P.adj <- round(q0.dunn.mock$P.adj,3)

q1.dunn = dunnTest(Hill_value ~ Labeler, data=q1.hill.mock, method="bh")
q1.dunn.mock <- q1.dunn$res
q1.dunn.mock <- separate(data=q1.dunn.mock, col = "Comparison", into=c("Group1", "Group2"), sep=" - ")
q1.dunn.mock$P.adj <- round(q1.dunn.mock$P.adj,3)

q2.dunn = dunnTest(Hill_value ~ Labeler, data=q2.hill.mock, method="bh")
q2.dunn.mock <- q2.dunn$res
q2.dunn.mock <- separate(data=q2.dunn.mock, col = "Comparison", into=c("Group1", "Group2"), sep=" - ")
q2.dunn.mock$P.adj <- round(q2.dunn.mock$P.adj,3)
