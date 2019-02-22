library(tidyverse)
library(reshape2)
library(vegan)
library(cowplot)
library(FSA)
library(Matrix)

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
rrarewithdrop <- function(x, sample) {
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

## Hull Numbers function applied to calculate diversity measures: observed OTUs, simpson, and shannon
hillnumber.function <- function(data, filter_exp, filter_exp2, filter_exp3) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  filter_exp_enq3 <- enquo(filter_exp3)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2) %>% filter(!!filter_exp_enq3)
  meta.df <- tmp.df %>% distinct(SeqID, Library, SampleType)
  Method <- tmp.df %>% distinct(Method)
  Filt <- tmp.df %>% distinct(Filt)
  RarefyType <- tmp.df %>% distinct(RarefyType)
  tmp.mat <- dcast(tmp.df, SeqID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SeqID
  SeqID <- row.names(tmp.mat)
  tmp.mat$SeqID <- NULL
  tmp.hill <- renyi(tmp.mat, scales = c(0,1,2), hill=TRUE)
  tmp <- data.frame(tmp.hill, SeqID, Method, Filt, RarefyType, row.names = NULL)
  colnames(tmp)[1:3] <- c("q=0", "q=1", "q=2")
  tmp <- gather(tmp, key="Hill_qType", value = "Hill_value", c('q=0', 'q=1', 'q=2'))
  merge(tmp, meta.df)
}

## calculate Hill Numbers for unrarefied and rarefied data
unrarefied.dada2.basic <- hillnumber.function(df, Method=="dada2", Filt=="basic", RarefyType=="unrarefied")
unrarefied.dada2.standard <- hillnumber.function(df, Method=="dada2", Filt=="standard", RarefyType=="unrarefied")
unrarefied.dada2.extra <- hillnumber.function(df, Method=="dada2", Filt=="extra", RarefyType=="unrarefied")
unrarefied.deblur.basic <- hillnumber.function(df, Method=="deblur", Filt=="basic", RarefyType=="unrarefied")
unrarefied.deblur.standard <- hillnumber.function(df, Method=="deblur", Filt=="standard", RarefyType=="unrarefied")
unrarefied.deblur.extra <- hillnumber.function(df, Method=="deblur", Filt=="extra", RarefyType=="unrarefied")
unrarefied.vsearch.basic <- hillnumber.function(df, Method=="vsearch", Filt=="basic", RarefyType=="unrarefied")
unrarefied.vsearch.standard <- hillnumber.function(df, Method=="vsearch", Filt=="standard", RarefyType=="unrarefied")
unrarefied.vsearch.extra <- hillnumber.function(df, Method=="vsearch", Filt=="extra", RarefyType=="unrarefied")

rarefied.dada2.basic <- hillnumber.function(df, Method=="dada2", Filt=="basic", RarefyType=="rarefied")
rarefied.dada2.standard <- hillnumber.function(df, Method=="dada2", Filt=="standard", RarefyType=="rarefied")
rarefied.dada2.extra <- hillnumber.function(df, Method=="dada2", Filt=="extra", RarefyType=="rarefied")
rarefied.deblur.basic <- hillnumber.function(df, Method=="deblur", Filt=="basic", RarefyType=="rarefied")
rarefied.deblur.standard <- hillnumber.function(df, Method=="deblur", Filt=="standard", RarefyType=="rarefied")
rarefied.deblur.extra <- hillnumber.function(df, Method=="deblur", Filt=="extra", RarefyType=="rarefied")
rarefied.vsearch.basic <- hillnumber.function(df, Method=="vsearch", Filt=="basic", RarefyType=="rarefied")
rarefied.vsearch.standard <- hillnumber.function(df, Method=="vsearch", Filt=="standard", RarefyType=="rarefied")
rarefied.vsearch.extra <- hillnumber.function(df, Method=="vsearch", Filt=="extra", RarefyType=="rarefied")

# rm(df)

## merge into single dataframe
all.guano.hill <- rbind(unrarefied.dada2.basic, unrarefied.dada2.standard, unrarefied.dada2.extra, 
                        unrarefied.deblur.basic, unrarefied.deblur.standard, unrarefied.deblur.extra, 
                        unrarefied.vsearch.basic, unrarefied.vsearch.standard, unrarefied.vsearch.extra,
                        rarefied.dada2.basic, rarefied.dada2.standard, rarefied.dada2.extra, 
                        rarefied.deblur.basic, rarefied.deblur.standard, rarefied.deblur.extra, 
                        rarefied.vsearch.basic, rarefied.vsearch.standard, rarefied.vsearch.extra)
## add Labeler to group data
all.guano.hill$Labeler <- paste(all.guano.hill$Method, all.guano.hill$Filt, all.guano.hill$Library, sep = "-")
all.guano.hill$Labeler <- as.factor(all.guano.hill$Labeler)

## cleanup:
rm(list=ls(pattern = "rarefied"))

## split mock by Hill q-value and RarefyType:
u.q0.hill.guano <- all.guano.hill %>% filter(Hill_qType == "q=0" & RarefyType == "unrarefied")
r.q0.hill.guano <- all.guano.hill %>% filter(Hill_qType == "q=0" & RarefyType == "rarefied")
u.q1.hill.guano <- all.guano.hill %>% filter(Hill_qType == "q=1" & RarefyType == "unrarefied")
r.q1.hill.guano <- all.guano.hill %>% filter(Hill_qType == "q=1" & RarefyType == "rarefied")
u.q2.hill.guano <- all.guano.hill %>% filter(Hill_qType == "q=2" & RarefyType == "unrarefied")
r.q2.hill.guano <- all.guano.hill %>% filter(Hill_qType == "q=2" & RarefyType == "rarefied")

## 3way anova, one per RarefyType and Hill Number
u.q0.anova <- aov(Hill_value ~ Library * Method * Filt, data=u.q0.hill.guano)
summary(u.q0.anova)
  ## all main effects and all interactions significant
  ## Method > Filt > Library for F value

r.q0.anova <- aov(Hill_value ~ Library * Method * Filt, data=r.q0.hill.guano)
summary(r.q0.anova)
  ## only the Library:Method:Filt interaction is NOT significant
  ## Filt >> Method or Library for F value

u.q1.anova <- aov(Hill_value ~ Library * Method * Filt, data=u.q1.hill.guano)
summary(u.q1.anova)
  ## all main effects and interactions significantly different;
  ## Filt > Library >> Method for F-value

r.q1.anova <- aov(Hill_value ~ Library * Method * Filt, data=r.q1.hill.guano)
summary(r.q1.anova)
  ## Method:Filt and 3-way interaction NOT significant
  ## Library and Method >> Filt

u.q2.anova <- aov(Hill_value ~ Library * Method * Filt, data=u.q2.hill.guano)
summary(u.q2.anova)
  ## All main efects significant; 3 way interaction NOT significant
  ## Filt > Library >> Method

r.q2.anova <- aov(Hill_value ~ Library * Method * Filt, data=r.q2.hill.guano)
summary(r.q2.anova)
## All main efects significant; 3 way interaction NOT significant
## Filt > Library >> Method


## dunn test and resulting pairwise correlation heatmap:
## big messy function to generate the data.frame for plotting a correlation matrix of the pairwise values from Dunn test
corrplot.function <- function(data, filter_exp) {
  tmp.dunn = dunnTest(Hill_value ~ Labeler, data=data, method="bh")
  tmp.dunn.df <- tmp.dunn$res
  tmp.dunn.df <- separate(data=tmp.dunn.df, col = "Comparison", into=c("Group1", "Group2"), sep=" - ")
  tmp.dunn.df$P.adj <- round(tmp.dunn.df$P.adj,2)
  tmp.dunn.mat <- dcast(tmp.dunn.df, Group1 ~ Group2, value.var = "P.adj", drop = FALSE)
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
  DesiredOrder <- c("dada2-basic-libA", "dada2-basic-libB", "dada2-basic-libC", "dada2-basic-libD","dada2-standard-libA", "dada2-standard-libB", "dada2-standard-libC", "dada2-standard-libD","dada2-extra-libA", "dada2-extra-libB", "dada2-extra-libC", "dada2-extra-libD",
                    "deblur-basic-libA", "deblur-basic-libB", "deblur-basic-libC", "deblur-basic-libD","deblur-standard-libA", "deblur-standard-libB", "deblur-standard-libC", "deblur-standard-libD","deblur-extra-libA", "deblur-extra-libB", "deblur-extra-libC", "deblur-extra-libD",
                    "vsearch-basic-libA", "vsearch-basic-libB", "vsearch-basic-libC", "vsearch-basic-libD","vsearch-standard-libA", "vsearch-standard-libB", "vsearch-standard-libC", "vsearch-standard-libD","vsearch-extra-libA", "vsearch-extra-libB", "vsearch-extra-libC", "vsearch-extra-libD")
  tmp.dunn.mat2 <- tmp.dunn.mat[DesiredOrder,DesiredOrder]
  tmp2 <- tmp.dunn.mat2 + 10
  tmp3 <- tril(tmp2)
  tmp3 <- tmp3 - 10
  tmp4 <- as.matrix(tmp3)
  tmp5 <- melt(tmp4, varnames = c("Group1", "Group2"), value.name = "P.adj") %>% filter(P.adj >= 0)
  tmp5$textval <- tmp5$P.adj
  tmp5$textval <- as.character(round(tmp5$textval,2))
  tmp5
}

## generate pairwise charts with self:self and combine in plots for rarefied and unrarefied
u.q0.dunn.guano <- corrplot.function(u.q0.hill.guano)
r.q0.dunn.guano <- corrplot.function(r.q0.hill.guano)
u.q1.dunn.guano <- corrplot.function(u.q1.hill.guano)
r.q1.dunn.guano <- corrplot.function(r.q1.hill.guano)
u.q2.dunn.guano <- corrplot.function(u.q2.hill.guano)
r.q2.dunn.guano <- corrplot.function(r.q2.hill.guano)

## combine rarefied/unrarefied data to determine which P.adj values differed with respect to significance
## defining difference as notable when one group significance is less than 0.05, and other is greater than 0.1
## q=0
colnames(u.q0.dunn.guano) <- paste("u", colnames(u.q0.dunn.guano), sep = ".")
colnames(r.q0.dunn.guano) <- paste("r", colnames(r.q0.dunn.guano), sep = ".")
all.q0.dunn.guano <- cbind(u.q0.dunn.guano, r.q0.dunn.guano)
all.q0.dunn.guano$SigLoss <- all.q0.dunn.guano$u.P.adj < 0.05 & all.q0.dunn.guano$r.P.adj > 0.1 ## "loss" is relative to action of rarefying the data
all.q0.dunn.guano$SigGain <- all.q0.dunn.guano$u.P.adj > 0.1 & all.q0.dunn.guano$r.P.adj < 0.05 ## "gain" is relative to action of rarefying

## q=1
colnames(u.q1.dunn.guano) <- paste("u", colnames(u.q1.dunn.guano), sep = ".")
colnames(r.q1.dunn.guano) <- paste("r", colnames(r.q1.dunn.guano), sep = ".")
all.q1.dunn.guano <- cbind(u.q1.dunn.guano, r.q1.dunn.guano)
all.q1.dunn.guano$SigLoss <- all.q1.dunn.guano$u.P.adj < 0.05 & all.q1.dunn.guano$r.P.adj > 0.1 ## "loss" is relative to action of rarefying the data
all.q1.dunn.guano$SigGain <- all.q1.dunn.guano$u.P.adj > 0.1 & all.q1.dunn.guano$r.P.adj < 0.05 ## "gain" is relative to action of rarefying

## q=2
colnames(u.q2.dunn.guano) <- paste("u", colnames(u.q2.dunn.guano), sep = ".")
colnames(r.q2.dunn.guano) <- paste("r", colnames(r.q2.dunn.guano), sep = ".")
all.q2.dunn.guano <- cbind(u.q2.dunn.guano, r.q2.dunn.guano)
all.q2.dunn.guano$SigLoss <- all.q2.dunn.guano$u.P.adj < 0.05 & all.q2.dunn.guano$r.P.adj > 0.1 ## "loss" is relative to action of rarefying the data
all.q2.dunn.guano$SigGain <- all.q2.dunn.guano$u.P.adj > 0.1 & all.q2.dunn.guano$r.P.adj < 0.05 ## "gain" is relative to action of rarefying


## plots; one per q-value
## save as 12_pairwiseDunn_guano_q0; export at 1000x900
ggplot() +
  geom_tile(data=all.q0.dunn.guano, aes(x=u.Group2, y=u.Group1, fill=u.P.adj), color="black", size=0.25, alpha=0.5) +
  #geom_tile(data=all.q0.dunn.guano %>% filter(SigGain == TRUE), aes(x=u.Group2, y=u.Group1, fill=u.P.adj), color="red", size=1) +
  #geom_tile(data=all.q0.dunn.guano %>% filter(SigLoss == TRUE), aes(x=u.Group2, y=u.Group1, fill=u.P.adj), color="blue", size=1) +
  geom_text(data=all.q0.dunn.guano, aes(x=u.Group2, y=u.Group1, label=u.textval, fontface="bold"), size=2, color="black") +
  geom_tile(data=all.q0.dunn.guano, aes(x=r.Group1, y=r.Group2, fill=r.P.adj), color="black", size=0.25, alpha=0.5) +
  geom_tile(data=all.q0.dunn.guano %>% filter(SigGain == TRUE), aes(x=r.Group1, y=r.Group2, fill=r.P.adj), color="blue", size=1, alpha=0.5) +
  geom_tile(data=all.q0.dunn.guano %>% filter(SigLoss == TRUE), aes(x=r.Group1, y=r.Group2, fill=r.P.adj), color="red", size=1, alpha=0.5) +
  geom_text(data=all.q0.dunn.guano, aes(x=r.Group1, y=r.Group2, label=r.textval, fontface="bold"), size=2, color="black") +
  scale_fill_viridis_c(direction = -1, option = "cividis", begin=0.3) +
  #scale_color_viridis_d(direction = 1, option = "magma") +
  guides(colour=FALSE) +
  #scale_y_discrete(position = "right") +
  #scale_x_discrete(position = "bottom") +
  labs(x="", y="", fill="BH-adjusted\np.value") +
  theme_devon() +
  theme(axis.text.x = element_text(angle=22.5, hjust=1, size = 7),
        axis.text.y = element_text(size=7), legend.position = "right",
        panel.background = element_blank(), panel.border = element_blank(), panel.grid = element_blank(),
        plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  theme(legend.position = "top")
#


## save as 12_pairwiseDunn_guano_q1; export at 1000x900
ggplot() +
  geom_tile(data=all.q1.dunn.guano, aes(x=u.Group2, y=u.Group1, fill=u.P.adj), color="black", size=0.25, alpha=0.5) +
  geom_text(data=all.q1.dunn.guano, aes(x=u.Group2, y=u.Group1, label=u.textval, fontface="bold"), size=2, color="black") +
  geom_tile(data=all.q1.dunn.guano, aes(x=r.Group1, y=r.Group2, fill=r.P.adj), color="black", size=0.25, alpha=0.5) +
  geom_tile(data=all.q1.dunn.guano %>% filter(SigGain == TRUE), aes(x=r.Group1, y=r.Group2, fill=r.P.adj), color="blue", size=1, alpha=0.5) +
  geom_tile(data=all.q1.dunn.guano %>% filter(SigLoss == TRUE), aes(x=r.Group1, y=r.Group2, fill=r.P.adj), color="red", size=1, alpha=0.5) +
  geom_text(data=all.q1.dunn.guano, aes(x=r.Group1, y=r.Group2, label=r.textval, fontface="bold"), size=2, color="black") +
  scale_fill_viridis_c(direction = -1, option = "cividis", begin=0.3) +
  guides(colour=FALSE) +
  labs(x="", y="", fill="BH-adjusted\np.value") +
  theme_devon() +
  theme(axis.text.x = element_text(angle=22.5, hjust=1, size = 7),
        axis.text.y = element_text(size=7), legend.position = "right",
        panel.background = element_blank(), panel.border = element_blank(), panel.grid = element_blank(),
        plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  theme(legend.position = "top")


## save as 12_pairwiseDunn_guano_q2; export at 1000x900
ggplot() +
  geom_tile(data=all.q2.dunn.guano, aes(x=u.Group2, y=u.Group1, fill=u.P.adj), color="black", size=0.25, alpha=0.5) +
  geom_text(data=all.q2.dunn.guano, aes(x=u.Group2, y=u.Group1, label=u.textval, fontface="bold"), size=2, color="black") +
  geom_tile(data=all.q2.dunn.guano, aes(x=r.Group1, y=r.Group2, fill=r.P.adj), color="black", size=0.25, alpha=0.5) +
  geom_tile(data=all.q2.dunn.guano %>% filter(SigGain == TRUE), aes(x=r.Group1, y=r.Group2, fill=r.P.adj), color="blue", size=1, alpha=0.5) +
  geom_tile(data=all.q2.dunn.guano %>% filter(SigLoss == TRUE), aes(x=r.Group1, y=r.Group2, fill=r.P.adj), color="red", size=1, alpha=0.5) +
  geom_text(data=all.q2.dunn.guano, aes(x=r.Group1, y=r.Group2, label=r.textval, fontface="bold"), size=2, color="black") +
  scale_fill_viridis_c(direction = -1, option = "cividis", begin=0.3) +
  guides(colour=FALSE) +
  labs(x="", y="", fill="BH-adjusted\np.value") +
  theme_devon() +
  theme(axis.text.x = element_text(angle=22.5, hjust=1, size = 7),
        axis.text.y = element_text(size=7), legend.position = "right",
        panel.background = element_blank(), panel.border = element_blank(), panel.grid = element_blank(),
        plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  theme(legend.position = "top")


## Summarizing Dunn comparisons a few different ways
## function to generate comparisons
corsummary.function <- function(data, qtype, raretype) {
  tmp <- data
  colnames(tmp) <- c("Group1", "Group2", "pval", "ignore")
  tmp <- separate(tmp, col=Group1, into=c("Method1", "Filt1", "Lib1"), sep="-")
  tmp <- separate(tmp, col=Group2, into=c("Method2", "Filt2", "Lib2"), sep="-")
  tmp$Hill_qval <- qtype
  tmp$Rare_type <- raretype
  tmp
}

tmp1 <- corsummary.function(u.q0.dunn.guano, qtype="q=0", raretype="unrarefied")
tmp2 <- corsummary.function(r.q0.dunn.guano, qtype="q=0", raretype="rarefied")
tmp3 <- corsummary.function(u.q1.dunn.guano, qtype="q=1", raretype="unrarefied")
tmp4 <- corsummary.function(r.q1.dunn.guano, qtype="q=1", raretype="rarefied")
tmp5 <- corsummary.function(u.q2.dunn.guano, qtype="q=2", raretype="unrarefied")
tmp6 <- corsummary.function(r.q2.dunn.guano, qtype="q=2", raretype="rarefied")
all.guano.dunn <- rbind(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6)

### How many significant differences for each Hill Number + RarefyType (note there are 630 non:self pairwise comparisons here...)
table.HillNum_by_RareType <- all.guano.dunn %>% group_by(Hill_qval, Rare_type) %>% filter(pval <= 0.05) %>% summarise(counts=n())

### How many significant differences for each Hill Number + RarefyType + FilterMethod? (there are 210 non:self pairwise comparisons per Method)
da.tmp <- all.guano.dunn %>% filter(grepl("dada2", Method1) | grepl("dada2", Method2)) %>% 
  filter(pval <= 0.05) %>%
  group_by(Hill_qval, Rare_type) %>% 
  summarise(counts=n()) %>%
  mutate(Method="dada2")

db.tmp <- all.guano.dunn %>% filter(grepl("deblur", Method1) | grepl("deblur", Method2)) %>% 
  filter(pval <= 0.05) %>%
  group_by(Hill_qval, Rare_type) %>% 
  summarise(counts=n()) %>%
  mutate(Method="deblur")

vs.tmp <- all.guano.dunn %>% filter(grepl("vsearch", Method1) | grepl("vsearch", Method2)) %>% 
  filter(pval <= 0.05) %>%
  group_by(Hill_qval, Rare_type) %>% 
  summarise(counts=n()) %>%
  mutate(Method="vsearch")

tmp.table <- rbind(da.tmp, db.tmp, vs.tmp)
table.HillNumb_by_RareType_and_Method <- spread(tmp.table, key = Hill_qval, value=counts)
