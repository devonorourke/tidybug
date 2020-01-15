library(tidyverse)
library(reshape2)
library(vegan)
library(phyloseq)
library(dunn.test)

## import data:
df <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/all.filtmethods.df.csv.gz")
mock <- df %>% filter(SampleType=="mock")
rm(df)

## rename mock samples 
mock$SeqID <- as.character(mock$SeqID)
mock$SeqID[which(mock$SeqID=="mockIM4p4L1")] <- "libA"
mock$SeqID[which(mock$SeqID=="mockIM4p4L2")] <- "libB"
mock$SeqID[which(mock$SeqID=="mockIM4p7L1")] <- "libC"
mock$SeqID[which(mock$SeqID=="mockIM4p7L2")] <- "libD"



betadist.function <- function(data, filter_exp, filter_exp2, betatest, binaryval) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  Method <- tmp.df %>% distinct(Method)
  Filt <- tmp.df %>% distinct(Filt)
  BetaType <- binaryval
  tmp.meta <- tmp.df %>% distinct(SeqID)
  tmp.mat <- dcast(tmp.df, SeqID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SeqID
  SeqID <- row.names(tmp.mat)
  tmp.mat$SeqID <- NULL
  tmp.betadiv <- vegdist(tmp.mat, betatest, binary = binaryval)
  tmp.dist <- data.frame(t(combn(rownames(tmp.mat),2)), as.numeric(tmp.betadiv))  ## pairwise distances
  colnames(tmp.dist) <- c("SeqID1", "SeqID2", "betadistance")
  tmp.out <- data.frame(tmp.dist, Method, Filt, BetaType)
  tmp.out <- merge(tmp.out, tmp.meta, by.x = "SeqID1", by.y = "SeqID")
}


betadist.morisita.dada2.basic <- betadist.function(mock, Method=="dada2", Filt=="basic", "morisita", FALSE)
betadist.morisita.dada2.standard <- betadist.function(mock, Method=="dada2", Filt=="standard", "morisita", FALSE)
betadist.morisita.dada2.extra <- betadist.function(mock, Method=="dada2", Filt=="extra", "morisita", FALSE)
betadist.morisita.deblur.basic <- betadist.function(mock, Method=="deblur", Filt=="basic", "morisita", FALSE)
betadist.morisita.deblur.standard <- betadist.function(mock, Method=="deblur", Filt=="standard", "morisita", FALSE)
betadist.morisita.deblur.extra <- betadist.function(mock, Method=="deblur", Filt=="extra", "morisita", FALSE)
betadist.morisita.vsearch.basic <- betadist.function(mock, Method=="vsearch", Filt=="basic", "morisita", FALSE)
betadist.morisita.vsearch.standard <- betadist.function(mock, Method=="vsearch", Filt=="standard", "morisita", FALSE)
betadist.morisita.vsearch.extra <- betadist.function(mock, Method=="vsearch", Filt=="extra", "morisita", FALSE)

tmp1 <- rbind(betadist.morisita.dada2.basic, betadist.morisita.dada2.standard, betadist.morisita.dada2.extra, betadist.morisita.deblur.basic, betadist.morisita.deblur.standard, betadist.morisita.deblur.extra, betadist.morisita.vsearch.basic, betadist.morisita.vsearch.standard, betadist.morisita.vsearch.extra)
tmp1$BetaType[which(tmp1$BetaType==FALSE)] <- "morisita"
rm(list=ls(pattern = "betadist.morisita*"))

betadist.bray.dada2.basic <- betadist.function(mock, Method=="dada2", Filt=="basic", "bray", FALSE)
betadist.dice.dada2.basic <- betadist.function(mock, Method=="dada2", Filt=="basic", "bray", TRUE)
betadist.bray.dada2.standard <- betadist.function(mock, Method=="dada2", Filt=="standard", "bray", FALSE)
betadist.dice.dada2.standard <- betadist.function(mock, Method=="dada2", Filt=="standard", "bray", TRUE)
betadist.bray.dada2.extra <- betadist.function(mock, Method=="dada2", Filt=="extra", "bray", FALSE)
betadist.dice.dada2.extra <- betadist.function(mock, Method=="dada2", Filt=="extra", "bray", TRUE)
betadist.bray.deblur.basic <- betadist.function(mock, Method=="deblur", Filt=="basic", "bray", FALSE)
betadist.dice.deblur.basic <- betadist.function(mock, Method=="deblur", Filt=="basic", "bray", TRUE)
betadist.bray.deblur.standard <- betadist.function(mock, Method=="deblur", Filt=="standard", "bray", FALSE)
betadist.dice.deblur.standard <- betadist.function(mock, Method=="deblur", Filt=="standard", "bray", TRUE)
betadist.bray.deblur.extra <- betadist.function(mock, Method=="deblur", Filt=="extra", "bray", FALSE)
betadist.dice.deblur.extra <- betadist.function(mock, Method=="deblur", Filt=="extra", "bray", TRUE)
betadist.bray.vsearch.basic <- betadist.function(mock, Method=="vsearch", Filt=="basic", "bray", FALSE)
betadist.dice.vsearch.basic <- betadist.function(mock, Method=="vsearch", Filt=="basic", "bray", TRUE)
betadist.bray.vsearch.standard <- betadist.function(mock, Method=="vsearch", Filt=="standard", "bray", FALSE)
betadist.dice.vsearch.standard <- betadist.function(mock, Method=="vsearch", Filt=="standard", "bray", TRUE)
betadist.bray.vsearch.extra <- betadist.function(mock, Method=="vsearch", Filt=="extra", "bray", FALSE)
betadist.dice.vsearch.extra <- betadist.function(mock, Method=="vsearch", Filt=="extra", "bray", TRUE)

tmp2 <- rbind(betadist.dice.vsearch.extra, betadist.dice.vsearch.standard, betadist.dice.vsearch.basic,
                         betadist.bray.vsearch.extra, betadist.bray.vsearch.standard, betadist.bray.vsearch.basic,
                         betadist.dice.dada2.extra, betadist.dice.dada2.standard, betadist.dice.dada2.basic,
                         betadist.bray.dada2.extra, betadist.bray.dada2.standard, betadist.bray.dada2.basic,
                         betadist.dice.deblur.extra, betadist.dice.deblur.standard, betadist.dice.deblur.basic,
                         betadist.bray.deblur.extra, betadist.bray.deblur.standard, betadist.bray.deblur.basic)
rm(list=ls(pattern = "betadist.dice*"))
rm(list=ls(pattern = "betadist.bray*"))
tmp2$BetaType[which(tmp2$BetaType==FALSE)] <- "bray"
tmp2$BetaType[which(tmp2$BetaType==TRUE)] <- "dice"
tmp3 <- rbind(tmp1, tmp2)
rm(tmp1, tmp2)
tmp3$RareType <- "unrarefied"


########### and for rarefied data
rare.betadist.function <- function(data, filter_exp, filter_exp2, betatest, binaryval) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  Method <- tmp.df %>% distinct(Method)
  Filt <- tmp.df %>% distinct(Filt)
  BetaType <- binaryval
  tmp.meta <- tmp.df %>% distinct(SeqID)
  tmp.mat <- dcast(tmp.df, SeqID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SeqID
  SeqID <- row.names(tmp.mat)
  tmp.mat$SeqID <- NULL
  tmp.mat2 <- rrarefy(tmp.mat, 5000)
  tmp.betadiv <- vegdist(tmp.mat2, betatest, binary = binaryval)
  tmp.dist <- data.frame(t(combn(rownames(tmp.mat2),2)), as.numeric(tmp.betadiv))  ## pairwise distances
  colnames(tmp.dist) <- c("SeqID1", "SeqID2", "betadistance")
  tmp.out <- data.frame(tmp.dist, Method, Filt, BetaType)
  tmp.out <- merge(tmp.out, tmp.meta, by.x = "SeqID1", by.y = "SeqID")
}

betadist.morisita.dada2.basic <- rare.betadist.function(mock, Method=="dada2", Filt=="basic", "morisita", FALSE)
betadist.morisita.dada2.standard <- rare.betadist.function(mock, Method=="dada2", Filt=="standard", "morisita", FALSE)
betadist.morisita.dada2.extra <- rare.betadist.function(mock, Method=="dada2", Filt=="extra", "morisita", FALSE)
betadist.morisita.deblur.basic <- rare.betadist.function(mock, Method=="deblur", Filt=="basic", "morisita", FALSE)
betadist.morisita.deblur.standard <- rare.betadist.function(mock, Method=="deblur", Filt=="standard", "morisita", FALSE)
betadist.morisita.deblur.extra <- rare.betadist.function(mock, Method=="deblur", Filt=="extra", "morisita", FALSE)
betadist.morisita.vsearch.basic <- rare.betadist.function(mock, Method=="vsearch", Filt=="basic", "morisita", FALSE)
betadist.morisita.vsearch.standard <- rare.betadist.function(mock, Method=="vsearch", Filt=="standard", "morisita", FALSE)
betadist.morisita.vsearch.extra <- rare.betadist.function(mock, Method=="vsearch", Filt=="extra", "morisita", FALSE)

tmp1 <- rbind(betadist.morisita.dada2.basic, betadist.morisita.dada2.standard, betadist.morisita.dada2.extra, betadist.morisita.deblur.basic, betadist.morisita.deblur.standard, betadist.morisita.deblur.extra, betadist.morisita.vsearch.basic, betadist.morisita.vsearch.standard, betadist.morisita.vsearch.extra)
tmp1$BetaType[which(tmp1$BetaType==FALSE)] <- "morisita"
rm(list=ls(pattern = "betadist.morisita*"))

betadist.bray.dada2.basic <- rare.betadist.function(mock, Method=="dada2", Filt=="basic", "bray", FALSE)
betadist.dice.dada2.basic <- rare.betadist.function(mock, Method=="dada2", Filt=="basic", "bray", TRUE)
betadist.bray.dada2.standard <- rare.betadist.function(mock, Method=="dada2", Filt=="standard", "bray", FALSE)
betadist.dice.dada2.standard <- rare.betadist.function(mock, Method=="dada2", Filt=="standard", "bray", TRUE)
betadist.bray.dada2.extra <- rare.betadist.function(mock, Method=="dada2", Filt=="extra", "bray", FALSE)
betadist.dice.dada2.extra <- rare.betadist.function(mock, Method=="dada2", Filt=="extra", "bray", TRUE)
betadist.bray.deblur.basic <- rare.betadist.function(mock, Method=="deblur", Filt=="basic", "bray", FALSE)
betadist.dice.deblur.basic <- rare.betadist.function(mock, Method=="deblur", Filt=="basic", "bray", TRUE)
betadist.bray.deblur.standard <- rare.betadist.function(mock, Method=="deblur", Filt=="standard", "bray", FALSE)
betadist.dice.deblur.standard <- rare.betadist.function(mock, Method=="deblur", Filt=="standard", "bray", TRUE)
betadist.bray.deblur.extra <- rare.betadist.function(mock, Method=="deblur", Filt=="extra", "bray", FALSE)
betadist.dice.deblur.extra <- rare.betadist.function(mock, Method=="deblur", Filt=="extra", "bray", TRUE)
betadist.bray.vsearch.basic <- rare.betadist.function(mock, Method=="vsearch", Filt=="basic", "bray", FALSE)
betadist.dice.vsearch.basic <- rare.betadist.function(mock, Method=="vsearch", Filt=="basic", "bray", TRUE)
betadist.bray.vsearch.standard <- rare.betadist.function(mock, Method=="vsearch", Filt=="standard", "bray", FALSE)
betadist.dice.vsearch.standard <- rare.betadist.function(mock, Method=="vsearch", Filt=="standard", "bray", TRUE)
betadist.bray.vsearch.extra <- rare.betadist.function(mock, Method=="vsearch", Filt=="extra", "bray", FALSE)
betadist.dice.vsearch.extra <- rare.betadist.function(mock, Method=="vsearch", Filt=="extra", "bray", TRUE)

tmp2 <- rbind(betadist.dice.vsearch.extra, betadist.dice.vsearch.standard, betadist.dice.vsearch.basic,
              betadist.bray.vsearch.extra, betadist.bray.vsearch.standard, betadist.bray.vsearch.basic,
              betadist.dice.dada2.extra, betadist.dice.dada2.standard, betadist.dice.dada2.basic,
              betadist.bray.dada2.extra, betadist.bray.dada2.standard, betadist.bray.dada2.basic,
              betadist.dice.deblur.extra, betadist.dice.deblur.standard, betadist.dice.deblur.basic,
              betadist.bray.deblur.extra, betadist.bray.deblur.standard, betadist.bray.deblur.basic)
rm(list=ls(pattern = "betadist.dice*"))
rm(list=ls(pattern = "betadist.bray*"))
tmp2$BetaType[which(tmp2$BetaType==FALSE)] <- "bray"
tmp2$BetaType[which(tmp2$BetaType==TRUE)] <- "dice"
tmp4 <- rbind(tmp1, tmp2)
rm(tmp1, tmp2)
tmp4$RareType <- "rarefied"

## merge rarefied and unrarefied data
all.betadist.df <- rbind(tmp3, tmp4)
rm(tmp3, tmp4)

## theme for plot
theme_devon <- function () { 
  theme_bw(base_size=12, base_family="Avenir") %+replace% 
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill="transparent", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA)
    )
}

## labeler for plots
all.betadist.df$LegendLabel <- paste(all.betadist.df$BetaType, all.betadist.df$Filt, sep=" - ")
## levels for plots
all.betadist.df$RareType <- factor(all.betadist.df$RareType, levels = c("unrarefied", "rarefied"))
all.betadist.df$Filt <- factor(all.betadist.df$Filt, levels = c("basic", "standard", "extra"))
all.betadist.df$LegendLabel <- factor(all.betadist.df$LegendLabel,
                                      levels=c("dice - basic", "dice - standard", "dice - extra",
                                               "bray - basic", "bray - standard", "bray - extra",
                                               "morisita - basic", "morisita - standard", "morisita - extra"))

## extra label for geom_text
#notrun: all.betadist.df$pwiseID <- paste(all.betadist.df$SeqID1, all.betadist.df$SeqID2, sep = ":")

## get max value for plotted texts
#notrun: max.dice.standard <- all.betadist.df %>% filter(LegendLabel=="dice - standard") %>% summarise(betadistance=max(betadistance))
#notrun: max.dice.standard <- as.character(max.dice.standard)
#notrun: max.bray.standard <- all.betadist.df %>% filter(LegendLabel=="bray - standard") %>% summarise(betadistance=max(betadistance))
#notrun: max.bray.standard <- as.character(max.bray.standard)
#notrun: max.morisita.standard <- all.betadist.df %>% filter(LegendLabel=="morisita - standard") %>% summarise(betadistance=max(betadistance))
#notrun: max.morisita.standard <- as.character(max.morisita.standard)


## color and shape sets:
pal3 <- rep(c('#9f9244', '#6c42b8', '#628a47'),4)
shape4 <- c(rep(15,3),rep(0,3),rep(2,3),rep(1,3))


## plot; save as 13_mock_betadist_boxplot; export at 800x800
#ggplot(all.betadist.df, aes(x=LegendLabel, y=betadistance, color=LegendLabel, shape=LegendLabel, label=pwiseID)) +
ggplot(all.betadist.df, aes(x=LegendLabel, y=betadistance, color=LegendLabel, shape=LegendLabel)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.7, position=position_jitterdodge(jitter.width = 0.2)) +
  #geom_label_repel(data=all.betadist.df %>% filter(LegendLabel=="bray - standard" & betadistance == max.bray.standard), nudge_x = 0.8) +
  #geom_label_repel(data=all.betadist.df %>% filter(LegendLabel=="dice - standard" & betadistance == max.dice.standard), nudge_x = 0.8) +
  #geom_label_repel(data=all.betadist.df %>% filter(LegendLabel=="morisita - standard" & betadistance == max.morisita.standard), nudge_y = 0.2, nudge_x = 0.2) +
  facet_grid( Method ~ RareType) +
  scale_color_manual(values=pal3) +
  scale_shape_manual(values=shape4) +
  labs(title="", x="", y="distance", color="", shape="") +
  theme_devon() +
  theme(legend.position="top", axis.text.x = element_text(angle=22.5, hjust=1, size=9)) +
  guides(col=guide_legend(nrow = 3))

# write to disk:
write.csv(all.betadist.df, file="~/Repos/tidybug/data/text_tables/mock_betadist.out.csv", quote=FALSE, row.names = FALSE)


## ANOVA estimate of variation
tmp <- aov(betadistance ~ Method * Filt * RareType * BetaType, data = all.betadist.df)
summary(tmp)

## Kruskal-Wallis
all.betadist.df$Labeler <- paste(all.betadist.df$Method, all.betadist.df$Filt, all.betadist.df$BetaType, all.betadist.df$RareType, sep="-")
kruskal.test(betadistance ~ Labeler, data=all.betadist.df)
all.betadist.df$Labeler <- as.factor(all.betadist.df$Labeler)

## Dunn's test
tmp.dunn = dunnTest(betadistance ~ Labeler, data=all.betadist.df, method="bh")
mock.dunn.df <- (tmp.dunn$res)


