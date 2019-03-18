library(tidyverse)
library(reshape2)
library(vegan)
library(ggpubr)

## import data:
df <- read_csv("~/Repos/tidybug/data/text_tables/all.filtmethods.df.csv.gz")
meta <- read_csv("~/Repos/tidybug/data/metadata/meta_nomock.csv")
mock <- df %>% filter(SampleType == "mock")
mock$Labeler <- paste(mock$Method, mock$Filt, sep="-")
mock$bigID <- paste(mock$SeqID, mock$Method, mock$Filt, sep="-")
mock.meta <- mock %>% group_by(Method, Filt) %>% distinct(bigID)

## calculate distances for each Method and Filt group with RAREFIED data:
rrarewithdrop <- 
  function(x, sample) 
  {
    rrarefy(x[rowSums(x) >= sample, , drop = FALSE], sample)
  }

rare.betadist.function <- function(data, filter_exp, filter_exp2, betatest, binaryval) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  Method <- tmp.df %>% distinct(Method)
  Filt <- tmp.df %>% distinct(Filt)
  tmp.mat <- dcast(tmp.df, SeqID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SeqID
  SeqID <- row.names(tmp.mat)
  tmp.mat$SeqID <- NULL
  tmp.mat2 <- data.frame(rrarewithdrop(tmp.mat, 5000))
  tmp.betadiv <- vegdist(tmp.mat2, betatest, binary = binaryval)
  tmp.dist <- data.frame(t(combn(rownames(tmp.mat2),2)), as.numeric(tmp.betadiv))  ## pairwise distances
  colnames(tmp.dist) <- c("SeqID1", "SeqID2", "betadistance")
  data.frame(tmp.dist, Method, Filt, betatest, binaryval)
}


betadist.bray.dada2.basic <- rare.betadist.function(mock, Method=="dada2", Filt=="basic", "bray", FALSE)
betadist.dice.dada2.basic <- rare.betadist.function(mock, Method=="dada2", Filt=="basic", "bray", TRUE)
betadist.mori.dada2.basic <- rare.betadist.function(mock, Method=="dada2", Filt=="basic", "morisita", FALSE)
betadist.bray.dada2.standard <- rare.betadist.function(mock, Method=="dada2", Filt=="standard", "bray", FALSE)
betadist.dice.dada2.standard <- rare.betadist.function(mock, Method=="dada2", Filt=="standard", "bray", TRUE)
betadist.mori.dada2.standard <- rare.betadist.function(mock, Method=="dada2", Filt=="standard", "morisita", FALSE)
betadist.bray.dada2.extra <- rare.betadist.function(mock, Method=="dada2", Filt=="extra", "bray", FALSE)
betadist.dice.dada2.extra <- rare.betadist.function(mock, Method=="dada2", Filt=="extra", "bray", TRUE)
betadist.mori.dada2.extra <- rare.betadist.function(mock, Method=="dada2", Filt=="extra", "morisita", FALSE)

betadist.bray.deblur.basic <- rare.betadist.function(mock, Method=="deblur", Filt=="basic", "bray", FALSE)
betadist.dice.deblur.basic <- rare.betadist.function(mock, Method=="deblur", Filt=="basic", "bray", TRUE)
betadist.mori.deblur.basic <- rare.betadist.function(mock, Method=="deblur", Filt=="basic", "morisita", FALSE)
betadist.bray.deblur.standard <- rare.betadist.function(mock, Method=="deblur", Filt=="standard", "bray", FALSE)
betadist.dice.deblur.standard <- rare.betadist.function(mock, Method=="deblur", Filt=="standard", "bray", TRUE)
betadist.mori.deblur.standard <- rare.betadist.function(mock, Method=="deblur", Filt=="standard", "morisita", FALSE)
betadist.bray.deblur.extra <- rare.betadist.function(mock, Method=="deblur", Filt=="extra", "bray", FALSE)
betadist.dice.deblur.extra <- rare.betadist.function(mock, Method=="deblur", Filt=="extra", "bray", TRUE)
betadist.mori.deblur.extra <- rare.betadist.function(mock, Method=="deblur", Filt=="extra", "morisita", FALSE)

betadist.bray.vsearch.basic <- rare.betadist.function(mock, Method=="vsearch", Filt=="basic", "bray", FALSE)
betadist.dice.vsearch.basic <- rare.betadist.function(mock, Method=="vsearch", Filt=="basic", "bray", TRUE)
betadist.mori.vsearch.basic <- rare.betadist.function(mock, Method=="vsearch", Filt=="basic", "morisita", FALSE)
betadist.bray.vsearch.standard <- rare.betadist.function(mock, Method=="vsearch", Filt=="standard", "bray", FALSE)
betadist.dice.vsearch.standard <- rare.betadist.function(mock, Method=="vsearch", Filt=="standard", "bray", TRUE)
betadist.mori.vsearch.standard <- rare.betadist.function(mock, Method=="vsearch", Filt=="standard", "morisita", FALSE)
betadist.bray.vsearch.extra <- rare.betadist.function(mock, Method=="vsearch", Filt=="extra", "bray", FALSE)
betadist.dice.vsearch.extra <- rare.betadist.function(mock, Method=="vsearch", Filt=="extra", "bray", TRUE)
betadist.mori.vsearch.extra <- rare.betadist.function(mock, Method=="vsearch", Filt=="extra", "morisita", FALSE)

all.betadist.df <- rbind(betadist.dice.vsearch.extra, betadist.dice.vsearch.standard, betadist.dice.vsearch.basic,
                         betadist.bray.vsearch.extra, betadist.bray.vsearch.standard, betadist.bray.vsearch.basic,
                         betadist.mori.vsearch.extra, betadist.mori.vsearch.standard, betadist.mori.vsearch.basic,
                         betadist.dice.dada2.extra, betadist.dice.dada2.standard, betadist.dice.dada2.basic,
                         betadist.bray.dada2.extra, betadist.bray.dada2.standard, betadist.bray.dada2.basic,
                         betadist.mori.dada2.extra, betadist.mori.dada2.standard, betadist.mori.dada2.basic,
                         betadist.dice.deblur.extra, betadist.dice.deblur.standard, betadist.dice.deblur.basic,
                         betadist.bray.deblur.extra, betadist.bray.deblur.standard, betadist.bray.deblur.basic,
                         betadist.mori.deblur.extra, betadist.mori.deblur.standard, betadist.mori.deblur.basic)

row.names(all.betadist.df) <- NULL

## replace factor with character type
all.betadist.df$betatest <- as.character(all.betadist.df$betatest)

# replace $betatest "bray" that is "dice"
all.betadist.df$betatest <- ifelse(grepl(TRUE, all.betadist.df$binaryval), gsub("bray", "dice", all.betadist.df$betatest), all.betadist.df$betatest)

## cleanup:
rm(list=ls(pattern = "betadist.dice*"))
rm(list=ls(pattern = "betadist.bray*"))
rm(list=ls(pattern = "betadist.mori*"))

#write.csv(all.betadist.df, file="~/Repos/tidybug/data/text_tables/betadist_MOCK.out.csv", quote=FALSE, row.names = FALSE)

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
#all.betadist.df$Labeler <- paste(all.betadist.df$betatest, all.betadist.df$Filt, all.betadist.df$MonthComp, sep="-")
all.betadist.df$LegendLabel <- paste(all.betadist.df$betatest, all.betadist.df$Filt, sep="\n")
## levels for plots
all.betadist.df$Filt <- factor(all.betadist.df$Filt, levels = c("basic", "standard", "extra"))
all.betadist.df$LegendLabel <- factor(all.betadist.df$LegendLabel,
                                      levels=c("dice\nbasic", "dice\nstandard", "dice\nextra",
                                               "bray\nbasic",  "bray\nstandard","bray\nextra",
                                               "morisita\nbasic", "morisita\nstandard", "morisita\nextra"))

## color and shape sets:
pal9 <- c(rep('#9f9244', 3), rep('#6c42b8', 3), rep('#628a47',3))
shape9 <- rep(c(0,1,2),3)

## plot; save as 13_mock_betadist_boxplot; export at 900x800
ggplot(all.betadist.df, aes(x=LegendLabel, y=betadistance, color=LegendLabel, shape=LegendLabel)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.5, position=position_jitterdodge(jitter.width = 0.3)) +
  theme_devon() +
  scale_color_manual(values=pal9) +
  scale_shape_manual(values=shape9) +
  facet_grid(Method ~ .) +
  labs(title="", x="", y="distance", color="", shape="") +
  theme(legend.position="top", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  guides(color=guide_legend(nrow = 1, title.position = "top"))
