library(tidyverse)
library(reshape2)
library(vegan)
library(ggpubr)

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
  BetaType <- binaryval
  tmp.meta1 <- tmp.df %>% distinct(SeqID, MonthStart, WOY)
  tmp.meta2 <- tmp.df %>% distinct(SeqID, MonthStart)
  colnames(tmp.meta2)[2] <- "MonthComp"
  tmp.mat <- dcast(tmp.df, SeqID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SeqID
  SeqID <- row.names(tmp.mat)
  tmp.mat$SeqID <- NULL
  tmp.mat2 <- rrarewithdrop(tmp.mat, 5000)
  tmp.betadiv <- vegdist(tmp.mat2, betatest, binary = binaryval)
  tmp.dist <- data.frame(t(combn(rownames(tmp.mat2),2)), as.numeric(tmp.betadiv))  ## pairwise distances
  colnames(tmp.dist) <- c("SeqID1", "SeqID2", "betadistance")
  tmp.dist2 <- tmp.dist
  colnames(tmp.dist2) <- c("SeqID2", "SeqID1", "betadistance")
  tmp.dist <- rbind(tmp.dist, tmp.dist2)
  tmp.out <- data.frame(tmp.dist, Method, Filt, BetaType)
  tmp.out <- merge(tmp.out, tmp.meta1, by.x = "SeqID1", by.y = "SeqID")
  merge(tmp.out, tmp.meta2, by.x="SeqID2", by.y="SeqID")
}


betadist.bray.dada2.basic <- rare.betadist.function(fox, Method=="dada2", Filt=="basic", "bray", FALSE)
betadist.dice.dada2.basic <- rare.betadist.function(fox, Method=="dada2", Filt=="basic", "bray", TRUE)
betadist.bray.dada2.standard <- rare.betadist.function(fox, Method=="dada2", Filt=="standard", "bray", FALSE)
betadist.dice.dada2.standard <- rare.betadist.function(fox, Method=="dada2", Filt=="standard", "bray", TRUE)
betadist.bray.dada2.extra <- rare.betadist.function(fox, Method=="dada2", Filt=="extra", "bray", FALSE)
betadist.dice.dada2.extra <- rare.betadist.function(fox, Method=="dada2", Filt=="extra", "bray", TRUE)

betadist.bray.deblur.basic <- rare.betadist.function(fox, Method=="deblur", Filt=="basic", "bray", FALSE)
betadist.dice.deblur.basic <- rare.betadist.function(fox, Method=="deblur", Filt=="basic", "bray", TRUE)
betadist.bray.deblur.standard <- rare.betadist.function(fox, Method=="deblur", Filt=="standard", "bray", FALSE)
betadist.dice.deblur.standard <- rare.betadist.function(fox, Method=="deblur", Filt=="standard", "bray", TRUE)
betadist.bray.deblur.extra <- rare.betadist.function(fox, Method=="deblur", Filt=="extra", "bray", FALSE)
betadist.dice.deblur.extra <- rare.betadist.function(fox, Method=="deblur", Filt=="extra", "bray", TRUE)

betadist.bray.vsearch.basic <- rare.betadist.function(fox, Method=="vsearch", Filt=="basic", "bray", FALSE)
betadist.dice.vsearch.basic <- rare.betadist.function(fox, Method=="vsearch", Filt=="basic", "bray", TRUE)
betadist.bray.vsearch.standard <- rare.betadist.function(fox, Method=="vsearch", Filt=="standard", "bray", FALSE)
betadist.dice.vsearch.standard <- rare.betadist.function(fox, Method=="vsearch", Filt=="standard", "bray", TRUE)
betadist.bray.vsearch.extra <- rare.betadist.function(fox, Method=="vsearch", Filt=="extra", "bray", FALSE)
betadist.dice.vsearch.extra <- rare.betadist.function(fox, Method=="vsearch", Filt=="extra", "bray", TRUE)

all.betadist.df <- rbind(betadist.dice.vsearch.extra, betadist.dice.vsearch.standard, betadist.dice.vsearch.basic,
                       betadist.bray.vsearch.extra, betadist.bray.vsearch.standard, betadist.bray.vsearch.basic,
                       betadist.dice.dada2.extra, betadist.dice.dada2.standard, betadist.dice.dada2.basic,
                       betadist.bray.dada2.extra, betadist.bray.dada2.standard, betadist.bray.dada2.basic,
                       betadist.dice.deblur.extra, betadist.dice.deblur.standard, betadist.dice.deblur.basic,
                       betadist.bray.deblur.extra, betadist.bray.deblur.standard, betadist.bray.deblur.basic)
row.names(all.betadist.df) <- NULL
all.betadist.df$BetaType[which(all.betadist.df$BetaType==TRUE)] <- "dice"
all.betadist.df$BetaType[which(all.betadist.df$BetaType==FALSE)] <- "bray"

rm(list=ls(pattern = "betadist.dice*"))
rm(list=ls(pattern = "betadist.bray*"))

#write.csv(all.betadist.df, file="~/Repos/tidybug/data/text_tables/betadist.out.csv", quote=FALSE, row.names = FALSE)

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
#all.betadist.df$Labeler <- paste(all.betadist.df$BetaType, all.betadist.df$Filt, all.betadist.df$MonthComp, sep="-")
all.betadist.df$LegendLabel <- paste(all.betadist.df$BetaType, all.betadist.df$Filt, sep=" - ")
## levels for plots
all.betadist.df$Filt <- factor(all.betadist.df$Filt, levels = c("basic", "standard", "extra"))
all.betadist.df$MonthStart <- factor(all.betadist.df$MonthStart, levels = c("April", "May", "September", "October"))
all.betadist.df$MonthComp <- factor(all.betadist.df$MonthComp, levels = c("April", "May", "September", "October"))

all.betadist.df$LegendLabel <- factor(all.betadist.df$LegendLabel,
                                      levels=c("dice - basic","bray - basic","dice - standard",
                                               "bray - standard","dice - extra","bray - extra"))


## color and shape sets:
pal6 <- c('#9f9244', '#ebdb8e', '#6c42b8', '#c8b2e8', '#628a47', '#a9d190')
shape6 <- c(0,0,1,1,2,2)

## plot
april <- ggplot(all.betadist.df %>% filter(MonthStart=="April"), aes(x=LegendLabel, y=betadistance, color=LegendLabel, shape=LegendLabel)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.15, position=position_jitterdodge(jitter.width = 0.3)) +
  facet_grid(Method ~ MonthStart + MonthComp) +
  theme_devon() +
  scale_color_manual(values=pal6) +
  scale_shape_manual(values=shape6) +
  labs(title="", x="", y="distance", color="", shape="") +
  theme(legend.position="top", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  guides(color=guide_legend(nrow = 1, title.position = "top"))

may <- ggplot(all.betadist.df %>% filter(MonthStart=="May"), aes(x=LegendLabel, y=betadistance, color=LegendLabel, shape=LegendLabel)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.15, position=position_jitterdodge(jitter.width = 0.3)) +
  facet_grid(Method ~ MonthStart + MonthComp) +
  theme_devon() +
  scale_color_manual(values=pal6) +
  scale_shape_manual(values=shape6) +
  labs(title="", x="", y="distance", color="", shape="") +
  theme(legend.position="top", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  guides(color=guide_legend(nrow = 1, title.position = "top"))


sept <- ggplot(all.betadist.df %>% filter(MonthStart=="September"), aes(x=LegendLabel, y=betadistance, color=LegendLabel, shape=LegendLabel)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.15, position=position_jitterdodge(jitter.width = 0.3)) +
  facet_grid(Method ~ MonthStart + MonthComp) +
  theme_devon() +
  scale_color_manual(values=pal6) +
  scale_shape_manual(values=shape6) +
  labs(title="", x="", y="distance", color="", shape="") +
  theme(legend.position="top", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  guides(color=guide_legend(nrow = 1, title.position = "top"))

oct <- ggplot(all.betadist.df %>% filter(MonthStart=="October"), aes(x=LegendLabel, y=betadistance, color=LegendLabel, shape=LegendLabel)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.15, position=position_jitterdodge(jitter.width = 0.3)) +
  facet_grid(Method ~ MonthStart + MonthComp) +
  theme_devon() +
  scale_color_manual(values=pal6) +
  scale_shape_manual(values=shape6) +
  labs(title="", x="", y="distance", color="", shape="") +
  theme(legend.position="top", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  guides(color=guide_legend(nrow = 1, title.position = "top"))
  
##plot; save as 1_rarefied_betadist_boxplot; export at 1000x1000
ggarrange(april, may, sept, oct, common.legend = TRUE)
