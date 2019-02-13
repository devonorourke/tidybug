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
rm(SelectMonths)

## function
betadisp.function <- function(data, filter_exp, filter_exp2, betatest, binaryval) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  Method <- tmp.df %>% distinct(Method)
  Filt <- tmp.df %>% distinct(Filt)
  BetaType <- binaryval
  tmp.meta <- tmp.df %>% distinct(SeqID, MonthStart, WOY)
  groups <- as.character(tmp.meta$MonthStart)
  tmp.mat <- dcast(tmp.df, SeqID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SeqID
  SeqID <- row.names(tmp.mat)
  tmp.mat$SeqID <- NULL
  tmp.betadiv <- vegdist(tmp.mat, betatest, binary = binaryval)
  tmp.disp <- betadisper(tmp.betadiv, groups, type = c("median"))
  tmp.disp_sco <- as.data.frame(scores(tmp.disp, display = "sites"))
  tmp.disp_sco$SeqID <- row.names(tmp.disp_sco)
  row.names(tmp.disp_sco) <- NULL
  colnames(tmp.disp_sco) <- c("pc1end", "pc2end", "SeqID")
  tmp.disp_sco$Type <- "sample"
  tmp.disp_sco <- merge(tmp.disp_sco, tmp.meta)
  tmp.disp_centroids <- as.data.frame(scores(tmp.disp, display="centroids"))
  tmp.disp_centroids$SeqID <- row.names(tmp.disp_centroids)
  row.names(tmp.disp_centroids) <- NULL
  colnames(tmp.disp_centroids) <- c("pc1start", "pc2start", "MonthStart")
  tmp.disp.plot <- merge(tmp.disp_sco, tmp.disp_centroids)
  data.frame(tmp.disp.plot, Method, Filt, BetaType)
}

## datasets
bray.basic.dada2 <- betadisp.function(fox, Method=="dada2", Filt=="basic", "bray", TRUE)
dice.basic.dada2 <- betadisp.function(fox, Method=="dada2", Filt=="basic", "bray", FALSE)
bray.standard.dada2 <- betadisp.function(fox, Method=="dada2", Filt=="standard", "bray", TRUE)
dice.standard.dada2 <- betadisp.function(fox, Method=="dada2", Filt=="standard", "bray", FALSE)
bray.extra.dada2 <- betadisp.function(fox, Method=="dada2", Filt=="extra", "bray", TRUE)
dice.extra.dada2 <- betadisp.function(fox, Method=="dada2", Filt=="extra", "bray", FALSE)
bray.basic.deblur <- betadisp.function(fox, Method=="deblur", Filt=="basic", "bray", TRUE)
dice.basic.deblur <- betadisp.function(fox, Method=="deblur", Filt=="basic", "bray", FALSE)
bray.standard.deblur <- betadisp.function(fox, Method=="deblur", Filt=="standard", "bray", TRUE)
dice.standard.deblur <- betadisp.function(fox, Method=="deblur", Filt=="standard", "bray", FALSE)
bray.extra.deblur <- betadisp.function(fox, Method=="deblur", Filt=="extra", "bray", TRUE)
dice.extra.deblur <- betadisp.function(fox, Method=="deblur", Filt=="extra", "bray", FALSE)
bray.basic.vsearch <- betadisp.function(fox, Method=="vsearch", Filt=="basic", "bray", TRUE)
dice.basic.vsearch <- betadisp.function(fox, Method=="vsearch", Filt=="basic", "bray", FALSE)
bray.standard.vsearch <- betadisp.function(fox, Method=="vsearch", Filt=="standard", "bray", TRUE)
dice.standard.vsearch <- betadisp.function(fox, Method=="vsearch", Filt=="standard", "bray", FALSE)
bray.extra.vsearch <- betadisp.function(fox, Method=="vsearch", Filt=="extra", "bray", TRUE)
dice.extra.vsearch <- betadisp.function(fox, Method=="vsearch", Filt=="extra", "bray", FALSE)

betadisp.all.df <- rbind(bray.basic.dada2, dice.basic.dada2, bray.standard.dada2, dice.standard.dada2, bray.extra.dada2, dice.extra.dada2,
                     bray.basic.deblur, dice.basic.deblur, bray.standard.deblur, dice.standard.deblur, bray.extra.deblur, dice.extra.deblur,
                     bray.basic.vsearch, dice.basic.vsearch, bray.standard.vsearch, dice.standard.vsearch, bray.extra.vsearch, dice.extra.vsearch)
rm(bray.basic.dada2, dice.basic.dada2, bray.standard.dada2, dice.standard.dada2, bray.extra.dada2, dice.extra.dada2,
   bray.basic.deblur, dice.basic.deblur, bray.standard.deblur, dice.standard.deblur, bray.extra.deblur, dice.extra.deblur,
   bray.basic.vsearch, dice.basic.vsearch, bray.standard.vsearch, dice.standard.vsearch, bray.extra.vsearch, dice.extra.vsearch)
betadisp.all.df$BetaType[which(betadisp.all.df$BetaType==TRUE)] <- "bray"
betadisp.all.df$BetaType[which(betadisp.all.df$BetaType==FALSE)] <- "dice"

## plotting:
## define color palette:
pal4 <- c("#85c069","#a386d6","#d59247","#eda4ba")

## levels
betadisp.all.df$Filt <- factor(betadisp.all.df$Filt, levels = c("basic", "standard", "extra"))
betadisp.all.df$BetaType <- factor(betadisp.all.df$BetaType, levels = c("dice", "bray"))

## plot dispersions
#ggplot(dicedisp.deblur.extra.dfplot, aes(x=pc1end, y=pc2end, shape=MonthStart, color=MonthStart)) + 
ggplot(betadisp.all.df, aes(x=pc1end, y=pc2end, shape=MonthStart, color=MonthStart, label=WOY)) + 
  geom_point(data=betadisp.all.df, aes(x=pc1end, y=pc2end), size=0.1) +
  geom_text(data=betadisp.all.df, aes(x=pc1end, y=pc2end), size=3, show.legend = FALSE) +
  geom_segment(data=betadisp.all.df, aes(x=pc1end, y=pc2end, xend=pc1start, yend=pc2start), alpha=0.5) +
  scale_color_manual(values=pal4) + 
  facet_grid(Filt ~ Method + BetaType) +
  labs(title="", color="", x="PCoA1", y="PCoA2", label="", shape="") +
  theme_devon() +
  theme(legend.position = "top") +
  guides(shape = guide_legend(override.aes = list(size=3)))


## Unused...
anova(tmp.disp)
## basic plot
plot(tmp.disp)
## no labels, but added ellipses
plot(tmp.disp, label = T, hull = T, ellipse = T)
## boxplot
boxplot(tmp.disp)
## significane tests
TukeyHSD(tmp.disp)
