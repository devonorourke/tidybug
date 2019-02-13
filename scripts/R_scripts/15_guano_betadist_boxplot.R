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
  
}


betadist.bray.dada2.basic <- rare.betadist.function(fox, Method=="dada2", Filt=="basic", "bray", FALSE)
betadist.dice.dada2.basic <- rare.betadist.function(fox, Method=="dada2", Filt=="basic", "bray", TRUE)
betadist.bray.dada2.standard <- rare.betadist.function(fox, Method=="dada2", Filt=="standard", "bray", FALSE)
betadist.dice.dada2.standard <- rare.betadist.function(fox, Method=="dada2", Filt=="standard", "bray", TRUE)
betadist.bray.dada2.extra <- rare.betadist.function(fox, Method=="dada2", Filt=="extra", "bray", FALSE)
betadist.dice.dada2.extra <- rare.betadist.function(fox, Method=="dada2", Filt=="extra", "bray", TRUE)

tmp1 <- rbind(betadist.bray.dada2.basic, betadist.dice.dada2.basic, betadist.bray.dada2.standard, 
              betadist.dice.dada2.standard, betadist.bray.dada2.extra, betadist.dice.dada2.extra)
row.names(tmp1) <- NULL
tmp1$BetaType[which(tmp1$BetaType==TRUE)] <- "dice"
tmp1$BetaType[which(tmp1$BetaType==FALSE)] <- "bray"

rm(list=ls(pattern = "betadist.dice*"))
rm(list=ls(pattern = "betadist.bray*"))

betadist.morisita.dada2.basic <- rare.betadist.function(fox, Method=="dada2", Filt=="basic", "morisita", FALSE)
betadist.morisita.dada2.standard <- rare.betadist.function(fox, Method=="dada2", Filt=="standard", "morisita", FALSE)
betadist.morisita.dada2.extra <- rare.betadist.function(fox, Method=="dada2", Filt=="extra", "morisita", FALSE)

tmp2 <- rbind(betadist.morisita.dada2.basic, betadist.morisita.dada2.standard, betadist.morisita.dada2.extra)
tmp2$BetaType[which(tmp2$BetaType==FALSE)] <- "morisita"
rm(list=ls(pattern = "betadist.morisita*"))

all.betadist.df <- rbind(tmp1, tmp2)

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
all.betadist.df$LegendLabel <- paste(all.betadist.df$BetaType, all.betadist.df$Filt, sep=" - ")

## levels for plots
all.betadist.df$Filt <- factor(all.betadist.df$Filt, levels = c("basic", "standard", "extra"))
all.betadist.df$MonthStart <- factor(all.betadist.df$MonthStart, levels = c("April", "May", "September", "October"))
all.betadist.df$MonthComp <- factor(all.betadist.df$MonthComp, levels = c("April", "May", "September", "October"))
all.betadist.df$LegendLabel <- factor(all.betadist.df$LegendLabel,
                                      levels=c("dice - basic","bray - basic","morisita - basic",
                                               "dice - standard","bray - standard","morisita - standard",
                                               "dice - extra","bray - extra","morisita - extra"))

## color and shape sets:
pals6 <- c('gold', 'yellow', 'purple', 'pink', 'green', 'lightgreen')
pal6 <- c('#9f9244', '#ebdb8e', '#6c42b8', '#c8b2e8', '#628a47', '#a9d190')
shape9 <- c(15,15,15,0,0,0,2,2,2)

## plot
april <- ggplot(all.betadist.df %>% filter(MonthStart=="April"), aes(x=LegendLabel, y=betadistance, color=LegendLabel, shape=LegendLabel)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.15, position=position_jitterdodge(jitter.width = 0.3)) +
  facet_grid(MonthStart ~ MonthComp) +
  theme_devon() +
  #scale_color_manual(values=pal6) +
  scale_shape_manual(values=shape9) +
  labs(title="", x="", y="distance", color="", shape="") +
  theme(legend.position="top", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  guides(color=guide_legend(nrow = 1, title.position = "top"))

may <- ggplot(all.betadist.df %>% filter(MonthStart=="May"), aes(x=LegendLabel, y=betadistance, color=LegendLabel, shape=LegendLabel)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.15, position=position_jitterdodge(jitter.width = 0.3)) +
  facet_grid(MonthStart ~ MonthComp) +
  theme_devon() +
  #scale_color_manual(values=pal6) +
  scale_shape_manual(values=shape9) +
  labs(title="", x="", y="distance", color="", shape="") +
  theme(legend.position="top", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  guides(color=guide_legend(nrow = 1, title.position = "top"))


sept <- ggplot(all.betadist.df %>% filter(MonthStart=="September"), aes(x=LegendLabel, y=betadistance, color=LegendLabel, shape=LegendLabel)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.15, position=position_jitterdodge(jitter.width = 0.3)) +
  facet_grid(MonthStart ~ MonthComp) +
  theme_devon() +
  #scale_color_manual(values=pal6) +
  scale_shape_manual(values=shape9) +
  labs(title="", x="", y="distance", color="", shape="") +
  theme(legend.position="top", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  guides(color=guide_legend(nrow = 1, title.position = "top"))

oct <- ggplot(all.betadist.df %>% filter(MonthStart=="October"), aes(x=LegendLabel, y=betadistance, color=LegendLabel, shape=LegendLabel)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.15, position=position_jitterdodge(jitter.width = 0.3)) +
  facet_grid(MonthStart ~ MonthComp) +
  theme_devon() +
  #scale_color_manual(values=pal6) +
  scale_shape_manual(values=shape9) +
  labs(title="", x="", y="distance", color="", shape="") +
  theme(legend.position="top", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  guides(color=guide_legend(nrow = 1, title.position = "top"))
  
##plot; save as 1_rarefied_betadist_boxplot; export at 1000x1000
ggarrange(april, may, sept, oct, common.legend = TRUE)
