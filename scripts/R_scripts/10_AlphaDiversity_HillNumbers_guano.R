library(tidyverse)
library(reshape2)
library(vegan)
library(cowplot)

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
df <- df %>% filter(SampleType != "mock")

## Hull Numbers function applied to calculate diversity measures: observed OTUs, simpson, and shannon
hillnumber.function <- function(data, filter_exp, filter_exp2) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  Method <- tmp.df %>% distinct(Method)
  Filt <- tmp.df %>% distinct(Filt)
  tmp.mat <- dcast(tmp.df, SeqID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SeqID
  SeqID <- row.names(tmp.mat)
  tmp.mat$SeqID <- NULL
  tmp.hill <- renyi(tmp.mat, scales = c(0,1,2), hill=TRUE)
  tmp <- data.frame(tmp.hill, SeqID, Method, Filt, row.names = NULL)
  colnames(tmp)[1:3] <- c("q=0", "q=1", "q=2")
  gather(tmp, key="Hill_qType", value = "Hill_value", c('q=0', 'q=1', 'q=2'))
}

## calculate Hill Numbers for unrarefied data
dada2.basic <- hillnumber.function(df, Method=="dada2", Filt=="basic")
dada2.standard <- hillnumber.function(df, Method=="dada2", Filt=="standard")
dada2.extra <- hillnumber.function(df, Method=="dada2", Filt=="extra")
deblur.basic <- hillnumber.function(df, Method=="deblur", Filt=="basic")
deblur.standard <- hillnumber.function(df, Method=="deblur", Filt=="standard")
deblur.extra <- hillnumber.function(df, Method=="deblur", Filt=="extra")
vsearch.basic <- hillnumber.function(df, Method=="vsearch", Filt=="basic")
vsearch.standard <- hillnumber.function(df, Method=="vsearch", Filt=="standard")
vsearch.extra <- hillnumber.function(df, Method=="vsearch", Filt=="extra")

## merge into single dataframe
all.df.hill <- rbind(dada2.basic, dada2.standard, dada2.extra, deblur.basic, deblur.standard, deblur.extra, vsearch.basic, vsearch.standard, vsearch.extra)
rm(dada2.basic, dada2.standard, dada2.extra, deblur.basic, deblur.standard, deblur.extra, vsearch.basic, vsearch.standard, vsearch.extra)
#rm(df)
## set the levels
all.df.hill$Filt <- factor(all.df.hill$Filt, levels = c("basic", "standard", "extra"))
## add labeler to make facet labels easier to understand
all.df.hill$Labeler <- paste(all.df.hill$Method, all.df.hill$Hill_qType, sep="  ")

## split df by Hill q-value:
q0.hill.df <- all.df.hill %>% filter(Hill_qType == "q=0")
q1.hill.df <- all.df.hill %>% filter(Hill_qType == "q=1")
q2.hill.df <- all.df.hill %>% filter(Hill_qType == "q=2")

## generate palette for 6 colors following plot5 color scheme but altering hue:
pal3 <- c('#9f9244', '#6c42b8', '#628a47', '#a9d190')

## make separate plots per qtype, then merge into single figure
plot.q0 <- ggplot(q0.hill.df, aes(x=Filt, y=Hill_value, color=Filt, group=Filt, shape=Filt)) +
  scale_color_manual(values = pal3) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.2, position=position_jitterdodge(jitter.width = 0.3)) +
  scale_y_continuous(limits = c(0,400)) +
  facet_grid( ~ Labeler) +
  labs(title = "", x="", y="sequence variants", color="", shape="") +
  theme_devon() + theme(legend.position = "top", axis.text.x = element_blank(), axis.ticks.x = element_blank())

plot.q1 <- ggplot(q1.hill.df, aes(x=Filt, y=Hill_value, color=Filt, group=Filt, shape=Filt)) +
  scale_color_manual(values = pal3) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.2, position=position_jitterdodge(jitter.width = 0.3)) +
  facet_grid( ~ Labeler) +
  labs(x="", y="sequence variant equivalents", color="", shape="") +
  theme_devon() + theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())

plot.q2 <- ggplot(q2.hill.df, aes(x=Filt, y=Hill_value, color=Filt, group=Filt, shape=Filt)) +
  scale_color_manual(values = pal3) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.2, position=position_jitterdodge(jitter.width = 0.3)) +
  facet_grid( ~ Labeler) +
  labs(x="", y="sequence variant equivalents", color="", shape="") +
  theme_devon() + theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())

## plot; save as 10_guanoAlpha_HillNumbers; export at 1000x1000
plot_grid(plot.q0, plot.q1, plot.q2, ncol=1, rel_heights = c(1.3, 1, 1))
