## script generates Figure3 in paper

library(tidyverse)
library(scales)
library(ggridges)
library(ggpubr)

# create theme function for all plots
# theme function for custom plot style:
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
df <- df %>% filter(SampleType != "mock")

## generate 3 color palette to distinguish between filtering pipelines:
pal3 <- c('#9f9244', '#6c42b8', '#628a47', '#a9d190')

## reset the levels for plot
df$Filt <- factor(df$Filt,levels = c("basic", "standard", "extra"))
df$Method <- factor(df$Method, levels = c("dada2", "deblur", "vsearch"))

## generate 3 plots to combine
ridge.ba <- ggplot(data = df %>% filter(Filt=="basic"), aes(y = Method, x = Reads, fill=Method)) +
  geom_density_ridges(scale=5, alpha=0.7) +
  scale_x_continuous(trans = "log2",
                     limits = c(1,1000),
                     breaks = c(1, 4, 16, 64, 256)) +
  facet_wrap(~ Filt) +
  scale_fill_manual(values=pal3) +
  labs(x="", y="", fill="") +
  theme_devon() +
  theme(legend.position = "none", 
        axis.text.x = element_blank(), axis.ticks.x = element_blank())


ridge.st <- ggplot(data = df %>% filter(Filt=="standard"), aes(y = Method, x = Reads, fill=Method)) +
  geom_density_ridges(scale=5, alpha=0.7) +
  scale_x_continuous(trans = "log2",
                     limits = c(1,1000),
                     breaks = c(1, 4, 16, 64, 256)) +
  facet_wrap(~ Filt) +
  scale_fill_manual(values=pal3) +
  labs(x="", y="Fraction of observations", fill="") +
  theme_devon() +
  theme(legend.position = "none", 
        axis.text.x = element_blank(), axis.ticks.x = element_blank())


ridge.ex <- ggplot(data = df %>% filter(Filt=="extra"), aes(y = Method, x = Reads, fill=Method)) +
  geom_density_ridges(scale=5, alpha=0.7) +
  scale_x_continuous(trans = "log2",
                     limits = c(1,1000),
                     breaks = c(1, 4, 16, 64, 256)) +
  facet_wrap(~ Filt) +
  scale_fill_manual(values=pal3) +
  labs(x="Number of sequences per ASV", y="", fill="") +
  theme_devon() + 
  theme(legend.position = "none")

## save as 3_figure_GuanoSeqCounts_byFilterMethod; export at 600x600
ggarrange(ridge.ba, ridge.st, ridge.ex, nrow=3)


################################################################################
## Alternatively, just use jitter plot and histograms to show number of reads across all guano data
################################################################################

## save as s2_figure_GuanoSeqCounts_perLibrary_byFilterMethod_Violin; export at 750x750
ggplot(df, aes(x=Library, y=Reads, fill=Method)) +
  #geom_jitter(alpha = 0.05) +
  geom_violin() +
  scale_fill_manual(values=pal3) +
  facet_grid(Filt ~ Method) +
  theme_devon() +
  theme(legend.position = "none") +
  scale_y_continuous(trans="log2", labels = comma) +
  labs(x="", y="sequence counts")

## save as 3_figure_GuanoSeqCounts_comboLibs_byFilterMethod_Violin; export at 666x222
ggplot(df, aes(x=Method, y=Reads, fill=Method)) +
  #geom_jitter(alpha = 0.01) +
  geom_violin() +
  scale_fill_manual(values=pal3) +
  #geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ Filt) +
  theme_devon() + 
  theme(legend.position = "none") +
  scale_y_continuous(trans="log2", labels = comma) +
  labs(x = "", y="sequence counts")


################################################################################
## can also illustrate the per-ASV scatterplot of read depth and frequency of occurrence
################################################################################

tmp <- df %>% 
  group_by(Method, Filt, HashID) %>% 
  summarise(nSamples = n_distinct(SeqID),
            nReads = sum(Reads))

da_basic <- df %>% filter(Method == "dada2" & Filt == "basic")

tmp$Filt <- factor(tmp$Filt, levels=c("basic", "standard", "extra"))

## save as s3_figure_GuanoSeqAndOccurrenceCounts_comboLibs_scatterplot; export at 750x750
ggplot(tmp, aes(x=nSamples, y=nReads)) + 
  geom_point(alpha = 0.25) + 
  facet_grid(Filt ~ Method) +
  theme_devon() +
  scale_y_continuous(labels=comma, breaks = c(0, 1000000, 2000000)) +
  labs(x="Samples", y="Sequence counts")
