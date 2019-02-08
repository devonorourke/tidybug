# Plots used in paper

library(tidyverse)
library(scales)
library(ggridges)

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
pipepal3 <- c('#9f9244', '#6c42b8', '#628a47', '#a9d190')

## generate summaries for number of ASVs per sample, per pipeline and filtering method 
df_hashCounts <- df %>% group_by(Method, Filt, SeqID, Library) %>% summarise(HashCounts=n())

## reset the levels for plot
df_hashCounts$Filt <- factor(df_hashCounts$Filt,levels = c("basic", "standard", "extra"))
df_hashCounts$Method <- factor(df_hashCounts$Method,levels = c("dada2", "deblur", "vsearch"))

## plot; save as 6_figure_GuanoASVCounts_byFilterMethod
ggplot(data = df_hashCounts, aes(x = Method, y = HashCounts, color=Method)) +
  geom_hline(yintercept = 50, linetype="dotted", color="firebrick", size=1) +
  geom_jitter(alpha=0.55, width = 0.25) +
  facet_grid(Filt ~ Library) +
  scale_y_continuous(trans = "log2") +
  scale_color_manual(values=pipepal3) +
  labs(title="", x="", y="sequence variants per sample", color="") +
  theme_devon() +
  theme(legend.position="top", legend.text = element_text(size = 12), axis.text.x = element_text(angle=22.5, hjust = 1))

