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

## reset the levels for plot
df$Filt <- factor(df$Filt,levels = c("basic", "standard", "extra"))
df$Method <- factor(df$Method, levels = c("dada2", "deblur", "vsearch"))

## plot; save as 5_figure_GuanoSeqCounts_byFilterMethod
ggplot(data = df, aes(y = Method, x = Reads, label = Alias, fill=Method)) +
  geom_density_ridges(scale=3, alpha=0.7) +
  scale_x_continuous(labels = comma, trans = "log2") +
  facet_grid(Filt ~ Library) +
  scale_fill_manual(values=pipepal3) +
  labs(title="", x="sequence counts", y="", fill="") +
  theme_devon() +
  theme(legend.position="top", legend.text = element_text(size = 12), axis.text.x = element_text(angle=22.5, hjust = 1))

        