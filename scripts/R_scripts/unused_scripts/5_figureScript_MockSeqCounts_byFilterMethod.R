# Plots used in paper

library(tidyverse)
library(scales)
library(ggridges)
library(cowplot)

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
df <- df %>% filter(SampleType == "mock")

## generate 3 color palette to distinguish between filtering pipelines:
pal3 <- c('#9f9244', '#6c42b8', '#628a47', '#a9d190')

## reset the levels for plot
df$Filt <- factor(df$Filt,levels = c("basic", "standard", "extra"))
df$Method <- factor(df$Method, levels = c("dada2", "deblur", "vsearch"))

## generate plots to stitch together
ridge.da <- ggplot(data = df %>% filter(Method=="dada2"), aes(y = Filt, x = Reads, fill=Filt)) +
  geom_density_ridges(scale=1, alpha=0.7,
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = '|', point_size = 3, point_alpha = 1) +
  #scale_x_continuous(labels = comma, trans = "log2") +
  #xlim(0,100) +
  facet_grid(Method ~ Library, scales = "free_x") +
  scale_fill_manual(values=pal3) +
  labs(title="", x="", y="", fill="") +
  theme_devon() +
  theme(legend.position="top", axis.text.x = element_text(angle=22.5, hjust = 1),
        legend.text = element_text(size = 12))


ridge.db <- ggplot(data = df %>% filter(Method=="deblur"), aes(y = Filt, x = Reads, fill=Filt)) +
  geom_density_ridges(scale=1, alpha=0.7,
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = '|', point_size = 3, point_alpha = 1) +
  #scale_x_continuous(labels = comma, trans = "log2") +
  #xlim(0,100) +
  facet_grid(Method ~ Library, scales = "free_x") +
  scale_fill_manual(values=pal3) +
  labs(x="", y="frequency of observed sequence count") +
  theme_devon() +
  theme(legend.position="none", axis.text.x = element_text(angle=22.5, hjust = 1),
        strip.background.x = element_blank(), strip.text.x = element_blank())


ridge.vs <- ggplot(data = df %>% filter(Method=="vsearch"), aes(y = Filt, x = Reads, fill=Filt)) +
  geom_density_ridges(scale=1, alpha=0.7,
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = '|', point_size = 3, point_alpha = 1) +
  #scale_x_continuous(labels = comma, trans = "log2") +
  #xlim(0,100) +
  facet_grid(Method ~ Library, scales = "free_x") +
  scale_fill_manual(values=pal3) +
  labs(x="sequence counts", y="", fill="") +
  theme_devon() +
  theme(legend.position="none", axis.text.x = element_text(angle=22.5, hjust = 1),
        strip.background.x = element_blank(), strip.text.x = element_blank())

## plot; save as 5_figure_MockSeqCounts_byFilterMethod
plot_grid(ridge.da, ridge.db, ridge.vs, ncol=1, rel_heights = c(1.7,1.1,1.3))
