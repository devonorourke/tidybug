setwd("~/Repos/tidybug/data/databases/")
library(tidyverse)

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


## import data
porter_lengths <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/porter.lengths.txt.gz", delim = "\t")
porter_lengths$Dataset <- "Porter"
porter_lengths <- porter_lengths[order(porter_lengths$length),]

palmer_lengths <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/palmer.lengths.txt.gz", delim = "\t")
palmer_lengths$Dataset <- "Palmer"
palmer_lengths <- palmer_lengths[order(palmer_lengths$length),]

derep_lengths <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/derep.lengths.txt.gz", delim = "\t")
derep_lengths$Dataset <- "derep"
palmer_lengths <- palmer_lengths[order(palmer_lengths$length),]

raw_lengths <- read_delim(file = "https://github.com/devonorourke/tidybug/raw/master/data/databases/raw.lengths.txt.gz", delim = "\t", col_names = FALSE)
colnames(raw_lengths)[1] <- "length"
raw_lengths$Dataset <- "raw"


## combine
all_lengths <- rbind(porter_lengths, palmer_lengths, derep_lengths, raw_lengths)
rm(porter_lengths, palmer_lengths, derep_lengths, raw_lengths)

## summarize frequency of lengths observed, and normalize by total number of observations
obs_sumry <- all_lengths %>% 
  group_by(Dataset) %>%
  summarise(nObs=n()) 

len_sumry <- all_lengths %>% 
  group_by(Dataset, length) %>% 
  summarise(lengthCounts=n())

all_sumry <- merge(len_sumry, obs_sumry)  %>%
  mutate(pObs=lengthCounts/nObs)

rm(all_lengths, len_sumry, obs_sumry)

## write summary to disk:
write.csv(all_sumry, file = "length_summary.csv", row.names = FALSE, quote = FALSE)

## order plot levels:
all_sumry$Dataset <- factor(all_sumry$Dataset, levels=c("raw", "derep", "Palmer", "Porter"))
## color palette for plot
#notrun: vpal4 <- viridis(4, option = "plasma", end=0.85)
#notrun: pal4 <- c("#8A09A5FF", "#0D0887FF", "#DA5B6AFF", "#FEBA2CFF")

## plot; save as db_1_lengths; export at 800x400
ggplot(all_sumry, aes(x=length, y=pObs)) +
  geom_line() +
  facet_grid(Dataset ~ .) +
  scale_x_continuous(limits = c(0,2000)) +
  scale_y_continuous(limits = c(0, 0.4), breaks = c(0, 0.2, 0.4)) +
  theme_devon() +
  labs(x="length of sequence", y="fraction of sequences", title = "")
  