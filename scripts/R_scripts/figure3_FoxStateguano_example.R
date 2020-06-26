library(tidyverse)
library(ggpubr)

df <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/all.filtmethods.df.csv.gz")
meta <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/metadata/meta_nomock.csv") %>%
  select(-StudyID, -SampleType, -Library)
df <- merge(df, meta)
SelectMonths <- c("April", "May", "September", "October")
fox <- df %>%
  filter(Site=="FOX",
         SampleType=="sample",
         MonthStart %in% SelectMonths) %>%
  mutate(Labeler=paste(Method, Filt, MonthStart, sep="-"),
         bigID=paste(SeqID, Method, Filt, MonthStart, sep = "-"))

foxplot <- fox %>%
  group_by(SeqID, Method, Filt, MonthStart) %>%
  summarise(Richness = n_distinct(HashID))

ggplot(foxplot, aes(Richness, group=MonthStart)) +
  geom_boxplot(outlier.color = NA) +
  facet_wrap(Method ~ Filt) +
  coord_flip()

## reset axis factors for plot
foxplot$MonthStart <- factor(foxplot$MonthStart, levels = c("April", "May", "September", "October"))
foxplot$Filt <- factor(foxplot$Filt, levels = c("basic", "standard", "extra"))


my_comparisons <- list( c("April", "May"), c("April", "September"), c("April", "October"),
                        c("May", "September"), c("May", "October"),
                        c("September", "October"))

ggboxplot(foxplot, x = "MonthStart", y = "Richness",
          facet.by = c("Method", "Filt")) +
  stat_compare_means(comparisons = my_comparisons, hide.ns=TRUE) +
  stat_compare_means(label.y = 250)  +
  scale_y_continuous(breaks = c(0, 100, 200)) +
  labs(x="", y="ASVs observed")
