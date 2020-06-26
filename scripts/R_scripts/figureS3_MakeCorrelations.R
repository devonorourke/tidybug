library(tidyverse)
library(ggpubr)

df <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/all.filtmethods.df.csv.gz")

## gather summary info of nASVs among different denoising/filtering pipelines
## note we're retaining only those that have common values across all pipelines (754 samples)
sumrydf <- df %>%
  filter(SampleType == "sample") %>%
  group_by(Method, Filt, SeqID) %>%
  summarise(nASVs=n()) %>%
  mutate(Grouper = paste(Method, Filt, sep="_")) %>%
  ungroup() %>%
  select(-Method, -Filt) %>%
  pivot_wider(names_from = Grouper,
              values_from = nASVs) %>%
  filter(complete.cases(.))

## further filter data into three groups for plotting:
basicdat <- sumrydf %>% select(-ends_with("standard"), -ends_with("extra"))
standdat <- sumrydf %>% select(-ends_with("basic"), -ends_with("extra"))
extradat <- sumrydf %>% select(-ends_with("basic"), -ends_with("standard"))

rm(df, sumrydf)


## plot each subset of data, then stitch together
bp1 <- ggscatter(basicdat, x = "dada2_basic", y = "deblur_basic", add = "reg.line",
                 color = "black", shape = 21, size = 0.75) +
  stat_cor(label.x = 3, label.y = 525) +
  stat_regline_equation(label.x = 3, label.y = 575) +
  labs(x="ASVs dada2-basic", y="ASVs deblur-basic") +
  theme(axis.text = element_text(size=10), axis.title.x = element_text(color="gold")) +
  scale_x_continuous(limits = c(0,615)) + scale_y_continuous(limits = c(0,615))

bp2 <- ggscatter(basicdat, x = "dada2_basic", y = "vsearch_basic", add = "reg.line",
                 color = "black", shape = 21, size = 0.75) +
  labs(x="ASVs dada2-basic", y="ASVs vsearch-basic") +
  stat_cor(label.x = 3, label.y = 525) +
  stat_regline_equation(label.x = 3, label.y = 575) +
  theme(axis.text = element_text(size=10)) + scale_x_continuous(limits = c(0,615)) + scale_y_continuous(limits = c(0,615))
  
  bp3 <- ggscatter(basicdat, x = "deblur_basic", y = "vsearch_basic", add = "reg.line",
                 color = "black", shape = 21, size = 0.75) +
  labs(x="ASVs deblur-basic", y="ASVs vsearch-basic") +
  stat_cor(label.x = 3, label.y = 525) +
  stat_regline_equation(label.x = 3, label.y = 575) +
  theme(axis.text = element_text(size=10)) + scale_x_continuous(limits = c(0,615)) + scale_y_continuous(limits = c(0,615))

sp1 <- ggscatter(standdat, x = "dada2_standard", y = "deblur_standard", add = "reg.line",
                 color = "black", shape = 21, size = 0.75) +
  stat_cor(label.x = 3, label.y = 450) +
  stat_regline_equation(label.x = 3, label.y = 490) +
  labs(x="ASVs dada2-standard", y="ASVs deblur-standard") +
  theme(axis.text = element_text(size=10)) + scale_x_continuous(limits = c(0,615)) + scale_y_continuous(limits = c(0,615))

sp2 <- ggscatter(standdat, x = "dada2_standard", y = "vsearch_standard", add = "reg.line",
                 color = "black", shape = 21, size = 0.75) +
  labs(x="ASVs dada2-standard", y="ASVs vsearch-standard") +
  stat_cor(label.x = 3, label.y = 450) +
  stat_regline_equation(label.x = 3, label.y = 490) +
  theme(axis.text = element_text(size=10)) + scale_x_continuous(limits = c(0,615)) + scale_y_continuous(limits = c(0,615))

sp3 <- ggscatter(standdat, x = "deblur_standard", y = "vsearch_standard", add = "reg.line",
                 color = "black", shape = 21, size = 0.75) +
  labs(x="ASVs deblur-standard", y="ASVs vsearch-standard") +
  stat_cor(label.x = 3, label.y = 450) +
  stat_regline_equation(label.x = 3, label.y = 490) +
  theme(axis.text = element_text(size=10)) + scale_x_continuous(limits = c(0,615)) + scale_y_continuous(limits = c(0,615))

ep1 <- ggscatter(extradat, x = "dada2_extra", y = "deblur_extra", add = "reg.line",
                 color = "black", shape = 21, size = 0.75) +
  stat_cor(label.x = 3, label.y = 175) +
  stat_regline_equation(label.x = 3, label.y = 190) +
  labs(x="ASVs dada2-extra", y="ASVs deblur-extra") +
  theme(axis.text = element_text(size=10)) + scale_x_continuous(limits = c(0,200)) + scale_y_continuous(limits = c(0,200))

ep2 <- ggscatter(extradat, x = "dada2_extra", y = "vsearch_extra", add = "reg.line",
                 color = "black", shape = 21, size = 0.75) +
  labs(x="ASVs dada2-extra", y="ASVs vsearch-extra") +
  stat_cor(label.x = 3, label.y = 175) +
  stat_regline_equation(label.x = 3, label.y = 190) +
  theme(axis.text = element_text(size=10)) + scale_x_continuous(limits = c(0,200)) + scale_y_continuous(limits = c(0,200))
  
  ep3 <- ggscatter(extradat, x = "deblur_extra", y = "vsearch_extra", add = "reg.line",
                 color = "black", shape = 21, size = 0.75) +
  labs(x="ASVs deblur-extra", y="ASVs vsearch-extra") +
  stat_cor(label.x = 3, label.y = 175) +
  stat_regline_equation(label.x = 3, label.y = 190) +
  theme(axis.text = element_text(size=10)) + scale_x_continuous(limits = c(0,200)) + scale_y_continuous(limits = c(0,200))

ggarrange(bp1, bp2, bp3, sp1, sp2, sp3, ep1, ep2, ep3, ncol = 3, nrow=3)
#ggsave("asvRichness_correlation.png", width = 15, height=15, units = "cm")
