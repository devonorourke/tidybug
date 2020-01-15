##runonce: install.packages("iNEXT")
library(iNEXT)
library(tidyverse)
library(reshape2)
library(scales)
## import mock data:
df <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/all.filtmethods.df.csv.gz")
mock <- df %>% filter(SampleType == "mock")
rm(df)

## inext function
inext.function <- function(data, filter_exp, filter_exp2) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  Method <- tmp.df %>% distinct(Method)
  Filt <- tmp.df %>% distinct(Filt)
  tmp.mat <- dcast(tmp.df, HashID ~ SeqID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$HashID
  tmp.mat$HashID <- NULL
  tmp.inext <- iNEXT(tmp.mat, q=c(0), datatype = "abundance")
  df <- fortify(tmp.inext, type=1)
  data.frame(df, Method, Filt) %>% select(-datatype, -plottype, -order)
}

## generate data for plot
i.dada2.basic <- inext.function(mock, Method=="dada2", Filt=="basic")
i.dada2.standard <- inext.function(mock, Method=="dada2", Filt=="standard")
i.dada2.extra <- inext.function(mock, Method=="dada2", Filt=="extra")
i.deblur.basic <- inext.function(mock, Method=="deblur", Filt=="basic")
i.deblur.standard <- inext.function(mock, Method=="deblur", Filt=="standard")
i.deblur.extra <- inext.function(mock, Method=="deblur", Filt=="extra")
i.vsearch.basic <- inext.function(mock, Method=="vsearch", Filt=="basic")
i.vsearch.standard <- inext.function(mock, Method=="vsearch", Filt=="standard")
i.vsearch.extra <- inext.function(mock, Method=="vsearch", Filt=="extra")

## combine:
mock.inext.all <- rbind(i.dada2.basic, i.dada2.standard, i.dada2.extra, i.deblur.basic, i.deblur.standard, i.deblur.extra, i.vsearch.basic, i.vsearch.standard, i.vsearch.extra)
rm(i.dada2.basic, i.dada2.standard, i.dada2.extra, i.deblur.basic, i.deblur.standard, i.deblur.extra, i.vsearch.basic, i.vsearch.standard, i.vsearch.extra)


ggplot(mock.inext.all, aes(x=x, y=y, color=Filt, shape=Filt)) +
  geom_point(data = mock.inext.all %>% filter(method=="observed")) +
  geom_line(data = mock.inext.all %>% filter(method=="interpolated")) +
  geom_line(data = mock.inext.all %>% filter(method=="extrapolated"), linetype="dashed") +
  facet_grid(Method ~ site) +
  ylim(0,150) + 
  geom_ribbon(data = mock.inext.all, mapping=aes(x=x, ymin=y.lwr, ymax=y.upr, fill=Filt), alpha=0.1) +
  scale_x_continuous(labels=comma, trans = "log2") +
  theme_devon()

  head(mock.inext.all)
