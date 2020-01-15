##runonce: install.packages("iNEXT")
library(iNEXT)
library(tidyverse)
library(reshape2)
library(scales)
library(ggpubr)
library(phyloseq)
## import mock data:
df <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/all.filtmethods.df.csv.gz")
mock <- df %>% filter(SampleType == "mock")

## generating two datasets: one that is rarefied and one that is not
## 1. For rarefied data, we transform each Method & Filt dataset into matrix and rarefy with Phyloseq..
## ..then run iNext specaccum to generate plot
## 2. For nonrarefied data, we do the same thing except we don't rarefy in Phyloseq

## plot function
theme_devon <- function () { 
  theme_bw(base_size=12, base_family="Avenir") %+replace% 
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill="transparent", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA)
    )
}


################################################################################
## generating data and plots for rarefied first
################################################################################

## function with rarefy
inext.function.wrare <- function(data, filter_exp, filter_exp2) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  Method <- tmp.df %>% distinct(Method)
  Filt <- tmp.df %>% distinct(Filt)
  tmp.mat <- dcast(tmp.df, HashID ~ SeqID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$HashID
  tmp.mat$HashID <- NULL
  tmp_otutbl <- otu_table(tmp.mat, taxa_are_rows = TRUE)
  tmp_rphy <- rarefy_even_depth(tmp_otutbl, sample.size = 5000, replace = FALSE, rngseed = 123)
  tmp_rphymat = as(otu_table(tmp_rphy), "matrix")
  tmp.inext <- iNEXT(tmp_rphymat, q=c(0,2), datatype = "abundance") %>% fortify(., type=1)
  tmp_plot <- data.frame(tmp.inext, Method, Filt) %>% select(-datatype, -plottype)
}

## generate data for plot
i.dada2.basic <- inext.function.wrare(mock, Method=="dada2", Filt=="basic")
i.dada2.standard <- inext.function.wrare(mock, Method=="dada2", Filt=="standard")
i.dada2.extra <- inext.function.wrare(mock, Method=="dada2", Filt=="extra")
i.deblur.basic <- inext.function.wrare(mock, Method=="deblur", Filt=="basic")
i.deblur.standard <- inext.function.wrare(mock, Method=="deblur", Filt=="standard")
i.deblur.extra <- inext.function.wrare(mock, Method=="deblur", Filt=="extra")
i.vsearch.basic <- inext.function.wrare(mock, Method=="vsearch", Filt=="basic")
i.vsearch.standard <- inext.function.wrare(mock, Method=="vsearch", Filt=="standard")
i.vsearch.extra <- inext.function.wrare(mock, Method=="vsearch", Filt=="extra")

mock.inext.all <- rbind(i.dada2.basic, i.dada2.standard, i.dada2.extra, i.deblur.basic, i.deblur.standard, i.deblur.extra, i.vsearch.basic, i.vsearch.standard, i.vsearch.extra)
rm(i.dada2.basic, i.dada2.standard, i.dada2.extra, i.deblur.basic, i.deblur.standard, i.deblur.extra, i.vsearch.basic, i.vsearch.standard, i.vsearch.extra)

## relabel mock names to match LibA-D like all other plots:
mock.inext.all$site <- as.character(mock.inext.all$site)
mock.inext.all$site[which(mock.inext.all$site=="mockIM4p4L1")] = "libA"
mock.inext.all$site[which(mock.inext.all$site=="mockIM4p4L2")] = "libB"
mock.inext.all$site[which(mock.inext.all$site=="mockIM4p7L1")] = "libC"
mock.inext.all$site[which(mock.inext.all$site=="mockIM4p7L2")] = "libD"

## reset levels for plot 
mock.inext.all$Filt <- factor(mock.inext.all$Filt,levels = c("basic", "standard", "extra"))

## color palette
libPal <- c('#CA3542', '#849FAD', '#F68930', '#BFB083')

## creating two plots; one for q=0, one for q=2

p1 <- ggplot(mock.inext.all %>% filter(order == 0), aes(x=x, y=y, color=site)) +
  geom_hline(yintercept = 24, color="gray50", size=0.75, linetype="dotted") +
  geom_line(data = mock.inext.all %>% filter(method=="interpolated" & order == 0)) +
  geom_line(data = mock.inext.all %>% filter(method=="extrapolated" & order == 0), linetype="dashed") +
  geom_ribbon(data = mock.inext.all %>% filter(order == 0), mapping=aes(x=x, ymin=y.lwr, ymax=y.upr), alpha=0.1) +
  geom_point(data = mock.inext.all %>% filter(method=="observed" & order == 0)) +
  facet_grid(Filt ~ Method) +
  scale_color_manual(values=libPal) +
  scale_x_continuous(labels=comma, breaks = c(0, 5000, 10000)) +
  scale_y_continuous(breaks = c(20, 40, 60, 80), limits = c(0,90)) +
  labs(x="sequence counts", y="sequence variants", color="") +
  theme_devon() + 
  theme(legend.position="top", axis.text.x = element_text(angle=22.5, hjust=1))

p2 <- ggplot(mock.inext.all %>% filter(order == 2), aes(x=x, y=y, color=site)) +
  geom_hline(yintercept = 24, color="gray50", size=0.75, linetype="dotted") +
  geom_line(data = mock.inext.all %>% filter(method=="interpolated" & order == 2)) +
  geom_line(data = mock.inext.all %>% filter(method=="extrapolated" & order == 2), linetype="dashed") +
  geom_ribbon(data = mock.inext.all %>% filter(order == 2), mapping=aes(x=x, ymin=y.lwr, ymax=y.upr), alpha=0.1) +
  geom_point(data = mock.inext.all %>% filter(method=="observed" & order == 2)) +
  facet_grid(Filt ~ Method) +
  scale_color_manual(values=libPal) +
  scale_x_continuous(labels=comma, breaks = c(0, 5000, 10000)) +
  scale_y_continuous(breaks = c(20, 40, 60, 80), limits = c(0,90)) +
  labs(x="sequence counts", y="sequence variant equivalents", color="") +
  theme_devon() + 
  theme(legend.position="top", axis.text.x = element_text(angle=22.5, hjust=1))

## plot; save as 4_figure_mock_specaccum_byFiltMethod; export at 1000x500
ggarrange(p1, NULL, p2, common.legend = TRUE, widths = c(1.2, 0.2, 1.2), ncol = 3, labels = c("A", "", "B"))

################################################################################
## generating data and plots for nonrarefied
################################################################################


inext.function.norare <- function(data, filter_exp, filter_exp2) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  Method <- tmp.df %>% distinct(Method)
  Filt <- tmp.df %>% distinct(Filt)
  tmp.mat <- dcast(tmp.df, HashID ~ SeqID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$HashID
  tmp.mat$HashID <- NULL
  tmp.inext <- iNEXT(tmp.mat, q=c(0,2), datatype = "abundance") %>% fortify(., type=1)
  tmp_plot <- data.frame(tmp.inext, Method, Filt) %>% select(-datatype, -plottype)
}


## generate data for plot
i.dada2.basic.norare <- inext.function.norare(mock, Method=="dada2", Filt=="basic")
i.dada2.standard.norare <- inext.function.norare(mock, Method=="dada2", Filt=="standard")
i.dada2.extra.norare <- inext.function.norare(mock, Method=="dada2", Filt=="extra")
i.deblur.basic.norare <- inext.function.norare(mock, Method=="deblur", Filt=="basic")
i.deblur.standard.norare <- inext.function.norare(mock, Method=="deblur", Filt=="standard")
i.deblur.extra.norare <- inext.function.norare(mock, Method=="deblur", Filt=="extra")
i.vsearch.basic.norare <- inext.function.norare(mock, Method=="vsearch", Filt=="basic")
i.vsearch.standard.norare <- inext.function.norare(mock, Method=="vsearch", Filt=="standard")
i.vsearch.extra.norare <- inext.function.norare(mock, Method=="vsearch", Filt=="extra")

mock.inext.all.norare <- rbind(i.dada2.basic.norare, i.dada2.standard.norare, i.dada2.extra.norare, i.deblur.basic.norare, i.deblur.standard.norare, i.deblur.extra.norare, i.vsearch.basic.norare, i.vsearch.standard.norare, i.vsearch.extra.norare)
rm(i.dada2.basic.norare, i.dada2.standard.norare, i.dada2.extra.norare, i.deblur.basic.norare, i.deblur.standard.norare, i.deblur.extra.norare, i.vsearch.basic.norare, i.vsearch.standard.norare, i.vsearch.extra.norare)

## relabel mock names to match LibA-D like all other plots:
mock.inext.all.norare$site <- as.character(mock.inext.all.norare$site)
mock.inext.all.norare$site[which(mock.inext.all.norare$site=="mockIM4p4L1")] = "libA"
mock.inext.all.norare$site[which(mock.inext.all.norare$site=="mockIM4p4L2")] = "libB"
mock.inext.all.norare$site[which(mock.inext.all.norare$site=="mockIM4p7L1")] = "libC"
mock.inext.all.norare$site[which(mock.inext.all.norare$site=="mockIM4p7L2")] = "libD"

## reset levels for plot 
mock.inext.all.norare$Filt <- factor(mock.inext.all.norare$Filt,levels = c("basic", "standard", "extra"))

p3 <- ggplot(mock.inext.all.norare %>% filter(order == 0), aes(x=x, y=y, color=site)) +
  geom_hline(yintercept = 24, color="gray50", size=0.75, linetype="dotted") +
  geom_line(data = mock.inext.all.norare %>% filter(method=="interpolated" & order == 0)) +
  geom_line(data = mock.inext.all.norare %>% filter(method=="extrapolated" & order == 0), linetype="dashed") +
  geom_ribbon(data = mock.inext.all.norare %>% filter(order == 0), mapping=aes(x=x, ymin=y.lwr, ymax=y.upr), alpha=0.1) +
  geom_point(data = mock.inext.all.norare %>% filter(method=="observed" & order == 0), size=2) +
  facet_grid(Filt ~ Method) +
  scale_color_manual(values=libPal) +
  scale_x_continuous(labels=comma, breaks = c(0, 1000000, 2000000)) +
  scale_y_continuous(breaks = c(4, 32, 256), limits = c(1, 630), trans = "log2") +
  labs(x="sequence counts", y="sequence variants", color = "") +
  theme_devon() + 
  theme(legend.position="top", axis.text.x = element_text(angle=22.5, hjust=1))

p4 <- ggplot(mock.inext.all.norare %>% filter(order == 2), aes(x=x, y=y, color=site)) +
  geom_hline(yintercept = 24, color="gray50", size=0.75, linetype="dotted") +
  geom_line(data = mock.inext.all.norare %>% filter(method=="interpolated" & order == 2)) +
  geom_line(data = mock.inext.all.norare %>% filter(method=="extrapolated" & order == 2), linetype="dashed") +
  geom_ribbon(data = mock.inext.all.norare %>% filter(order == 2), mapping=aes(x=x, ymin=y.lwr, ymax=y.upr), alpha=0.1) +
  geom_point(data = mock.inext.all.norare %>% filter(method=="observed" & order == 2), size=2) +
  facet_grid(Filt ~ Method) +
  scale_color_manual(values=libPal) +
  scale_x_continuous(labels=comma, breaks = c(0, 1000000, 2000000)) +
  scale_y_continuous(breaks = c(4, 32, 256), limits = c(1, 630), trans = "log2") +
  labs(x="sequence counts", y="sequence variant equivalents", color="") +
  theme_devon() + 
  theme(legend.position="top", axis.text.x = element_text(angle=22.5, hjust=1))

## plot; save as s4_figure_mock_specaccum_byFiltMethod_norarefy; export at 1000x500
ggarrange(p3, NULL, p4, common.legend = TRUE, widths = c(1.2, 0.2, 1.2), ncol = 3, labels = c("A", "", "B"))

################################################################################

## plot all of them?
#notrun: ggarrange(p1, p2, p3, p4, common.legend = TRUE, labels = c("A", "B", "C", "D"))