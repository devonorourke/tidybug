## runonce: install.packages("goeveg")
library(tidyverse)
library(reshape2)
library(vegan)
library(goeveg)

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
meta <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/metadata/meta_nomock.csv/")
meta <- meta %>% select(-StudyID, -SampleType, -Library)
df <- merge(df, meta)
rm(meta)
fox <- df %>% filter(Site=="FOX")
rm(df)
## only selected months
SelectMonths <- c("April", "May", "September", "October")
fox <- fox %>% filter(MonthStart %in% SelectMonths)
rm(SelectMonths)


## rarefy function, per Method, per Filt
rrarewithdrop <- function(x, sample) {
  rrarefy(x[rowSums(x) >= sample, , drop = FALSE], sample)
}

rarefy.function <- function(data, filter_exp, filter_exp2) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  tmp.mat <- dcast(tmp.df, SeqID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SeqID
  tmp.mat$SeqID <- NULL
  tmp.mat2 <- data.frame(rrarewithdrop(tmp.mat, 5000))
}


## screeplot function, adapted from goeveg
screeplot.function <- function (matrix, distance = "bray", k = 6, trymax = 20, autotransform = FALSE, binary = FALSE, Method, Filt, BetaType) 
{
  stress <- 0
  for (i in 1:k) {
    nmds_i <- metaMDS(matrix, distance = distance, k = i, 
                      trymax = trymax, engine = "monoMDS", autotransform = autotransform, binary = binary)
    stress[i] <- nmds_i$stress
  }
  stress <- print(stress)
  dims <- seq(1:k)
  df <- data.frame(dims, stress)
  df <- df %>% mutate(Method = Method) %>% mutate(Filt = Filt) %>% mutate(BetaType = BetaType)
}


## gather data for scree plot
mat.dada2.basic <- rarefy.function(fox, Method=="dada2", Filt=="basic")
mat.dada2.standard <- rarefy.function(fox, Method=="dada2", Filt=="standard")
mat.dada2.extra <- rarefy.function(fox, Method=="dada2", Filt=="extra")
mat.deblur.basic <- rarefy.function(fox, Method=="deblur", Filt=="basic")
mat.deblur.standard <- rarefy.function(fox, Method=="deblur", Filt=="standard")
mat.deblur.extra <- rarefy.function(fox, Method=="deblur", Filt=="extra")
mat.vsearch.basic <- rarefy.function(fox, Method=="vsearch", Filt=="basic")
mat.vsearch.standard <- rarefy.function(fox, Method=="vsearch", Filt=="standard")
mat.vsearch.extra <- rarefy.function(fox, Method=="vsearch", Filt=="extra")

## screeplot data
scree.dada2.basic.dice <- screeplot.function(mat.dada2.basic, binary = TRUE, Method = "dada2", Filt = "basic", BetaType="dice")
scree.dada2.basic.bray <- screeplot.function(mat.dada2.basic, Method = "dada2", Filt = "basic", BetaType="bray")
scree.dada2.basic.mori <- screeplot.function(mat.dada2.basic, distance = "morisita", Method = "dada2", Filt = "basic", BetaType="mori")
scree.dada2.standard.dice <- screeplot.function(mat.dada2.standard, binary = TRUE, Method = "dada2", Filt = "standard", BetaType="dice")
scree.dada2.standard.bray <- screeplot.function(mat.dada2.standard, Method = "dada2", Filt = "standard", BetaType="bray")
scree.dada2.standard.mori <- screeplot.function(mat.dada2.standard, distance = "morisita", Method = "dada2", Filt = "standard", BetaType="mori")
scree.dada2.extra.dice <- screeplot.function(mat.dada2.extra, binary = TRUE, Method = "dada2", Filt = "extra", BetaType="dice")
scree.dada2.extra.bray <- screeplot.function(mat.dada2.extra, Method = "dada2", Filt = "extra", BetaType="bray")
scree.dada2.extra.mori <- screeplot.function(mat.dada2.extra, distance = "morisita", Method = "dada2", Filt = "extra", BetaType="mori")

scree.deblur.basic.dice <- screeplot.function(mat.deblur.basic, binary = TRUE, Method = "deblur", Filt = "basic", BetaType="dice")
scree.deblur.basic.bray <- screeplot.function(mat.deblur.basic, Method = "deblur", Filt = "basic", BetaType="bray")
scree.deblur.basic.mori <- screeplot.function(mat.deblur.basic, distance = "morisita", Method = "deblur", Filt = "basic", BetaType="mori")
scree.deblur.standard.dice <- screeplot.function(mat.deblur.standard, binary = TRUE, Method = "deblur", Filt = "standard", BetaType="dice")
scree.deblur.standard.bray <- screeplot.function(mat.deblur.standard, Method = "deblur", Filt = "standard", BetaType="bray")
scree.deblur.standard.mori <- screeplot.function(mat.deblur.standard, distance = "morisita", Method = "deblur", Filt = "standard", BetaType="mori")
scree.deblur.extra.dice <- screeplot.function(mat.deblur.extra, binary = TRUE, Method = "deblur", Filt = "extra", BetaType="dice")
scree.deblur.extra.bray <- screeplot.function(mat.deblur.extra, Method = "deblur", Filt = "extra", BetaType="bray")
scree.deblur.extra.mori <- screeplot.function(mat.deblur.extra, distance = "morisita", Method = "deblur", Filt = "extra", BetaType="mori")

scree.vsearch.basic.dice <- screeplot.function(mat.vsearch.basic, binary = TRUE, Method = "vsearch", Filt = "basic", BetaType="dice")
scree.vsearch.basic.bray <- screeplot.function(mat.vsearch.basic, Method = "vsearch", Filt = "basic", BetaType="bray")
scree.vsearch.basic.mori <- screeplot.function(mat.vsearch.basic, distance = "morisita", Method = "vsearch", Filt = "basic", BetaType="mori")
scree.vsearch.standard.dice <- screeplot.function(mat.vsearch.standard, binary = TRUE, Method = "vsearch", Filt = "standard", BetaType="dice")
scree.vsearch.standard.bray <- screeplot.function(mat.vsearch.standard, Method = "vsearch", Filt = "standard", BetaType="bray")
scree.vsearch.standard.mori <- screeplot.function(mat.vsearch.standard, distance = "morisita", Method = "vsearch", Filt = "standard", BetaType="mori")
scree.vsearch.extra.dice <- screeplot.function(mat.vsearch.extra, binary = TRUE, Method = "vsearch", Filt = "extra", BetaType="dice")
scree.vsearch.extra.bray <- screeplot.function(mat.vsearch.extra, Method = "vsearch", Filt = "extra", BetaType="bray")
scree.vsearch.extra.mori <- screeplot.function(mat.vsearch.extra, distance = "morisita", Method = "vsearch", Filt = "extra", BetaType="mori")

## combine all scree plot data
scree_df <- do.call("rbind", mget(ls(pattern = "scree\\..*")))

## set levels
scree_df$Filt <- factor(scree_df$Filt, levels = c("basic", "standard", "extra"))
scree_df$BetaType <- factor(scree_df$BetaType, levels = c("dice", "bray", "mori"))

## plot; save as 15_FOX_screeplot; export at 800x800
ggplot(data = scree_df, aes(x=dims, y=stress, color=Filt)) + 
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "firebrick", alpha=0.5) +
  geom_line(alpha=0.5) +
  geom_point() + 
  scale_color_manual(values = c('#9f9244', '#6c42b8', '#628a47')) +
  #scale_x_continuous(breaks = c(seq(1:6))) +
  facet_grid(Method ~ BetaType)  +
  labs(x="dimensions", y="stress", color = "Filtering\nParameter") +
  theme_devon() +
  theme(legend.position="top")

#rm(list=ls(pattern = "scree\\..*"))

## calculate distances and run NMDS for each Filt, Method, Betatest 
nmdsfunction.function <- function(data, test, binaryval, Method, Filt, BetaType) {
  tmp.dist <- vegdist(data, method = test, binary = binaryval)
  tmp.nmds <- metaMDS(tmp.dist, k=3, try=30)
  tmp.nmds.df <- data.frame(tmp.nmds$points)
  tmp.nmds.df$SeqID <- row.names(tmp.nmds.df)
  tmp.nmds.df %>% mutate(Method = Method) %>% mutate(Filt = Filt) %>% mutate(BetaType = BetaType)
}

nmds.dada2.basic.dice <- nmdsfunction.function(mat.dada2.basic, "bray", TRUE, "dada2", "basic", "dice")
nmds.dada2.basic.bray <- nmdsfunction.function(mat.dada2.basic, "bray", FALSE, "dada2", "basic", "bray")
nmds.dada2.basic.mori <- nmdsfunction.function(mat.dada2.basic, "morisita", FALSE, "dada2", "basic", "mori")
nmds.dada2.standard.dice <- nmdsfunction.function(mat.dada2.standard, "bray", TRUE, "dada2", "standard", "dice")
nmds.dada2.standard.bray <- nmdsfunction.function(mat.dada2.standard, "bray", FALSE, "dada2", "standard", "bray")
nmds.dada2.standard.mori <- nmdsfunction.function(mat.dada2.standard, "morisita", FALSE, "dada2", "standard", "mori")
nmds.dada2.extra.dice <- nmdsfunction.function(mat.dada2.extra, "bray", TRUE, "dada2", "extra", "dice")
nmds.dada2.extra.bray <- nmdsfunction.function(mat.dada2.extra, "bray", FALSE, "dada2", "extra", "bray")
nmds.dada2.extra.mori <- nmdsfunction.function(mat.dada2.extra, "morisita", FALSE, "dada2", "extra", "mori")

nmds.deblur.basic.dice <- nmdsfunction.function(mat.deblur.basic, "bray", TRUE, "deblur", "basic", "dice")
nmds.deblur.basic.bray <- nmdsfunction.function(mat.deblur.basic, "bray", FALSE, "deblur", "basic", "bray")
nmds.deblur.basic.mori <- nmdsfunction.function(mat.deblur.basic, "morisita", FALSE, "deblur", "basic", "mori")
nmds.deblur.standard.dice <- nmdsfunction.function(mat.deblur.standard, "bray", TRUE, "deblur", "standard", "dice")
nmds.deblur.standard.bray <- nmdsfunction.function(mat.deblur.standard, "bray", FALSE, "deblur", "standard", "bray")
nmds.deblur.standard.mori <- nmdsfunction.function(mat.deblur.standard, "morisita", FALSE, "deblur", "standard", "mori")
nmds.deblur.extra.dice <- nmdsfunction.function(mat.deblur.extra, "bray", TRUE, "deblur", "extra", "dice")
nmds.deblur.extra.bray <- nmdsfunction.function(mat.deblur.extra, "bray", FALSE, "deblur", "extra", "bray")
nmds.deblur.extra.mori <- nmdsfunction.function(mat.deblur.extra, "morisita", FALSE, "deblur", "extra", "mori")

nmds.vsearch.basic.dice <- nmdsfunction.function(mat.vsearch.basic, "bray", TRUE, "vsearch", "basic", "dice")
nmds.vsearch.basic.bray <- nmdsfunction.function(mat.vsearch.basic, "bray", FALSE, "vsearch", "basic", "bray")
nmds.vsearch.basic.mori <- nmdsfunction.function(mat.vsearch.basic, "morisita", FALSE, "vsearch", "basic", "mori")
nmds.vsearch.standard.dice <- nmdsfunction.function(mat.vsearch.standard, "bray", TRUE, "vsearch", "standard", "dice")
nmds.vsearch.standard.bray <- nmdsfunction.function(mat.vsearch.standard, "bray", FALSE, "vsearch", "standard", "bray")
nmds.vsearch.standard.mori <- nmdsfunction.function(mat.vsearch.standard, "morisita", FALSE, "vsearch", "standard", "mori")
nmds.vsearch.extra.dice <- nmdsfunction.function(mat.vsearch.extra, "bray", TRUE, "vsearch", "extra", "dice")
nmds.vsearch.extra.bray <- nmdsfunction.function(mat.vsearch.extra, "bray", FALSE, "vsearch", "extra", "bray")
nmds.vsearch.extra.mori <- nmdsfunction.function(mat.vsearch.extra, "morisita", FALSE, "vsearch", "extra", "mori")

## combine all scree plot data
nmds_df <- do.call("rbind", mget(ls(pattern = "nmds\\..*")))
rm(list=ls(pattern = "nmds\\..*"))

## join metadata of MonthStart
tmp.fox.meta <- fox %>% distinct(SeqID, MonthStart)
nmds_df <- merge(nmds_df, tmp.fox.meta, by.x="SeqID", by.y="SeqID", all.x=TRUE)
rm(tmp.fox.meta)

## add color for MonthStart palette
monthpal <- c("#2b83ba", "#abdda4", "#fdae61", "#d7191c")

## set levels for plot
nmds_df$MonthStart <- factor(nmds_df$MonthStart, levels = c("April", "May", "September", "October"))
nmds_df$Filt <- factor(nmds_df$Filt, levels = c("basic", "standard", "extra"))

## plots are separated by BetaTest
## save as 16_FOX_nmds_dice; export at 700x700
ggplot(data=nmds_df %>% filter(BetaType=="dice"), 
       aes(x=MDS1, y=MDS2, color = MonthStart)) +
         geom_point() +
         facet_grid(Filt ~ Method) +
  scale_color_manual(values = monthpal) +
  labs(x='mds1', y='mds2', subtitle = "Dice-Sorensen index used to calculate distnaces", color = "") +
  theme_devon() + theme(legend.position = "top")

## save as 16_FOX_nmds_bray; export at 700x700
ggplot(data=nmds_df %>% filter(BetaType=="bray"), 
       aes(x=MDS1, y=MDS2, color = MonthStart)) +
  geom_point() +
  facet_grid(Filt ~ Method) +
  scale_color_manual(values = monthpal) +
  labs(x='mds1', y='mds2', subtitle = "Bray-Curtis index used to calculate distnaces", color = "") +
  theme_devon() + theme(legend.position = "top")

## save as 16_FOX_nmds_mori; export at 700x700
ggplot(data=nmds_df %>% filter(BetaType=="mori"), 
       aes(x=MDS1, y=MDS2, color = MonthStart)) +
  geom_point() +
  facet_grid(Filt ~ Method) +
  scale_color_manual(values = monthpal) +
  labs(x='mds1', y='mds2', subtitle = "Morisita-Horn index used to calculate distnaces", color = "") +
  theme_devon() + theme(legend.position = "top")
