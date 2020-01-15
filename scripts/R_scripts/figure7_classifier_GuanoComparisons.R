library(tidyverse)
library(reshape2)
library(ggpubr)

################################################################################
## 1. read in data
################################################################################

readerfunction <- function(urlpath, classifier) {
  tmp <- read_delim(file = urlpath, delim = "\t", col_names = FALSE) %>%
    rename(., HashID = X1, Taxon = X2) %>% 
    separate(., 
             col = Taxon, 
             into = c("kingdom", "phylum", "Class", "Order", "Family", "Genus", "Species"), 
             sep = ";") %>% 
    select(-kingdom, -phylum)
  tmp <- as.data.frame(apply(tmp, 2, function(y) (gsub(".__", "", y))))
  tmp <- as.data.frame(apply(tmp, 2, function(y) (gsub("^$|^ $", NA, y))))
  tmp <- as.data.frame(apply(tmp, 2, function(y) (gsub("Ambiguous taxa", NA, y))))
  tmp <- as.data.frame(apply(tmp, 2, function(y) (gsub("^NA$", NA, y))))
  tmp[] <- lapply(tmp, as.character)
  tmp %>% 
    gather(key = level, value = Taxon, -HashID) %>% 
    mutate(Classifier=classifier)
}

bl_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/guano_comps/blast/guano_blast_taxonomy.txt"
vs_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/guano_comps/vsearch/guano_vsearch_taxonomy.txt"
st_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/guano_comps/sintax/guano_sintax_taxonomy.txt"
nb_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/guano_comps/nbayes/guano_nbayes_taxonomy.txt"
bo_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/guano_comps/bold/guano_bold_taxonomy.txt"

bl_df <- readerfunction(bl_url, "blast")
vs_df <- readerfunction(vs_url, "vsearch")
st_df <- readerfunction(st_url, "sintax")
nb_df <- readerfunction(nb_url, "nbayes")
bo_df <- readerfunction(bo_url, "bold")

## combine df's
all_df <- rbind(bl_df, vs_df, st_df, nb_df, bo_df)
rm(bl_df, vs_df, st_df, nb_df, bo_df)
rm(bl_url, nb_url, st_url, vs_url, bo_url, readerfunction)

## generate dataset for total taxa with information for first of two plots
plot1_dat <- all_df %>% 
  filter(!is.na(Taxon)) %>% 
  group_by(Classifier, level) %>% 
  summarise(nTaxa=n()) %>% 
  filter(level != 'Class')

################################################################################
## 2. create per-taxa-Level data frames for Upset plots
################################################################################

plotdatfunction <- function(TaxaLevel) {
  tmp_df <- all_df %>% 
    filter(level == TaxaLevel) %>% 
    spread(key = Classifier, value = Taxon)
  match1_tmp <- tmp_df %>% 
    mutate(allMatch = case_when(blast == vsearch & blast == sintax & blast == nbayes & blast == bold ~ 'match')) %>% 
    select(HashID, allMatch)
  drop_match1 <- match1_tmp %>% filter(allMatch == "match") %>% pull(HashID)
  match2_tmp <- tmp_df %>% 
    filter(!HashID %in% drop_match1) %>% 
    mutate(Not_bo = case_when(blast == vsearch & blast == sintax & blast==nbayes & is.na(bold) ~ 'match',
                              blast == vsearch & blast == sintax & blast==nbayes & !is.na(bold) ~ 'diff')) %>% 
    mutate(Not_nb = case_when(blast == vsearch & blast == sintax & blast==bold & is.na(nbayes) ~ 'match',
                              blast == vsearch & blast == sintax & blast==bold & !is.na(nbayes) ~ 'diff')) %>% 
    mutate(Not_st = case_when(blast == vsearch & blast == nbayes & blast == bold & is.na(sintax) ~ 'match',
                              blast == vsearch & blast == nbayes & blast == bold & !is.na(sintax) ~ 'diff')) %>%   
    mutate(Not_vs = case_when(blast == sintax & blast == nbayes & blast == bold & is.na(vsearch) ~ 'match',
                              blast == sintax & blast == nbayes & blast == bold & !is.na(vsearch) ~ 'diff')) %>% 
    mutate(Not_bl = case_when(vsearch == sintax & vsearch == nbayes & vsearch == bold & is.na(blast) ~ 'match',
                              vsearch == sintax & vsearch == nbayes & vsearch == bold & !is.na(blast) ~ 'diff')) %>%
    select(-blast, -vsearch, -sintax, -nbayes, -bold, -level)
  
  drop_match2 <- match2_tmp %>% filter(!is.na(Not_bo) | !is.na(Not_nb) | !is.na(Not_st) | !is.na(Not_vs) | !is.na(Not_bl)) %>% pull(HashID)
  drop_match3 <- c(drop_match1, drop_match2)
  
  match3_tmp <- tmp_df %>% 
    filter(!HashID %in% drop_match3) %>% 
    mutate(bl_vs_st = case_when(blast == vsearch & blast == sintax & is.na(nbayes) &  is.na(bold) ~ 'match',
                                blast == vsearch & blast == sintax & !is.na(nbayes) & is.na(bold) ~ 'diff',
                                blast == vsearch & blast == sintax & is.na(nbayes) & !is.na(bold) ~ 'diff',
                                blast == vsearch & blast == sintax & !is.na(nbayes) & !is.na(bold) ~ 'diff')) %>% 
    mutate(bl_vs_nb = case_when(blast == vsearch & blast == nbayes & is.na(sintax) &  is.na(bold) ~ 'match',
                                blast == vsearch & blast == nbayes & !is.na(sintax) & is.na(bold) ~ 'diff',
                                blast == vsearch & blast == nbayes & is.na(sintax) & !is.na(bold) ~ 'diff',
                                blast == vsearch & blast == nbayes & !is.na(sintax) & !is.na(bold) ~ 'diff')) %>% 
    mutate(bl_vs_bo = case_when(blast == vsearch & blast == bold & is.na(nbayes) &  is.na(sintax) ~ 'match',
                                blast == vsearch & blast == bold & !is.na(nbayes) & is.na(sintax) ~ 'diff',
                                blast == vsearch & blast == bold & is.na(nbayes) & !is.na(sintax) ~ 'diff',
                                blast == vsearch & blast == bold & !is.na(nbayes) & !is.na(sintax) ~ 'diff')) %>% 
    mutate(bl_st_bo = case_when(blast == bold & blast == sintax & is.na(nbayes) &  is.na(vsearch) ~ 'match',
                                blast == bold & blast == sintax & !is.na(nbayes) & is.na(vsearch) ~ 'diff',
                                blast == bold & blast == sintax & is.na(nbayes) & !is.na(vsearch) ~ 'diff',
                                blast == bold & blast == sintax & !is.na(nbayes) & !is.na(vsearch) ~ 'diff')) %>% 
    mutate(bl_nb_st = case_when(blast == nbayes & blast == sintax & is.na(vsearch) &  is.na(bold) ~ 'match',
                                blast == nbayes & blast == sintax & !is.na(vsearch) & is.na(bold) ~ 'diff',
                                blast == nbayes & blast == sintax & is.na(vsearch) & !is.na(bold) ~ 'diff',
                                blast == nbayes & blast == sintax & !is.na(vsearch) & !is.na(bold) ~ 'diff')) %>% 
    mutate(bl_nb_bo = case_when(blast == nbayes & blast == bold & is.na(sintax) &  is.na(vsearch) ~ 'match',
                                blast == nbayes & blast == bold & !is.na(sintax) & is.na(vsearch) ~ 'diff',
                                blast == nbayes & blast == bold & is.na(sintax) & !is.na(vsearch) ~ 'diff',
                                blast == nbayes & blast == bold & !is.na(sintax) & !is.na(vsearch) ~ 'diff')) %>% 
    mutate(vs_st_nb = case_when(vsearch == nbayes & vsearch == sintax & is.na(blast) &  is.na(bold) ~ 'match',
                                vsearch == nbayes & vsearch == sintax & !is.na(blast) & is.na(bold) ~ 'diff',
                                vsearch == nbayes & vsearch == sintax & is.na(blast) & !is.na(bold) ~ 'diff',
                                vsearch == nbayes & vsearch == sintax & !is.na(blast) & !is.na(bold) ~ 'diff')) %>% 
    mutate(vs_st_bo = case_when(vsearch == bold & vsearch == sintax & is.na(blast) &  is.na(nbayes) ~ 'match',
                                vsearch == bold & vsearch == sintax & !is.na(blast) & is.na(nbayes) ~ 'diff',
                                vsearch == bold & vsearch == sintax & is.na(blast) & !is.na(nbayes) ~ 'diff',
                                vsearch == bold & vsearch == sintax & !is.na(blast) & !is.na(nbayes) ~ 'diff')) %>% 
    mutate(vs_nb_bo = case_when(vsearch == bold & vsearch == nbayes & is.na(blast) &  is.na(sintax) ~ 'match',
                                vsearch == bold & vsearch == nbayes & !is.na(blast) & is.na(sintax) ~ 'diff',
                                vsearch == bold & vsearch == nbayes & is.na(blast) & !is.na(sintax) ~ 'diff',
                                vsearch == bold & vsearch == nbayes & !is.na(blast) & !is.na(sintax) ~ 'diff')) %>% 
    mutate(st_nb_bo = case_when(sintax == nbayes & bold == sintax & is.na(vsearch) &  is.na(blast) ~ 'match',
                                sintax == nbayes & bold == sintax & !is.na(vsearch) & is.na(blast) ~ 'diff',
                                sintax == nbayes & bold == sintax & is.na(vsearch) & !is.na(blast) ~ 'diff',
                                sintax == nbayes & bold == sintax & !is.na(vsearch) & !is.na(blast) ~ 'diff')) %>% 
    select(-blast, -vsearch, -sintax, -nbayes, -bold, -level)
  
  drop_match3 <- match3_tmp %>% 
    filter(!is.na(bl_vs_st) | !is.na(bl_vs_nb) | !is.na(bl_vs_bo) | !is.na(bl_st_bo) | !is.na(bl_nb_st) | 
             !is.na(bl_nb_bo) | !is.na(vs_st_nb) | !is.na(vs_st_bo) | !is.na(vs_nb_bo) | !is.na(st_nb_bo)) %>% 
    pull(HashID)
  drop_match4 <- c(drop_match1, drop_match2, drop_match3)
  
  match4_tmp <- tmp_df %>% 
    filter(!HashID %in% drop_match4) %>% 
    mutate(bl_vs = case_when(blast == vsearch & is.na(sintax) & is.na(nbayes) &  is.na(bold) ~ 'match',
                             blast == vsearch & !is.na(nbayes) & is.na(bold) & is.na(sintax) ~ 'diff',
                             blast == vsearch & is.na(nbayes) & !is.na(bold) & is.na(sintax) ~ 'diff',
                             blast == vsearch & is.na(nbayes) & is.na(bold) & !is.na(sintax) ~ 'diff')) %>% 
    mutate(bl_nb = case_when(blast == nbayes & is.na(sintax) &  is.na(bold) & is.na(vsearch) ~ 'match',
                             blast == nbayes & !is.na(sintax) & is.na(bold) & is.na(vsearch) ~ 'diff',
                             blast == nbayes & is.na(sintax) & !is.na(bold) & is.na(vsearch) ~ 'diff',
                             blast == nbayes & is.na(sintax) & is.na(bold) & !is.na(vsearch) ~ 'diff')) %>% 
    mutate(bl_bo = case_when(blast == bold & is.na(nbayes) &  is.na(sintax) & is.na(vsearch) ~ 'match',
                             blast == bold & !is.na(nbayes) & is.na(sintax) & is.na(vsearch) ~ 'diff',
                             blast == bold & is.na(nbayes) & !is.na(sintax) & is.na(vsearch) ~ 'diff',
                             blast == bold & is.na(nbayes) & is.na(sintax) & !is.na(vsearch) ~ 'diff')) %>% 
    mutate(bl_st = case_when(blast == sintax & is.na(nbayes) &  is.na(vsearch) & is.na(bold) ~ 'match',
                             blast == sintax & !is.na(nbayes) & is.na(vsearch) & is.na(bold) ~ 'diff',
                             blast == sintax & is.na(nbayes) & !is.na(vsearch) & is.na(bold) ~ 'diff',
                             blast == sintax & is.na(nbayes) & is.na(vsearch) & !is.na(bold) ~ 'diff')) %>% 
    mutate(vs_st = case_when(vsearch == sintax & is.na(blast) &  is.na(bold) & is.na(nbayes) ~ 'match',
                             vsearch == sintax & !is.na(blast) & is.na(bold) & is.na(nbayes) ~ 'diff',
                             vsearch == sintax & is.na(blast) & !is.na(bold) & is.na(nbayes) ~ 'diff',
                             vsearch == sintax & is.na(blast) & is.na(bold) & !is.na(nbayes) ~ 'diff')) %>% 
    mutate(vs_nb = case_when(vsearch == nbayes & is.na(blast) &  is.na(bold) & is.na(sintax) ~ 'match',
                             vsearch == nbayes & !is.na(blast) & is.na(bold) & is.na(sintax) ~ 'diff',
                             vsearch == nbayes & is.na(blast) & !is.na(bold) & is.na(sintax) ~ 'diff',
                             vsearch == nbayes & is.na(blast) & is.na(bold) & !is.na(sintax) ~ 'diff')) %>%
    mutate(vs_bo = case_when(vsearch == bold & is.na(blast) &  is.na(sintax) & is.na(nbayes) ~ 'match',
                             vsearch == bold & !is.na(blast) & is.na(sintax) & is.na(nbayes) ~ 'diff',
                             vsearch == bold & is.na(blast) & !is.na(sintax) & is.na(nbayes) ~ 'diff',
                             vsearch == bold & is.na(blast) & is.na(sintax) & !is.na(nbayes) ~ 'diff')) %>%
    mutate(st_nb = case_when(sintax == nbayes & is.na(blast) &  is.na(bold) & is.na(vsearch) ~ 'match',
                             sintax == nbayes & !is.na(blast) & is.na(bold) & is.na(vsearch) ~ 'diff',
                             sintax == nbayes & is.na(blast) & !is.na(bold) & is.na(vsearch) ~ 'diff',
                             sintax == nbayes & is.na(blast) & is.na(bold) & !is.na(vsearch) ~ 'diff')) %>% 
    mutate(st_bo = case_when(sintax == bold & is.na(blast) &  is.na(nbayes) & is.na(vsearch) ~ 'match',
                             sintax == bold & !is.na(blast) & is.na(nbayes) & is.na(vsearch) ~ 'diff',
                             sintax == bold & is.na(blast) & !is.na(nbayes) & is.na(vsearch) ~ 'diff',
                             sintax == bold & is.na(blast) & is.na(nbayes) & !is.na(vsearch) ~ 'diff')) %>% 
    mutate(nb_bo = case_when(nbayes == bold & is.na(blast) &  is.na(sintax) & is.na(vsearch) ~ 'match',
                             nbayes == bold & !is.na(blast) & is.na(sintax) & is.na(vsearch) ~ 'diff',
                             nbayes == bold & is.na(blast) & !is.na(sintax) & is.na(vsearch) ~ 'diff',
                             nbayes == bold & is.na(blast) & is.na(sintax) & !is.na(vsearch) ~ 'diff')) %>% 
    select(-blast, -vsearch, -sintax, -nbayes, -bold, -level)
  
  drop_match4 <- match4_tmp %>% 
    filter(!is.na(bl_vs) | !is.na(bl_nb) | !is.na(bl_bo) | !is.na(bl_st) | !is.na(vs_st) | 
             !is.na(vs_nb) | !is.na(vs_bo) | !is.na(st_nb) | !is.na(st_bo) | !is.na(nb_bo)) %>% 
    pull(HashID)
  drop_match5 <- c(drop_match1, drop_match2, drop_match3, drop_match4)
  
  match5_tmp <- tmp_df %>% 
    filter(!HashID %in% drop_match5) %>% 
    mutate(uniq_bl = case_when(!is.na(blast) & is.na(vsearch) & is.na(sintax) & is.na(nbayes) & is.na(bold) ~ 'match')) %>% 
    mutate(uniq_vs = case_when(!is.na(vsearch) & is.na(blast) & is.na(sintax) & is.na(nbayes) & is.na(bold) ~ 'match')) %>% 
    mutate(uniq_st = case_when(!is.na(sintax) & is.na(blast) & is.na(vsearch) & is.na(nbayes) & is.na(bold) ~ 'match')) %>% 
    mutate(uniq_nb = case_when(!is.na(nbayes) & is.na(blast) & is.na(vsearch) & is.na(sintax) & is.na(bold) ~ 'match')) %>% 
    mutate(uniq_bo = case_when(!is.na(bold) & is.na(blast) & is.na(vsearch) & is.na(sintax) & is.na(nbayes) ~ 'match')) %>% 
    select(-level, -blast, -vsearch, -sintax, -nbayes, -bold)
  
  match_joind <- left_join(match1_tmp, match2_tmp) %>% 
    left_join(., match3_tmp) %>% 
    left_join(., match4_tmp) %>% 
    left_join(., match5_tmp)
  
  tmp_out <- as.data.frame(match_joind %>% 
                             gather(key = "MatchType", value = "Value", -HashID, na.rm = TRUE) %>% 
                             group_by(MatchType, Value) %>% 
                             summarise(counts=n()))
#  indexthis <- tmp_out$Value == "diff"
#  tmp_out$counts[indexthis] <- -abs(tmp_out$counts[indexthis])
  tmp_out %>% mutate(Level = TaxaLevel)
}

species_dat <- plotdatfunction("Species")
genus_dat <- plotdatfunction("Genus")
family_dat <- plotdatfunction("Family")
order_dat <- plotdatfunction("Order")
class_dat <- plotdatfunction("Class")

all_dat <- rbind(species_dat, genus_dat, family_dat, order_dat, class_dat)
rm(species_dat, genus_dat, family_dat, order_dat, class_dat)

all_dat <- all_dat %>% 
  group_by(Level) %>% 
  mutate(TotalCounts = sum(counts)) %>% 
  mutate(FracCounts = counts/TotalCounts)

################################################################################
## 3. plots
################################################################################
## theme for plot
theme_devon <- function () { 
  theme_bw(base_size=12, base_family="Avenir") %+replace% 
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill="transparent", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA)
    )
}


## going to limit our plots to just Class, Family, and Species rank (rather than plotting all 5...)
plotTaxaLevels <- c("Order", "Family", "Genus", "Species") ## make sure to change Levels being set as factors below if this is changed!
plot_dat <- all_dat %>% filter(Level %in% plotTaxaLevels)

## set levels
plot_dat$Level <- factor(plot_dat$Level, levels = c("Order", "Family", "Genus", "Species"))  ## change this if you update plotTaxaLevels term!

## main plot
ptop <- ggplot(plot_dat %>% filter(FracCounts >= 0.03), 
               aes(x = reorder(MatchType, -counts, sum), 
                   y = counts)) +
  geom_bar(stat="identity", position=position_dodge(width = 1)) +
  facet_wrap(~Level, nrow=4, scales = "free_y") +
  scale_y_continuous(breaks = scales::pretty_breaks(3), limits = c(-800, NA)) +
  theme_devon() +
  labs(x="", y="shared taxa names\n") +
  geom_text(aes(x=reorder(MatchType, -counts, sum),
                y=-200,
                label=counts),
            data=plot_dat %>% filter(FracCounts >= 0.03),
            size = 3, 
            vjust = 1) +
  theme(strip.text = element_text(size = 13),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(1,1,0,1), "lines"))

## dot plot below
classifier_plotdat <- c(rep("blast", 9), rep("vsearch", 9), rep("nbayes", 9), rep("sintax", 9), rep("bold", 9))
matchtype_plotdat <- rep(c("allMatch", "Not_bo", "st_nb", "uniq_st", "uniq_nb", 
                           "uniq_bo", "bl_vs_st", "st_nb_bo", "Not_nb"),5)
testplotdat <- data.frame(classifier_plotdat, matchtype_plotdat)
testplotdat$dot <- c(1,1,0,0,0,0,1,0,1,      # blast
                     1,1,0,0,0,0,1,0,1,      # vsearch
                     1,1,1,0,0,1,0,1,0,      # nbayes
                     1,1,1,1,0,0,1,1,1,      # sintax
                     1,0,0,0,1,0,0,1,1)      # bold

testplotdat$dot2 <- ifelse(testplotdat$dot==0, NA, 1)
testplotdat$dot2 <- as.character(testplotdat$dot2)

testplotdat$classifier_plotdat <- factor(testplotdat$classifier_plotdat, 
                                         levels = c("bold", "sintax", "nbayes", "vsearch", "blast"))

testplotdat$matchtype_plotdat <- factor(testplotdat$matchtype_plotdat, 
                                        levels = c("allMatch", "Not_bo", "st_nb", "uniq_st", "uniq_nb", 
                                                   "uniq_bo", "bl_vs_st", "st_nb_bo", "Not_nb"))

pbot <- ggplot(testplotdat, 
               aes(x=matchtype_plotdat, y=classifier_plotdat, color=dot2, group=matchtype_plotdat)) +
  geom_point(size = 3) +
  geom_line(data = testplotdat, 
            aes(x=matchtype_plotdat, 
                y=classifier_plotdat, 
                group=matchtype_plotdat),
            color="black") +
  scale_color_manual(values = "gray20", na.value="gray75") +
  theme_devon() +
  labs(x="", y="", color="") +
  theme(legend.position = "none", 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(0,1,1,1), "lines"))


## stitch all these various plots together
pright <- ggarrange(ptop, pbot, nrow = 2, heights = c(1, .15), align = "v")

plot1_dat$level <- factor(plot1_dat$level, levels = c("Order", "Family", "Genus", "Species"))
plot1_dat$Classifier <- factor(plot1_dat$Classifier, levels = c("blast", "vsearch", "nbayes", "sintax", "bold"))

ptopleft <- ggplot(plot1_dat, aes(x=Classifier, y=nTaxa)) +
  geom_bar(stat="identity") +
  scale_y_continuous(limits = c(-1200,NA)) +
  facet_wrap(~ level, nrow=4) +
  labs(x="", y="taxa with names classified\n") +
  theme_devon() +
  theme(strip.text = element_text(size = 13),
        plot.margin = unit(c(1,1,0,1), "lines")) +
  geom_text(aes(x=Classifier,
                y=-400,
                label=nTaxa),
            data=plot1_dat,
            size = 3, 
            vjust = 1)
  

pleft <- ggarrange(ptopleft, NULL, nrow = 2, heights = c(1, .15), align = "v")

## finally, pull it all together
## save as 'figure7_classifierComps'; export at 900w x 750h
## the x-axis labels of Figure 7A and the grid/dot labels in Figure 7B were manually edited in Adobe Illustrator
ggarrange(pleft, pright, ncol = 2, widths = c(.75, 1), align = "h", labels = c("A", "B"))
