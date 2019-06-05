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


################################################################################
## 2. count instances in which a Classifier has information (non NA) per taxonomic Level
################################################################################

infoSumry <- all_df %>% 
  group_by(level, Classifier) %>% 
  summarise(withData = sum(!is.na(Taxon)),
            noData = sum(is.na(Taxon)),
            frac_withData = withData / (withData + noData),
            frac_noData = noData / (withData + noData))



################################################################################
## 3. plot
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
## set levels
infoSumry$level <- factor(infoSumry$level, levels = c('Class', 'Order', 'Family', 'Genus', 'Species'))
infoSumry$Classifier <- factor(infoSumry$Classifier, levels = c('blast', 'vsearch', 'nbayes', 'sintax', 'bold'))

## plot; save as 'cl_1_guano_nASVsClassifiedPerTaxon'; export at 400x800
ggplot(infoSumry, aes(x=Classifier, y=withData)) +
  geom_bar(stat="identity") +
  facet_grid(level ~ .) +
  #scale_y_continuous(breaks = c(0, .5, 1)) +
  theme_devon() +
  labs(x = "", y="Number of Taxa Classified\n")


################################################################################
## 3. data table of plot info
################################################################################
frac_out <- infoSumry %>% 
  mutate(frac_withData = round(frac_withData, 3)) %>% 
  select(Classifier, level, frac_withData) %>% 
  spread(level, frac_withData)

data_out <- infoSumry %>% 
  select(Classifier, level, withData) %>% 
  spread(level, withData)

write.table(frac_out, 
            file="~/Repos/tidybug/data/text_tables/classify_comps/guanoFracDataClassified.txt", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(data_out, 
            file="~/Repos/tidybug/data/text_tables/classify_comps/guanoDataClassified.txt", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)
