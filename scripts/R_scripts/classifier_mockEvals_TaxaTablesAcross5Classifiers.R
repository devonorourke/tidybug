## script generates tables that compare per Taxa Level the mock community classifications for all 5 Classifiers

library(tidyverse)

################################################################################
## import data
################################################################################


readerfunction <- function(urlpath, classifier) {
  tmp <- read_delim(file = urlpath, delim = "\t", col_names = FALSE) %>%
    rename(., IMalias = X1, Taxon = X2) %>% 
    separate(., 
             col = Taxon, 
             into = c("kingdom", "phylum", "Class", "Order", "Family", "Genus", "Species"), 
             sep = ";") %>% 
    select(-kingdom, -phylum)
  tmp <- as.data.frame(apply(tmp, 2, function(y) (gsub(".__", "", y))))
  tmp <- as.data.frame(apply(tmp, 2, function(y) (gsub("^$|^ $", NA, y))))
  tmp <- as.data.frame(apply(tmp, 2, function(y) (gsub("^NA$", NA, y))))
  tmp[] <- lapply(tmp, as.character)
  tmp %>% gather(key = level, value = Taxon, -IMalias) %>% mutate(Classifier=classifier)
}

exp_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/expected_mockData/mock_expected_taxa.txt"
bl_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/blast/mock_blast_observedTaxa.txt"
vs_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/vsearch/mock_vsearch_observedTaxa.txt"
st_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/sintax/mock_sintax_observedTaxa.txt"
nb_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/nbayes/mock_nbayes_observedTaxa.txt"
bo_url <- "https://github.com/devonorourke/tidybug/raw/master/data/classify_comps/mock_comps/observed_mockData/bold/boldAPI_mock_lca97.txt"

exp_df <- readerfunction(exp_url, "Expected")
bl_df <- readerfunction(bl_url, "blast")
vs_df <- readerfunction(vs_url, "vsearch")
st_df <- readerfunction(st_url, "sintax")
nb_df <- readerfunction(nb_url, "nbayes")
bo_df <- readerfunction(bo_url, "bold")

all_df <- rbind(exp_df, bl_df, vs_df, st_df, nb_df, bo_df)

################################################################################
## generate tables
################################################################################

TaxaTableFunction <- function(Level) {
  tmp <- all_df %>% 
    filter(level == Level) %>%
    select(-level) %>% 
    spread(key = Classifier, value = Taxon)
}

c_tmp <- TaxaTableFunction("Class")
write.table(c_tmp, file = "~/Repos/tidybug/data/text_tables/classify_comps/mock_ClassTaxaTable_all5Classifiers.txt", col.names = TRUE, quote=FALSE, row.names = FALSE)

o_tmp <- TaxaTableFunction("Order")
write.table(o_tmp, file = "~/Repos/tidybug/data/text_tables/classify_comps/mock_OrderTaxaTable_all5Classifiers.txt", col.names = TRUE, quote=FALSE, row.names = FALSE)

f_tmp <- TaxaTableFunction("Family")
write.table(f_tmp, file = "~/Repos/tidybug/data/text_tables/classify_comps/mock_FamilyTaxaTable_all5Classifiers.txt", col.names = TRUE, quote=FALSE, row.names = FALSE)

g_tmp <- TaxaTableFunction("Genus")
write.table(g_tmp, file = "~/Repos/tidybug/data/text_tables/classify_comps/mock_GenusTaxaTable_all5Classifiers.txt", col.names = TRUE, quote=FALSE, row.names = FALSE)

s_tmp <- TaxaTableFunction("Species")
write.table(s_tmp, file = "~/Repos/tidybug/data/text_tables/classify_comps/mock_SpeciesTaxaTable_all5Classifiers.txt", col.names = TRUE, quote=FALSE, row.names = FALSE)
