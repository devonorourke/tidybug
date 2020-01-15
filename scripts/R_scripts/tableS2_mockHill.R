## script for Table S2
## describes alpha diversity comparisons using Hill Numbers for MOCK data
## KW and Dunn's tests for comparisons among groups (denoiser software) at each filtering level

library(tidyverse)
library(reshape2)
library(vegan)
library(phyloseq)
library(FSA)

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
df <- df %>% filter(SampleType == "mock")

################################################################################
## data generation
################################################################################

## function to calculate Hill Numbers per Method + Filt (grouping all guano data among all libraries)
hill.function.wrare <- function(data, filter_exp, filter_exp2) {
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
  tmp_raremat = as(otu_table(tmp_rphy), "matrix")
  tmp.hill1 <- renyi(t(tmp_raremat), scales = c(0,1,2), hill=TRUE)
  tmp.hill2 <- data.frame(tmp.hill1, Method, Filt, Lib=row.names(tmp.hill1))
  colnames(tmp.hill2)[1:3] <- c("q=0", "q=1", "q=2")
  tmp_out <- gather(tmp.hill2, key="Hill_qType", value = "Hill_value", c("q=0", "q=1", "q=2"))
}

rarefied.dada2.basic <- hill.function.wrare(df, Method=="dada2", Filt=="basic")
rarefied.dada2.standard <- hill.function.wrare(df, Method=="dada2", Filt=="standard")
rarefied.dada2.extra <- hill.function.wrare(df, Method=="dada2", Filt=="extra")
rarefied.deblur.basic <- hill.function.wrare(df, Method=="deblur", Filt=="basic")
rarefied.deblur.standard <- hill.function.wrare(df, Method=="deblur", Filt=="standard")
rarefied.deblur.extra <- hill.function.wrare(df, Method=="deblur", Filt=="extra")
rarefied.vsearch.basic <- hill.function.wrare(df, Method=="vsearch", Filt=="basic")
rarefied.vsearch.standard <- hill.function.wrare(df, Method=="vsearch", Filt=="standard")
rarefied.vsearch.extra <- hill.function.wrare(df, Method=="vsearch", Filt=="extra")

## merge into single dataframe
all.mock.hill <- rbind(rarefied.dada2.basic, rarefied.dada2.standard, rarefied.dada2.extra, 
                             rarefied.deblur.basic, rarefied.deblur.standard, rarefied.deblur.extra, 
                             rarefied.vsearch.basic, rarefied.vsearch.standard, rarefied.vsearch.extra)
all.mock.hill$Method <- as.factor(all.mock.hill$Method)
rm(list=ls(pattern = "rarefied*"))

write.csv(all.mock.hill, 
          file="~/Documents/nau_projects/guano/mole_ecol_methods_paper/mock_hillvals_table.csv", 
          quote = FALSE, 
          row.names = FALSE)

################################################################################
## KW and Dunn's tests
################################################################################

## Kruskal-Wallis:
kw_hill_mock_function <- function(HillVal){
  basic_dat <- all.mock.hill %>% filter(Hill_qType == HillVal, Filt=="basic")
  stand_dat <- all.mock.hill %>% filter(Hill_qType == HillVal, Filt=="standard")
  extra_dat <- all.mock.hill %>% filter(Hill_qType == HillVal, Filt=="extra")
  list1 = kruskal.test(Hill_value ~ Method, data = basic_dat)
  list2 = kruskal.test(Hill_value ~ Method, data = stand_dat)
  list3 = kruskal.test(Hill_value ~ Method, data = extra_dat)
  list <- list(list1, list2, list3)
  tmp_df = data.frame(Filt = c("basic", "standard", "extra"),
                      HillVal = c(rep(HillVal, 3)),
                      Pvals = c(list[[1]][3]$p.value, list[[2]][3]$p.value, list[[3]][3]$p.value),
                      KWstat = c(list[[1]][1]$statistic, list[[2]][1]$statistic, list[[3]][1]$statistic))
  tmp_df %>% mutate(p.adj = p.adjust(tmp_df$Pvals, method="BH"))
}

q0_kw_df <- kw_hill_mock_function('q=0')
q1_kw_df <- kw_hill_mock_function('q=1')
q2_kw_df <- kw_hill_mock_function('q=2')

kw_df_stats <- rbind(q0_kw_df, q1_kw_df, q2_kw_df)
rm(q0_kw_df, q1_kw_df, q2_kw_df)

write.csv(kw_df_stats, 
          file="~/Documents/nau_projects/guano/mole_ecol_methods_paper/kw_mock_table.csv", 
          quote = FALSE, 
          row.names = FALSE)


## Dunn's test:
dunn_hill_mock_function <- function(HillVal){
  basic_dat <- all.mock.hill %>% filter(Hill_qType == HillVal, Filt=="basic")
  stand_dat <- all.mock.hill %>% filter(Hill_qType == HillVal, Filt=="standard")
  extra_dat <- all.mock.hill %>% filter(Hill_qType == HillVal, Filt=="extra")
  list1 = dunnTest(Hill_value ~ Method, data = basic_dat, method="bh")
  df1 = data.frame(Filt = rep("basic", 3), list1[[2]], Hill_Number = HillVal)
  list2 = dunnTest(Hill_value ~ Method, data = stand_dat, method="bh")
  df2 = data.frame(Filt = rep("standard", 3), list2[[2]], Hill_Number = HillVal)
  list3 = dunnTest(Hill_value ~ Method, data = extra_dat, method="bh")
  df3 = data.frame(Filt = rep("extra", 3), list3[[2]], Hill_Number = HillVal)
  rbind(df1, df2, df3)
}


q0_dunn_df <- dunn_hill_mock_function('q=0')
q1_dunn_df <- dunn_hill_mock_function('q=1')
q2_dunn_df <- dunn_hill_mock_function('q=2')

dunn_df_stats <- rbind(q0_dunn_df, q1_dunn_df, q2_dunn_df)
rm(q0_dunn_df, q1_dunn_df, q2_dunn_df)

write.csv(dunn_df_stats, 
          file="~/Documents/nau_projects/guano/mole_ecol_methods_paper/dunn_mock_table.csv", 
          quote = FALSE, 
          row.names = FALSE)
          

################################################################################
## Plotting diversity estimates of mock libraries
################################################################################
pal3 <- c('#9f9244', '#6c42b8', '#628a47', '#a9d190')

ggplot(all.guano.hill.rare,
       aes(x=Method, y=Hill_value, color=Method)) +
  geom_jitter(width = 0.2, alpha=0.7) +
  facet_grid(Hill_qType ~ Filt) +
  labs(x="", y="sequence variant equivalents\n") +
  scale_color_manual(values=pal3) +
  scale_y_continuous(breaks = c(0, 25, 50),
                     limits = c(0,60)) +
  theme_devon() +
  theme(legend.position="none",
        strip.text.x = element_text(size=13), strip.text.y = element_text(size=13),
        panel.spacing = unit(1.5, "lines"), panel.grid.major.x = element_blank() )

