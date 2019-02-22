library(tidyverse)
library(FSA)
library(vegan)
library(reshape2)
library(viridis)
library(ggcorrplot)
library(ggpubr)

## import data and select mock samples:
df <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/all.filtmethods.df.csv.gz")
## filter out mock samples:
df <- df %>% filter(SampleType != "mock") %>% select(-StudyID, -Alias)
df$Labeler <- paste(df$Method, df$Filt, df$Library, sep="-")

## create similar dataframe with rarefied data:
rrarewithdrop <- 
  function(x, sample) 
  {
    rrarefy(x[rowSums(x) >= sample, , drop = FALSE], sample)
  }

rare.function <- function(data, filter_exp, filter_exp2) {
  filter_exp_enq <- enquo(filter_exp)
  filter_exp_enq2 <- enquo(filter_exp2)
  tmp.df <- data %>% filter(!!filter_exp_enq) %>% filter(!!filter_exp_enq2)
  meta.df <- tmp.df %>% distinct(SeqID, Library, Method, Filt, SampleType, Labeler)
  tmp.mat <- dcast(tmp.df, SeqID ~ HashID, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SeqID
  SeqID <- row.names(tmp.mat)
  tmp.mat$SeqID <- NULL
  tmp.mat2 <- data.frame(rrarewithdrop(tmp.mat, 5000))
  tmp.mat2$SeqID <- row.names(tmp.mat2)
  tmp.df2 <- melt(data = tmp.mat2, id.vars = "SeqID",
                  variable.name = "HashID", 
                  value.name = "Reads") %>% 
    filter(Reads > 0) %>%
    mutate(HashID=sub('.', '', HashID))
  tmp.df2 <- merge(tmp.df2, meta.df)
}

## generate rarefied dataset:
tmp.dada2.basic <- rare.function(df, Method=="dada2", Filt=="basic")
tmp.dada2.standard <- rare.function(df, Method=="dada2", Filt=="standard")
tmp.dada2.extra <- rare.function(df, Method=="dada2", Filt=="extra")
tmp.deblur.basic <- rare.function(df, Method=="deblur", Filt=="basic")
tmp.deblur.standard <- rare.function(df, Method=="deblur", Filt=="standard")
tmp.deblur.extra <- rare.function(df, Method=="deblur", Filt=="extra")
tmp.vsearch.basic <- rare.function(df, Method=="vsearch", Filt=="basic")
tmp.vsearch.standard <- rare.function(df, Method=="vsearch", Filt=="standard")
tmp.vsearch.extra <- rare.function(df, Method=="vsearch", Filt=="extra")
rare.df <- rbind(tmp.dada2.basic, tmp.dada2.standard, tmp.dada2.extra,
                 tmp.deblur.basic, tmp.deblur.standard, tmp.deblur.extra,
                 tmp.vsearch.basic, tmp.vsearch.standard, tmp.vsearch.extra)
rm(list=ls(pattern="tmp.*"))


## combine both datasets and add Rarefy label
df$RarefyType <- "unrarefied"
rare.df$RarefyType <- "rarefied"
df <- rbind(df, rare.df)
rm(rare.df)

## generate summaries for number of ASVs per sample, per pipeline and filtering method and rarefying type
unrare.df.ASVcounts <- df %>% filter(RarefyType=="unrarefied") %>%group_by(Method, Filt, SeqID, Library, RarefyType) %>% summarise(ASVcounts=n())
rare.df.ASVcounts <- df %>% filter(RarefyType=="rarefied") %>%group_by(Method, Filt, SeqID, Library, RarefyType) %>% summarise(ASVcounts=n())
rm(df)

## 3way anova, one per RarefyType
unrare.anova <- aov(ASVcounts ~ Library * Method * Filt, data=unrare.df.ASVcounts)
summary(unrare.anova)
rare.anova <- aov(ASVcounts ~ Library * Method * Filt, data=rare.df.ASVcounts)
summary(rare.anova)
  ## note the 3way interaction is no longer significant

## normality check test
shapiro.test(unrare.df.ASVcounts %>% filter(Method=="dada2" & Filt=="basic" & RarefyType=="unrarefied") %>% pull(ASVcounts))   ## not normal
shapiro.test(rare.df.ASVcounts %>% filter(Method=="dada2" & Filt=="basic" & RarefyType=="rarefied") %>% pull(ASVcounts))   ## not normal
shapiro.test(unrare.df.ASVcounts %>% filter(Method=="deblur" & Filt=="extra" & RarefyType=="unrarefied") %>% pull(ASVcounts))   ## not normal
  ## caution when interpreting ANOVA - distributions not normally distributed!

## kruskal wallis test instead?
## add "Labeler" to group data by Library+Method+Filt
## reformat $Labeler from character to factor
unrare.df.ASVcounts$Labeler <- paste(unrare.df.ASVcounts$Method, unrare.df.ASVcounts$Filt, unrare.df.ASVcounts$Library, sep="-")
unrare.df.ASVcounts$Labeler <- as.factor(unrare.df.ASVcounts$Labeler)
rare.df.ASVcounts$Labeler <- paste(rare.df.ASVcounts$Method, rare.df.ASVcounts$Filt, rare.df.ASVcounts$Library, sep="-")
rare.df.ASVcounts$Labeler <- as.factor(rare.df.ASVcounts$Labeler)

kw.unrare <- kruskal.test(ASVcounts ~ Labeler, data=unrare.df.ASVcounts)
kw.rare <- kruskal.test(ASVcounts ~ Labeler, data=rare.df.ASVcounts)
kw.unrare   ## significant difference in 'Labeler' factor, but can't tell what...
kw.rare     ## significant difference here too

## dunn test and resulting pairwise correlation heatmap:
## big messy function to generate the data.frame for plotting a correlation matrix of the pairwise values from Dunn test
corrplot.function <- function(data, filter_exp) {
  tmp.dunn = dunnTest(ASVcounts ~ Labeler, data=data, method="bh")
  tmp.dunn.df <- tmp.dunn$res
  tmp.dunn.df <- separate(data=tmp.dunn.df, col = "Comparison", into=c("Group1", "Group2"), sep=" - ")
  tmp.dunn.df$P.adj <- round(tmp.dunn.df$P.adj,3)
  tmp.dunn.mat <- dcast(tmp.dunn.df, Group1 ~ Group2, value.var = "P.adj", drop = FALSE)
  row.names(tmp.dunn.mat) <- tmp.dunn.mat$Group1
  tmp.dunn.mat$Group1 <- NA
  tmp.dunn.mat <- as.matrix(tmp.dunn.mat)
  botrow <- rep(NA, length(colnames(tmp.dunn.mat)))
  tmp.dunn.mat <- rbind(tmp.dunn.mat, botrow)
  lengthcols <- length(colnames(tmp.dunn.mat))
  selectrow <- colnames(tmp.dunn.mat)[lengthcols]
  rownames(tmp.dunn.mat)[lengthcols] <- selectrow
  rownames(tmp.dunn.mat)[1] -> selectcol
  colnames(tmp.dunn.mat)[1] <- selectcol
  colnames(tmp.dunn.mat) <- rownames(tmp.dunn.mat)
  tmp.dunn.mat[lower.tri(tmp.dunn.mat)] = t(tmp.dunn.mat)[lower.tri(tmp.dunn.mat)]
  tmp.dunn.mat[is.na(tmp.dunn.mat)] <- 1
  tmp2 <- tmp.dunn.mat + 10
  tmp3 <- tril(tmp2)
  tmp3 <- tmp3 - 10
  tmp4 <- as.matrix(tmp3)
  melt(tmp4, varnames = c("Group1", "Group2"), value.name = "P.adj") %>% filter(P.adj >= 0)
}

## generate data
## rarefied data:
unrare.corr.data <- corrplot.function(unrare.df.ASVcounts, RarefyType=="unrarefied")

## unrarefied data:
rare.corr.data <- corrplot.function(rare.df.ASVcounts, RarefyType=="rarefied")

## heatmapt plots
## generate separate plots and combine with cowplot
## unrarefied
viridis(3,option="plasma")
p1 <- ggplot(unrare.corr.data, aes(x=Group1, y=Group2, fill=P.adj)) +
  geom_tile(color="black") +
  scale_y_discrete(position = "right") +
  theme(legend.position = "left") +
  labs(x="", y="", fill="BH-adjusted\np.value") +
  scale_fill_viridis_c(direction = -1, option = "plasma") +
  theme_devon() +
  theme(axis.text.x = element_text(angle=22.5, hjust=1, size = 7),
        axis.text.y = element_text(size=7), legend.position = "right",
        panel.background = element_blank(), panel.border = element_blank(), panel.grid = element_blank(),
        plot.margin = unit(c(.5,.5,.5,2.5), "cm"))

## rarefied  
p2 <- ggplot(rare.corr.data, aes(x=Group1, y=Group2, fill=P.adj)) +
  geom_tile(color="black") +
  scale_y_discrete(position = "right") +
  theme(legend.position = "left") +
  labs(x="", y="", fill="BH-adjusted\np.value") +
  scale_fill_viridis_c(direction = -1, option = "plasma") +
  theme_devon() +
  theme(axis.text.x = element_text(angle=22.5, hjust=1, size = 7),
        axis.text.y = element_text(size=7), legend.position = "right",
        panel.background = element_blank(), panel.border = element_blank(), panel.grid = element_blank(),
        plot.margin = unit(c(.5,.5,.5,2.5), "cm"))

## group together
#legend = get_legend(p1)
ggarrange(p1, p2, common.legend = TRUE)


#### 
#### Wilcox pairwise test not used in publication; code preserved
####

## not run:
#wilcox = pairwise.wilcox.test(df.asvspersample$HashCounts, 
#                              df.asvspersample$Labeler, 
#                          p.adjust.method="BH")
#wilcox.dist = as.dist(wilcox$p.value)
#wilcox.df <- dist2list(wilcox.dist)
#wilcox.df$value <- round(wilcox.df$value,3)

## function for `dist2list` from 'devtools::install_github("vmikk/metagMisc")'
#dist2list <- function (dist, tri=TRUE) {
#  if (!class(dist) == "dist") { stop("Error: The input data must be a dist object.\n") }
#  
#  dat <- as.data.frame(as.matrix(dist))
#  if (is.null(names(dat))) {
#    rownames(dat) <- paste(1:nrow(dat))
#  }
#  value <- stack(dat)$values
#  rnames <- rownames(dat)
#  namecol <- expand.grid(rnames, rnames)
#  colnames(namecol) <- c("col", "row")
#  res <- data.frame(namecol, value)
  
#  if(tri == TRUE){    # return only lower triangular part of dist
#    res <- res[-which(upper.tri(as.matrix(dist), diag = T)), ]
#  }
#  
#  return(res)
#}