## create geometric mean function
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


## import data and select mock samples:
df <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/all.filtmethods.df.csv.gz")
mock <- df %>% filter(SampleType == "mock")
## import mock matches (exact, partial, and misses)
mockExact <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/mock_community/mock.exactHits.txt", col_names = FALSE)
mockExact$Type <- "exact"
colnames(mockExact)[1] <- "HashID"
mockPartial <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/mock_community/mock.partialHits.txt", col_names = FALSE)
mockPartial$Type <- "partial"
colnames(mockPartial)[1] <- "HashID"
mockMiss <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/mock_community/mock.partialMisses.txt", col_names = FALSE)
mockMiss$Type <- "miss"
colnames(mockMiss)[1] <- "HashID"
mockHitTypes <- rbind(mockExact, mockPartial, mockMiss)
rm(mockExact, mockPartial, mockMiss)

## merge mockHitTypes $Type to mock object (by $HashID)
mock <- merge(mock, mockHitTypes)
rm(mockHitTypes)

## calculate read sums for "Miss" data
## one way: taking the max value observed
mockPartialSumry <- mock %>% 
  filter(Type == "miss") %>%
  group_by(Method, Library) %>% summarise(MaxFilt=max(Reads), MeanFilt=round(mean(Reads)), GeoFilt=round(gm_mean(Reads)))

