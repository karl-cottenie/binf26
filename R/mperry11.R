###### **********************
## Random Forest Activity
##
## Matthew Perry
##
## 2026-11-10
##
##***************************

## _ Packages used -------
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(viridis)
library(randomForest)
library(Biostrings)
# + scale_color/fill_viridis(discrete = T/F)
theme_set(theme_light())

dfBOLD <- read_tsv("../data/wasps.tsv")
View(dfBOLD)