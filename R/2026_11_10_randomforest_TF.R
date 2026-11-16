###### **********************
## Random Forest Meets TE's
##
## Tomas Fischer
##
## 2026-11-03
##
##***************************

## _ Packages used -------
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("as.factor", "base")
library(viridis)
library(randomForest)
library(Biostrings)
# + scale_color/fill_viridis(discrete = T/F)
theme_set(theme_light())

df_trich <- read.csv("../data/2026-10-26_t_spiralis_1.csv") # , locale=locale(encoding = "Latin1"
View(df_trich)

# df_trich_nt <- df_trich |> 
#   filter(!is.na(nucl_seq))|> 
#   arrange(desc(qry_width)) |> 
#   print()
# 
# rm(df_trich_nt)

# 
# df_trich |>
#   count(repeat_id) |>
#   arrange(desc(n))

#A-rich - 5198
#(AT)n - 1677

ls_largest_id <- c("A-rich", "(AT)n")

df_trich |>
  filter(repeat_id %in% ls_largest_id) |>
  View()

df_largest <- df_trich |> 
  filter(!is.na(nucl_seq)) |> 
  filter(qry_id == "A-rich" | qry_id == "(AT)n") |> 
  View()

