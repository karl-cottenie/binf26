##***************************
## title
##
## Kristian Gonzales
##
## 2026-11-10
##
##***************************

## _ Packages used -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(viridis)
library(Biostrings)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())

# Startup ends here

## _ Comment codes ------
# Coding explanations (#, often after the code, but not exclusively)
# Code organization (## XXXXX -----)
# Justification for a section of code ## XXX
# Dead end analyses because it did not work, or not pursuing this line of inquiry (but leave it in as a trace of it, to potentially solve this issue, or avoid making the same mistake in the future # (>_<) 
# Solutions/results/interpretations (#==> XXX)
# Reference to manuscript pieces, figures, results, tables, ... # (*_*)
# TODO items #TODO
# names for data frames (df_name), lists (ls_name), ... (Thanks Jacqueline May)snippet lib
library(title)


## Differentiate the two most common repeat_id groups? (Random forest)
#2026_11_10 Kristian Gonzales Active participation




df_te_tspiralis |> 
  group_by(repeat_id) |> 
  count() |> 
  arrange(desc(n))

ls_mstcmn <- c("A-rich","(AT)n")

df_mstcmn_rpid <- df_te_tspiralis |> 
  filter(repeat_id %in% ls_mstcmn)

df_rpid_pdiv <-  df_te_tspiralis |> 
  filter(repeat_id %in% ls_mstcmn)|> 
  select(repeat_id, perc_div)

ggplot(df_rpid_pdiv) +
  geom_boxplot(aes(x = repeat_id, y = perc_div))