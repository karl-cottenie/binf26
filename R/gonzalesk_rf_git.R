##***************************
## Random Forest Github
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



## Differentiate the two most common repeat_id groups? (Random forest)
#2026_11_10 Kristian Gonzales Active participation

usethis::git_remotes()
usethis::git_sitrep()
usethis::pr_push()



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

#test
count()