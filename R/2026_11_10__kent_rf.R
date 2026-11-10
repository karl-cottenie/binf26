##***************************
## random forest cont.

##
## Mikayla Kent
##
## 2026-11-10
##
##***************************

## _ Packages used -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(viridis)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())
library(Biostrings)
library(randomForest)

# Startup ends here

## _ Comment codes ------
# Coding explanations (#, often after the code, but not exclusively)
# Code organization (## XXXXX -----)
# Justification for a section of code ## XXX
# Dead end analyses because it did not work, or not pursuing this line of inquiry (but leave it in as a trace of it, to potentially solve this issue, or avoid making the same mistake in the future # (>_<) 
# Solutions/results/interpretations (#==> XXX)
# Reference to manuscript pieces, figures, results, tables, ... # (*_*)
# TODO items #TODO
# names for data frames (df_name), lists (ls_name), ...
##load data ----

t_spiralis <- read_csv("../Data/2026-10-26_t_spiralis2.csv")

## differentiating classes of TEs ----

ls_transposons <- c("DNA/TcMar-Tc4", "DNA/TcMar-Tc1", "DNA/hAT", "DNA/PIF-Harbinger", "DNA/TcMar-Mariner", "RC/Helitron", "DNA/hAT-Ac", "DNA/MULE-MuDR", "DNA/TcMar-Pogo", "DNA/CMC-Chapaev")

ls_retrotransposons <- c("PLE/Naiad", "PLE", "LINE/CR1", "LTR/Pao", "LTR/Gypsy", "LINE/RTE-BovB", "SINE/tRNA-RTE", "LTR/Ngaro", "SINE/tRNA-Deu-L2", "LTR/DIRS")

##adding TE type info to dataframe
df_te_t_spir <- t_spiralis |> 
  mutate(te_class = if_else(matching_class %in% ls_retrotransposons, "retro",
                            if_else(matching_class %in% ls_transposons, "dna", "other")),
         .after = matching_class) 

##remove transposons
df_t_spir_ro <- df_te_t_spir |> 
  filter(te_class %notin% ls_transposons)

##training model ----
#need some sort of predictor value
#kmer frequencies?

tri <- oligonucleotideFrequency(
  DNAStringSet(df_t_spir_ro$nucl_seq),
  width = 3)


#putting timers with te class data
df_te_tri <- bind_cols(
  df_t_spir_ro["te_class"],
  as.data.frame(tri))


#split into train/test

set.seed(120)

train_rows <- sample(
  seq_len(nrow(df_te_tri)),
  size = 0.8 * nrow(df_te_tri)
)

train <- df_te_tri[train_rows, ]
test  <- df_te_tri[-train_rows, ]
#does this ix my no identifyer problem?


#pull request


