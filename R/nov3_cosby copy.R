##***************************
## Nov 3
##
## Taya Cosby
##
## 2026-11-03
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
# names for data frames (df_name), lists (ls_name), ... (Thanks Jacqueline May)

# differentiate DNA transposons from other repetitive elements ----
df_te_t_spir <- read_csv("../binf26/data/2026-10-26_t_spiralis copy 2.csv")
fa_t_spir <- readDNAStringSet("../binf26/data/GCA_001447595.1_T1_ISS3_r1.0.fna copy 2")

ls_non_te = c( "Simple_repeat",
               "Low_complexity", "Satellite",
               "rRNA", "snRNA", "ARTEFACT",
               "Unknown")

ls_dna_te = df_te_t_spir |> 
  filter(matching_class %notin% ls_non_te) |> 
  pull(matching_class) |> 
  unique() 
ls_dna_te

ls_dna_te = ls_dna_te[str_detect(ls_dna_te, "^DNA")]

df_te_t_spir |> 
  mutate(te_class = if_else(matching_class %in% ls_non_te, "other repetitive",
                            # nesting if_else statements
                            if_else(matching_class %in% ls_dna_te, "dna", "retro")),
         .after = matching_class) |>
  View()

ls_contig_large = fa_t_spir[width(fa_t_spir) > 1000] |> 
  names()

df_te_seq = df_te_t_spir |> 
  filter(qry_id %in% ls_contig_large) |> 
  # rows dropped from 41526 to 40518
  as.data.frame() 

df_te_seq$nucl_seq2 = DNAStringSet(df_te_seq$nucl_seq)

df_oligo4 = oligonucleotideFrequency(df_te_seq$nucl_seq2, width = 4, as.prob = T)
df_oligo4

df_oligo4[, 1:5]

df_te_t_spir |> 
  filter(qry_id %in% ls_contig_large) |> 
  bind_cols(df_oligo4[, 1:5])

set.seed(217)

dfSeq2 <- df_te_t_spir %>%
  mutate(nucleotides2 = str_remove(nucl_seq, "^[-N]+")) %>%
  mutate(nucleotides2 = str_remove(nucleotides2, "[-N]+$")) %>%
  mutate(nucleotides2 = str_remove_all(nucleotides2, "-+")) %>%
  filter(str_count(nucleotides2, "N") <= (0.05 * str_count(nucl_seq)))
dfSeq2 <- as.data.frame(dfSeq2)
dfSeq2$nucleotides2 <- DNAStringSet(dfSeq2$nucleotides2)

dfSeq2 <- cbind(dfSeq2, as.data.frame(letterFrequency(dfSeq2$nucleotides2, letters = c("A", "C","G", "T"))))

dfValidation2 <- dfSeq2 %>%
  mutate(te_class = if_else(matching_class %in% ls_non_te, "other repetitive",
                            # nesting if_else statements
                            if_else(matching_class %in% ls_dna_te, "dna", "retro")),
         .after = matching_class) %>%
  group_by(te_class) %>%
  sample_n(20)

set.seed(13)

dfTraining2 <- dfSeq2 %>%
  filter(!qry_id %in% dfValidation2$qry_id) %>%
  mutate(te_class = if_else(matching_class %in% ls_non_te, "other repetitive",
                            # nesting if_else statements
                            if_else(matching_class %in% ls_dna_te, "dna", "retro")),
         .after = matching_class) %>%
  group_by(te_class) %>%
  sample_n(62)

table(dfTraining2$te_class)

dfSeq2$Aprop <- (dfSeq2$A) / (dfSeq2$A + dfSeq2$T + dfSeq2$C + dfSeq2$G)

dfSeq2$Tprop <- (dfSeq2$T) / (dfSeq2$A + dfSeq2$T + dfSeq2$C + dfSeq2$G)

dfSeq2$Gprop <- (dfSeq2$G) / (dfSeq2$A + dfSeq2$T + dfSeq2$C + dfSeq2$G)

gene_classifier <- randomForest::randomForest(x = dfTraining2[, 86:88], y = as.factor(dfTraining2$te_class), ntree = 50, importance = TRUE)

test <- randomForest::randomForest(x = te)

?randomForest()