###### **********************
## 2026_10_29_sloot_kmer
##
## Caleb Sloot
##
## 2026-10-29
##
##***************************

rm(list = ls())

## _ Packages used -------
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(viridis)
library(Biostrings)
theme_set(theme_light())

df_te_t_spir <- read_csv("../data/2026-10-26_t_spiralis.csv") 

fa_t_spir = readDNAStringSet("../data/GCA_001447595.1_T1_ISS3_r1.0.fna")

df_te_t_spir
fa_t_spir

fa_t_spir |> 
  names() |> 
  head()

df_te_t_spir |> 
  pull(qry_id) |> 
  head()

fa_t_spir |> 
  names() |> 
  str_sub(end = 14) |> 
  head()

names(fa_t_spir) = fa_t_spir |> 
  names() |> 
  str_sub(end = 14) 

names(fa_t_spir)

df_te_t_spir |> 
  filter(qry_id %notin% names(fa_t_spir))

# all the names are the same :)
# onto the next step!

ls_width <- fa_t_spir |> 
  width() |> 
  sort(decreasing = TRUE)
  
plot(ls_width |> 
       log10())

df_te_t_spir |> 
  arrange(desc(qry_width)) 

ls_DNA_transposons <- c("DNA/TcMar-Tc1", "Dna/TcMar-Mariner", "DNA/TcMar-Pogo", "DNA/TcMar-Tc4", "DNA/PIF-Harbinger", "DNA/MULE-MuDR", "DNA/hAT-Ac", "DNA/hAT", "DNA/CMC-Chapaev")

ls_retrotransposons <-  c("LTR/Pao", "LTR/Gypsy", "LTR/Ngaro", "LTR/DIRS", "LINE/CR1", "LINE/RTE-BovB", "SINE/tRNA-Deu-L2", "SINE/tRNA-RTE", "PLE/Naiad", "PLE")

ls_non_te <- c("simple_repeat", "Simple_repeat", "Low_complexity", "Satellite", "rRNA", "snRNA", "ARTEFACT", "Unknown") 

ls_large_contigs <- names(fa_t_spir)[(width(fa_t_spir) > 10000)]

df_te_t_spir_categ <- df_te_t_spir |> 
  mutate(categories = if_else(matching_class %in% ls_non_te, "Other repetitive elements",
                              if_else(matching_class %in% ls_DNA_transposons, "DNA transposons",
                                      "Retrotransposons")))
df_te_t_spir_categ |> 
  count(categories)

df_te_large <- df_te_t_spir_categ |> 
  filter(qry_id %in% ls_large_contigs)

names(df_te_large)

## Time to calculate 4mers

df_te_4mer <- df_te_large |> 
  bind_cols(
    oligonucleotideFrequency(
      DNAStringSet(df_te_large$nucl_seq),
      width = 4,
      as.prob = TRUE) |> 
      as_tibble())

df_te_4mer |> 
  select(nucl_seq, AAAA, AAAT)

df_te_4mer_long <- df_te_4mer |> 
  pivot_longer(cols = 18:273,
              names_to = "kmer",
              values_to = "frequency")

df_te_4mer_long |> 
  count(kmer)

df_te_4mer_long |> 
  filter(kmer %in% unique(kmer)[1:5]) |> 
  ggplot(aes(x = categories, y = frequency, fill = categories)) +
  geom_boxplot() +
  facet_wrap(~ kmer) +
  labs(title = "4-mer Frequencies",
       x = "Repetitive Element Category",
       y = "Frequency") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


# Find the two most common TE matching class groups

df_te_t_spir |> 
  filter(matching_class %notin% ls_non_te) |> 
  group_by(matching_class) |> 
  summarize(count = n()) |> 
  arrange(desc(count)) |> 
  head(2)
  
# Challenge: show all the 3-mers, not just the first 5

df_oligo3 = oligonucleotideFrequency(df_te_seq$nucl_seq2, width = 3, as.prob = T)
df_oligo3

df_te_t_spir_oligo3 = df_te_t_spir |> 
  filter(qry_id %in% ls_contig_large) |> 
  bind_cols(df_oligo4)

df_te_t_spir_oligo3 |>
  filter(qry_id %in% ls_contig_large) |>
  bind_cols(df_oligo3) |>
  mutate(te_class = if_else(matching_class %in% ls_non_te, "other",
                            if_else(matching_class %in% ls_dna_te, "dna", "retro")),
         .after = matching_class) |>
  select(te_class, AAA:TTT) |> # this step not necessary but makes it easier to see the result
  pivot_longer(AAA:TTT, names_to = "kmer", values_to = "perc") |> 
  ggplot(aes(x = kmer, y = perc)) + 
  # more separation if you use 
  #  ggplot(aes(x = kmer, y = perc^0.25)) +
  geom_boxplot() + 
  facet_grid(~ te_class)



