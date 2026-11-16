library(tidyverse)
library(randomForest)
library(Biostrings)

fa_t_spir = readDNAStringSet("../data/GCA_001447595.1_T1_ISS3_r1.0.fna")

df_te_t_spir = read_csv("../data/2026-10-26_t_spiralis.csv")


df_te_t_spir |> 
  mutate(qry_left2 = as.numeric(str_replace_all(qry_left, "\\(|\\)", "")),
         .after = qry_left) |> 
  mutate(in_repeat_end2 = as.numeric(str_replace_all(in_repeat_end, "\\(|\\)", "")),
         .after = in_repeat_end) |> View()
ls_contig_large = fa_t_spir[width(fa_t_spir) > 1000] |> 
  names()

df_te_seq = df_te_t_spir |> 
  filter(qry_id %in% ls_contig_large) |> 
  as.data.frame() 

df_te_seq$nucl_seq2 = DNAStringSet(df_te_seq$nucl_seq)

df_oligo4 = oligonucleotideFrequency(df_te_seq$nucl_seq2, width = 4, as.prob = T)
df_oligo4

df_oligo4[, 1:5]

f_te_t_spir |> 
  filter(qry_id %in% ls_contig_large) |> 
  bind_cols(df_oligo4[, 1:5])

df_te_t_spir |> 
  pull(matching_class) |> 
  unique()

ls_non_te = c( "Simple_repeat",
               "Low_complexity", "Satellite",
               "rRNA", "snRNA", "ARTEFACT",
               "Unknown")

ls_dna_te = df_te_t_spir |> 
  filter(matching_class %notin% ls_non_te) |> 
  pull(matching_class) |> 
  unique() 

ls_dna_te = ls_dna_te[str_detect(ls_dna_te, "^DNA")]

#k-mer from last class
df_te_t_spir |> 
  filter(qry_id %in% ls_contig_large) |> 
  bind_cols(df_oligo4[, 1:5]) |>
  mutate(te_class = if_else(matching_class %in% ls_non_te, "other",
                            if_else(matching_class %in% ls_dna_te, "dna", "retro")),
         .after = matching_class) |>
  select(te_class, AAAA:AACA) |>
  pivot_longer(AAAA:AACA, names_to = "kmer", values_to = "perc") |>
  ggplot(aes(x = kmer, y = perc)) + 
  geom_boxplot() + 
  facet_grid(~ te_class)


#randomForest class activity ----

df_te_t_spir |> 
  filter(qry_id %in% ls_contig_large) |> 
  bind_cols(df_oligo4[, 1:5]) |>
  mutate(te_class = if_else(matching_class %in% ls_non_te, "other",
                            if_else(matching_class %in% ls_dna_te, "dna", "retro")),
         .after = matching_class) |>
  select(te_class, AAAA:AACA) |>
  pivot_longer(AAAA:AACA, names_to = "kmer", values_to = "perc") |> view()
