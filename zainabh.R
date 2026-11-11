##***************************
## kmer analysis + random forest
##
## Hala zainab
##
## 2026-10-29
##
##***************************

## _ Packages used -------
rm(list = ls())
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(Biostrings)
library(viridis)
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

## load in data -------

fa_t_spir = readDNAStringSet("../data/GCA_001447595.1_T1_ISS3_r1.0.fna")

fa_t_spir

df_te_t_spir = read_csv("../data/2026-10-26_t_spiralis.csv")

df_te_t_spir

## cleaning up two data sets -----

## _ first issue, from the instructions, different identifying names for contigs
fa_t_spir |> names() |> head()

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

fa_t_spir
names(fa_t_spir)

df_te_t_spir |> 
  filter(qry_id %notin% names(fa_t_spir))
# ===> yay!!!! all the names are the same

## _ second issue, from Livpreet's question, values in between parentheses ----

df_te_t_spir |> View()
#====> issue in 2 columns: qry_left and in_repeat_end

# contrary to what I said in class, you should fix the problem here!!!

df_te_t_spir |> 
  mutate(qry_left2 = as.numeric(str_replace_all(qry_left, "\\(|\\)", "")),
         # the \\( is necessary because ( is a special character))
         .after = qry_left) |> 
  mutate(in_repeat_end2 = as.numeric(str_replace_all(in_repeat_end, "\\(|\\)", "")),
         .after = in_repeat_end) |> View()

#====> depending on downstream analyses, you can either 
# overwrite the original column, and then save the dataframe
df_te_t_spir =  df_te_t_spir |> 
  mutate(qry_left = as.numeric(str_replace_all(qry_left, "\\(|\\)", "")),
         # notice the ary_left instead of qry_left2 from line 73
         .after = qry_left) |> 
  mutate(in_repeat_end = as.numeric(str_replace_all(in_repeat_end, "\\(|\\)", "")),
         .after = in_repeat_end)
# WARNING: the () indicate if the TE is on the + or C direction (matching_repeat column), and this might influence what you want to do later in your analysis!!! So use this option with caution

# the other option is to save the code from line 72 as as

## large contigs ------
ls_width = fa_t_spir |> 
  width() 

ls_width |> 
  plot()

ls_width |> log10() |> 
  plot()
#====> inflection point around 1000 bp
# anything bigger than 1000 is "large"

ls_contig_large = fa_t_spir[width(fa_t_spir) > 1000] |> 
  names()

## creating the ggplot from the assignment ------

## _filering large contigs -----
df_te_seq = df_te_t_spir |> 
  filter(qry_id %in% ls_contig_large) |> 
  # rows dropped from 41526 to 40518
  as.data.frame() 

## _obtain 4-mer frequencies -----
df_te_seq$nucl_seq2 = DNAStringSet(df_te_seq$nucl_seq)

df_oligo4 = oligonucleotideFrequency(df_te_seq$nucl_seq2, width = 4, as.prob = T)
df_oligo4

df_oligo4[, 1:5]

df_te_t_spir |> 
  filter(qry_id %in% ls_contig_large) |> 
  bind_cols(df_oligo4[, 1:5]) 

## _ creating a new variable with 3 categories -----

# reminder, from the tidyverse class
df_te_t_spir |> 
  pull(matching_class) |> 
  unique()
# copy-paste to Copilot and ask for repetitive elements that are not DNA or retrotransposons

ls_non_te = c( "Simple_repeat",
               "Low_complexity", "Satellite",
               "rRNA", "snRNA", "ARTEFACT",
               "Unknown")
ls_non_te

# from a later tidyverse class
ls_dna_te = df_te_t_spir |> 
  filter(matching_class %notin% ls_non_te) |> 
  pull(matching_class) |> 
  unique() 
ls_dna_te

ls_dna_te = ls_dna_te[str_detect(ls_dna_te, "^DNA")]

# creating the new variable
df_te_t_spir <- df_te_t_spir |> 
  mutate(te_class = if_else(matching_class %in% ls_non_te, "other",
                            # nesting if_else statements
                            if_else(matching_class %in% ls_dna_te, "dna", "retro")),
         .after = matching_class) 
# .after so that I can easily compare whether the new variable is created correctly

## _ putting it all together -------

df_te_t_spir |> 
  filter(qry_id %in% ls_contig_large) |> 
  bind_cols(df_oligo4[, 1:5]) |> 
  mutate(te_class = if_else(matching_class %in% ls_non_te, "other",
                            # nesting if_else statements
                            if_else(matching_class %in% ls_dna_te, "dna", "retro")),
         .after = matching_class) |> 
  select(te_class, AAAA:AACA) |> # this step not necessary but makes it easier to see the result
  pivot_longer(AAAA:AACA, names_to = "kmer", values_to = "perc") |> 
  ggplot(aes(x = kmer, y = perc)) + 
  # more separation if you use 
  #  ggplot(aes(x = kmer, y = perc^0.25)) +
  geom_boxplot() + 
  facet_grid(~ te_class)

## _ other way of doing it -------

# all in one go
# but difficult to get this done as you are coding
# the previous code has a lot more checks and balances
# and I changed the ggplot a little to provide another example of ggplot

df_te_t_spir |> 
  mutate(te_class = if_else(matching_class %in% ls_non_te, "other",
                            # nesting if_else statements
                            if_else(matching_class %in% ls_dna_te, "dna", "retro")),
         .after = matching_class) |> 
  bind_cols(oligonucleotideFrequency(DNAStringSet(df_te_t_spir$nucl_seq),
                                     width = 4, as.prob = T)[,1:5]) |> 
  filter(qry_id %in% ls_contig_large) |> 
  pivot_longer(AAAA:AACA, names_to = "kmer", values_to = "perc") |> 
  ggplot(aes(x = te_class, y = perc)) + 
  geom_boxplot() + 
  facet_grid(kmer ~ .)

## new challenges -----

### _ show all the 3-mers, not just the first 5 -------

df_te_seq$nucl_seq2 = DNAStringSet(df_te_seq$nucl_seq)

df_oligo3 = oligonucleotideFrequency(
  df_te_seq$nucl_seq2,
  width = 3,
  as.prob = TRUE
)

df_oligo3

df_te_t_spir |>
  filter(qry_id %in% ls_contig_large) |>
  bind_cols(df_oligo3) |>
  mutate(
    te_class = if_else(
      matching_class %in% ls_non_te, "other",
      if_else(matching_class %in% ls_dna_te, "dna", "retro")
    ),
    .after = matching_class
  ) |>
  select(te_class, all_of(colnames(df_oligo3))) |>
  pivot_longer(
    -te_class,
    names_to = "kmer",
    values_to = "perc"
  ) |>
  ggplot(aes(x = kmer, y = perc)) +
  geom_boxplot() +
  facet_grid(~ te_class) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

df_te_t_spir |>
  filter(qry_id %in% ls_contig_large) |>
  bind_cols(df_oligo3) |>
  mutate(
    te_class = if_else(
      matching_class %in% ls_non_te, "other",
      if_else(matching_class %in% ls_dna_te, "dna", "retro")
    )
  ) |>
  pivot_longer(
    cols = all_of(colnames(df_oligo3)),
    names_to = "kmer",
    values_to = "perc"
  ) |>
  ggplot(aes(x = te_class, y = perc)) +
  geom_boxplot() +
  facet_wrap(~ kmer, ncol = 8)


### -

## Can randomforest and k-mer approach differentiate the two most common TE matching_class groups?

df_te_t_spir |>
  count(matching_class, sort = TRUE)


top2_classes <- df_te_t_spir |>
  count(matching_class, sort = TRUE) |>
  dplyr::slice(1:2) |>
  pull(matching_class)

top2_classes

df_rf <- df_te_t_spir |>
  filter(qry_id %in% ls_contig_large) |>
  filter(matching_class %in% top2_classes)

dna_seqs <- DNAStringSet(df_rf$nucl_seq)

df_kmer <- oligonucleotideFrequency(
  dna_seqs,
  width = 3,
  as.prob = TRUE
)

df_rf <- bind_cols(
  df_rf,
  as.data.frame(df_kmer)
)

class(dna_seqs)
dim(df_kmer)
colnames(df_kmer)


library(randomForest)

rf_te <- randomForest(
  x = df_kmer,
  y = base::as.factor(df_rf$matching_class),
  ntree = 500,
  importance = TRUE
)

rf_te

## answer: Can a random forest using k-mer frequencies differentiate the two most common TE matching_class groups in Trichinella spiralis?

# Yes.

#The random forest classifier achieved an OOB error rate of 6.18%, which corresponds to an estimated classification accuracy of ~93.8. This means that the two most common transposable element classes have distinct trinucleotide sequence composition patterns that can be detected using 3-mer frequencies alone

