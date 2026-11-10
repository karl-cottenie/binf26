## _ random forest challenge: retro vs other -------

library(randomForest)


## keep only retrotransposons and other repetitive elements -------

df_rf <- df_te_all4 |>
  filter(te_class == "retro" | te_class == "other")

# check classes
unique(df_rf$te_class)

table(df_rf$te_class)


## create a row ID -------

df_rf <- df_rf |>
  mutate(row_id = row_number())

## create validation dataset -------

set.seed(217)

dfValidation <- df_rf |>
  group_by(te_class) |>
  sample_n(200)

# check validation sample size
table(dfValidation$te_class)


## create training dataset -------

set.seed(13)

dfTraining <- df_rf |>
  filter(!row_id %in% dfValidation$row_id) |>
  group_by(te_class) |>
  sample_n(600)

# check training sample size
table(dfTraining$te_class)


## identify the 4-mer columns -------

which(names(dfTraining) == "AAAA")

which(names(dfTraining) == "TTTT")