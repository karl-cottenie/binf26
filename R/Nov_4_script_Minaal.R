## class assignment random forest 


df_trich = read_csv("../Downloads/trichinella_spiralis_te (2).csv")

df_trich |> 
  group_by(matching_class) |> 
  count() |> 
  arrange(desc(n)) |>
  print (n = 28)

ls_twocommon <- c("LTR/Pao", "PLE/Naiad")

df_trich_twocommon <- df_trich |>
  filter(matching_class %in% ls_twocommon) #|>
  count(matching_class, name = "n") |>
  mutate(species = "T. spiralis")

df_trich_twocommon |> view()

table(df_trich_twocommon$matching_class)


# Setting seed for reproducibility
set.seed(217)

dfValidation <- df_trich_twocommon |>
  group_by(matching_class) |>
  sample_n(20)

# Checking sample size for each TE group
table(dfValidation$matching_class)


# Creating a unique identifier for each annotation
df_trich_twocommon <- df_trich_twocommon %>%
  mutate(annotation_id = row_number())

# Creating training dataset
set.seed(13)

dfTraining <- df_trich_twocommon %>%
  filter(!annotation_id %in% dfValidation$annotation_id) %>%
  group_by(matching_class) %>%
  sample_n(150)

# Checking sample size for each TE group
table(dfTraining$matching_class)

names(dfTraining)


# Building a classifier to separate LTR/Pao and PLE/Naiad
# using numerical TE annotation characteristics.

TE_classifier <- randomForest::randomForest( x = dfTraining[, c("sw_score", "perc_div", "perc_del", "perc_insert", "qry_width", "in_repeat_start", "in_repeat_end")], y = base::as.factor(dfTraining$matching_class), ntree = 50, importance = TRUE)

# Examining classifier results
TE_classifier




#### =

# Relative importance of each feature
TE_classifier$importance

# Number of times each row was left out and predicted
TE_classifier$oob.times

# Out-of-bag error rate
TE_classifier$err.rate

# Confusion matrix
TE_classifier$confusion

# Fraction of votes for each TE group
TE_classifier$votes



# Plot variable importance
varImpPlot(TE_classifier)




# Building a second classifier using only two features
TE_classifier2 <- randomForest(
  x = dfTraining[, c("sw_score", "perc_div")],
  y = base::as.factor(dfTraining$matching_class),
  ntree = 50,
  importance = TRUE
)

# Examine results
TE_classifier2






# Predicting TE groups using validation dataset
predictValidation <- predict(
  TE_classifier,
  dfValidation[, c("sw_score", "perc_div", "perc_del",
                   "perc_insert", "qry_width",
                   "in_repeat_start", "in_repeat_end")]
)

# Looking at predictions
predictValidation

# Checking prediction format
class(predictValidation)
length(predictValidation)

# Creating confusion matrix
table(
  observed = dfValidation$matching_class,
  predicted = predictValidation
)
