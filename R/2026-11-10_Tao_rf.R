## Data import & filtering ---- 
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(viridis)
# + scale_color/fill_viridis(discrete = T/F)
theme_set(theme_light())
library(randomForest)
library(BiocManager)
library(Biostrings)


#Differentiate the two most common TE matching_class groups?
df_trich_te <- read_csv("../Data/2026-10-26_t_spiralis.csv")
df_trich_te$matching_class |> 
  unique() 
df_trich_te |> 
  count(matching_class) |> 
  arrange(desc(n))

summary(df_trich_te)
df_2_most_common_te <- df_trich_te %>%
  filter(!is.na(nucl_seq)) %>%
  filter(matching_class == "LTR/Pao" | matching_class == "PLE/Naiad")

unique(df_2_most_common_te$matching_class)
sum(is.na(df_2_most_common_te$nucl_seq))

dfSeq <- df_2_most_common_te %>%
  mutate(nucl_seq2 = str_remove(nucl_seq, "^[-N]+")) %>%
  mutate(nucl_seq2 = str_remove(nucl_seq2, "[-N]+$")) %>%
  mutate(nucl_seq2 = str_remove_all(nucl_seq2, "-+")) %>%
  filter(str_count(nucl_seq2, "N") <= (0.05 * str_count(nucl_seq)))

## no N in trich data so this doesn't do anything I think 
nrow(df_2_most_common_te) - nrow(dfSeq)

q1_ltr <- quantile(nchar(dfSeq$nucl_seq2[dfSeq$matching_class == "LTR/Pao"]), probs = 0.25, na.rm = TRUE)
q1_ltr
q3_ltr <- quantile(nchar(dfSeq$nucl_seq2[dfSeq$matching_class == "LTR/Pao"]), probs = 0.75, na.rm = TRUE)
q3_ltr

q1_ple <- quantile(nchar(dfSeq$nucl_seq2[dfSeq$matching_class == "PLE/Naiad"]), probs = 0.25, na.rm = TRUE)
q1_ple
q3_ple <- quantile(nchar(dfSeq$nucl_seq2[dfSeq$matching_class == "PLE/Naiad"]), probs = 0.75, na.rm = TRUE)
q3_ple

dfSeq <- dfSeq %>%
  filter(((str_count(nucl_seq2) >= q1_ltr & str_count(nucl_seq2) <= q3_ltr & matching_class == "LTR/Pao") | (str_count(nucl_seq2) >= q1_ple & str_count(nucl_seq2) <= q3_ple & matching_class == "PLE/Naiad"))) 

dim(dfSeq)
unique(dfSeq$matching_class)
table(dfSeq$matching_class)
sum(is.na(dfSeq$nucl_seq2))
summary(str_count(dfSeq$nucl_seq2[dfSeq$matching_class == "LTR/Pao"]))
summary(str_count(dfSeq$nucl_seq2[dfSeq$matching_class == "PLE/Naiad"]))

rm(dfSeqCompare)

## Sequence features ---- 
dfSeq <- as.data.frame(dfSeq)
dfSeq$nucl_seq2 <- DNAStringSet(dfSeq$nucl_seq2)
class(dfSeq$nucl_seq2)

dfSeq <- cbind(dfSeq, as.data.frame(letterFrequency(dfSeq$nucl_seq2, letters = c("A", "C","G", "T"))))
dfSeq

dfSeq$Aprop <- (dfSeq$A) / (dfSeq$A + dfSeq$T + dfSeq$C + dfSeq$G)

dfSeq$Tprop <- (dfSeq$T) / (dfSeq$A + dfSeq$T + dfSeq$C + dfSeq$G)

dfSeq$Gprop <- (dfSeq$G) / (dfSeq$A + dfSeq$T + dfSeq$C + dfSeq$G)

dfSeq <- cbind(dfSeq, as.data.frame(dinucleotideFrequency(dfSeq$nucl_seq2, as.prob = TRUE)))

dfSeq <- cbind(dfSeq, as.data.frame(trinucleotideFrequency(dfSeq$nucl_seq2, as.prob = TRUE)))
dfSeq

#Training random forest model

dfSeq <- dfSeq %>%
  mutate(qry_id = row_number())
table(dfSeq$matching_class)


dfSeq$nucl_seq2 <- as.character(dfSeq$nucl_seq2)
table(dfSeq$matching_class)

set.seed(217)
dfValidation <- dfSeq %>%
  group_by(matching_class) %>%
  sample_n(20)

set.seed(13)
dfTraining <- dfSeq %>%
  filter(!qry_id %in% dfValidation$qry_id) %>%
  group_by(matching_class) %>%
  sample_n(88)

dfSeq
names(dfTraining)

gene_classifier <- randomForest::randomForest(x = dfTraining[, 22:24], y = as.factor(dfTraining$matching_class), ntree = 50, importance = TRUE)
gene_classifier

## Checking prediction accuracy ---- 
predictValidation <- predict(gene_classifier, dfValidation[, c(11, 22:24)])
predictValidation
class(predictValidation)
length(predictValidation)
table(observed = dfValidation$matching_class, predicted = predictValidation)




