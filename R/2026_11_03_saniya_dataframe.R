###### **********************
## CLASS 9 - SOFTWARE TOOLS - SEQUENCE CLASSIFICATION AND RANDOM FOREST V1
##
## Saniya Chadha
##
## 2026-11-03
##
##***************************

# Can a Random Forest using base proportions differentiate wasp taxonomy at the genus level (Cotesia vs. Campoletis)?

## _ Packages used -------
library(tidyverse)    # For data wrangling and tables
library(viridis)      # For color palettes
library(randomForest) # For machine learning classification
library(Biostrings)   # For biological string operations

# Set conflict preference for filter
conflicted::conflict_prefer("filter", "dplyr")


##_ Data used -------
dfBOLD <- read_tsv("../data/wasps.tsv")


##_filter data set for COI-5p -----
dfCOI28S <- dfBOLD %>%
  filter(!is.na(nucleotides)) %>% 
  # Code throws away rows that are missing DNA sequence data 
  filter(markercode == "COI-5P")


##_trim the data set now to remove gaps, ambiguous bases (N), internal gaps ------
dfSeq <- dfCOI28S %>%
  mutate(nucleotides2 = str_remove(nucleotides, "^[-N]+")) %>%
  # Removes leading gaps (-) or ambiguous bases (N) from start (^) of DNA sequence
  mutate(nucleotides2 = str_remove(nucleotides2, "[-N]+$")) %>%
  # Removes the same from the end of the sequence
  mutate(nucleotides2 = str_remove_all(nucleotides2, "-+")) %>%
  # Removes internal dash alignment gaps (-) converting it into a pure string 
  filter(str_count(nucleotides2, "N") <= (0.05 * str_count(nucleotides))) 
# Count how many ambiguous N letters are left. If more than 5% contain Ns, it is removed.


##_calculate the length variability -----
q1 <- quantile(nchar(dfSeq$nucleotides2), probs = 0.25, na.rm = TRUE)
q1 # q1 is the 25th percentile (short end), q3 is the 75th percentile (long end)
q3 <- quantile(nchar(dfSeq$nucleotides2), probs = 0.75, na.rm = TRUE)
q3


##_only keep sequences that are in the middle 50% length range ----
dfSeq <- dfSeq %>%
  filter(str_count(nucleotides2) >= q1 & str_count(nucleotides2) <= q3)


##_convert to DNAstringset for biostrings package ------
dfSeq <- as.data.frame(dfSeq)
dfSeq$nucleotides2 <- DNAStringSet(dfSeq$nucleotides2)
class(dfSeq$nucleotides2)


##_ calculate the absolute frequencies of each letter ---- 
dfSeq <- cbind(dfSeq, as.data.frame(letterFrequency(dfSeq$nucleotides2, letters = c("A", "C","G", "T"))))


##_ convert counts to proportions ---- 
dfSeq$Aprop <- (dfSeq$A) / (dfSeq$A + dfSeq$T + dfSeq$C + dfSeq$G)
dfSeq$Tprop <- (dfSeq$T) / (dfSeq$A + dfSeq$T + dfSeq$C + dfSeq$G)
dfSeq$Gprop <- (dfSeq$G) / (dfSeq$A + dfSeq$T + dfSeq$C + dfSeq$G)

dfSeq$nucleotides2 <- as.character(dfSeq$nucleotides2)

##_ filter for comparison of cotesia and campoletis ------
dfTaxonomy <- dfSeq %>%
  filter(genus_name == "Cotesia" | genus_name == "Campoletis")


##_ convert genus name back to character data type-----
dfTaxonomy$genus_name <- as.character(dfTaxonomy$genus_name)


##_ view -------
table(dfTaxonomy$genus_name) |> view()


##_ samples to be set aside for later when validating ----- 
set.seed(1104)
dfValidation <- dfTaxonomy %>%
  group_by(genus_name) %>%
  sample_n(20)


##_ create a training set for data -----
set.seed(1104)
dfTraining <- dfTaxonomy %>%
  filter(!processid %in% dfValidation$processid) %>%
  group_by(genus_name) %>%
  sample_n(50)


##_ train the model --- 
# Changed index positions 86:88 to explicit text names to protect against structural column shifts
gene_classifier <- randomForest::randomForest(
  x = dfTraining[, c("Aprop", "Tprop", "Gprop")], 
  y =  base::as.factor(dfTraining$genus_name), #use specifc base model
  ntree = 50, 
  importance = TRUE
)

##_ view model performance summary ---
gene_classifier

##_ validate the model on unseen testing data ---
# Fixed variable reference mapping errors to evaluate dfValidation properly
predictValidation <- predict(gene_classifier, dfValidation[, c("Aprop", "Tprop", "Gprop")])

final_table <- table(observed = dfValidation$genus_name, predicted = predictValidation)
final_table

#not sure if anything is right here
