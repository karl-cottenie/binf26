

##make randomforest that differentiates DNA transposons from other types of repetitive elements

## load in data -------

df_te_t_spir = read_csv("../data/2026-10-26_t_spiralis.csv")

matching_unique <- df_te_t_spir |> pull(matching_class) |> unique() 

ls_dna_te <- matching_unique[str_detect(matching_unique, "^DNA")]




####4 - CALCULATING SEQUENCE FEATURES----
df_te_t_spir <- as.data.frame(df_te_t_spir)
df_te_t_spir$nucl_seq_1 <- DNAStringSet(df_te_t_spir$nucl_seq)

class(df_te_t_spir$nucl_seq_1)


#Calculating the nucleotide frequencies and appending onto our dataframe using cbind(). Here, we are using letterFrequency, rather than alphabetFrequency, as we don't want the Ns (remember we may still have up to 5% Ns in a given sequence, given our above filtering steps), and other uncertain characters included at this stage. (We could use a more stringent filter on Ns if we have a large sample size, such as in a sensitivity analysis or future study.)
df_te_t_spir <- cbind(df_te_t_spir, as.data.frame(letterFrequency(df_te_t_spir$nucl_seq_1, letters = c("A", "C","G", "T"))))

#Look at dfSeq in viewer. Verify that the new columns are appended to the end of our data frame.

#Adding A, T, and G proportions in relation to total nucleotides (C would then be the remainder and so isn't an independent variable once we know A, T, and G) as these are helpful features less influenced by sequence length variation than are counts.
df_te_t_spir$Aprop <- (df_te_t_spir$A) / (df_te_t_spir$A + df_te_t_spir$T + df_te_t_spir$C + df_te_t_spir$G)

df_te_t_spir$Tprop <- (df_te_t_spir$T) / (df_te_t_spir$A + df_te_t_spir$T + df_te_t_spir$C + df_te_t_spir$G)

df_te_t_spir$Gprop <- (df_te_t_spir$G) / (df_te_t_spir$A + df_te_t_spir$T + df_te_t_spir$C + df_te_t_spir$G)

#After each batch of additions, check out df_te_t_spir in the viewer. Look at the new columns appended on the end, and check to ensure that is what we want. Note that above we were doing basic arithmetic (vectorized fashion R style!).

#Adding dinucleotide frequency (k-mers of length 2), here using proportions as that let's us account for variability in sequence lengths.
df_te_t_spir <- cbind(df_te_t_spir, as.data.frame(dinucleotideFrequency(df_te_t_spir$nucl_seq_1, as.prob = TRUE)))

#Let's check the viewer again.

#Adding trinucleotide frequency (k-mers of length 3)
df_te_t_spir <- cbind(df_te_t_spir, as.data.frame(trinucleotideFrequency(df_te_t_spir$nucl_seq_1, as.prob = TRUE)))

#turn DNAStringSet back into vector
df_te_t_spir$nucl_seq_1 <- as.vector(df_te_t_spir$nucl_seq_1)

# creating the new variable
df_te_t_spir <- df_te_t_spir |> 
  mutate(te_class = if_else(matching_class %in% ls_dna_te, "dna", "other"),
         .after = matching_class)  |> mutate(unique_id = row.names(df_te_t_spir))

#identify 25% of sample for validation
table(df_te_t_spir$te_class) #dna = 120, other = 41406

120*0.25


set.seed(1)

dfValidation <- df_te_t_spir %>%
  group_by(te_class) %>%
  sample_n(30)

#identify number of samples for training

set.seed(2)
dfTraining <- df_te_t_spir %>%
  filter(!unique_id %in% dfValidation$unique_id) %>%
  group_by(te_class) %>%
  sample_n(90)

colnames(df_te_t_spir[,25:40])
         
te_dna_classifier <- randomForest::randomForest(x = dfTraining[, 25:40], y = base::as.factor(dfTraining$te_class), ntree = 50, importance = TRUE)

#Let's look at results.
te_dna_classifier
#Excellent performance!

#We can also dive into the specifics, such as by looking at the following. If we type the name of our classifier and $, then we can use arrows to scroll through the available options. See the documentation for the randomForest() function for more detailed information.

#relative importance of each feature
te_dna_classifier$importance

#number of times each row in the data set was left "out of bag" and predicted
te_dna_classifier$oob.times

#out of bag error rate based upon all prior trees built up to that iteration. (i.e. error rate for observations left out of that iteration and predicted based upon all decision trees built up to that point). Here, we can see we had no errors from the start!
te_dna_classifier$err.rate

#confusion matrix. (true positives, true negatives, false positives, and false negatives)
te_dna_classifier$confusion

#for each row, fraction of votes from the random forest for each class. In this example, there is agreement. (Contrast this with the outcome for a more difficult classification problem.)
te_dna_classifier$votes







