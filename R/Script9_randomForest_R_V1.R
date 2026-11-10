###### **********************
## CLASS 9 - SOFTWARE TOOLS - SEQUENCE CLASSIFICATION AND RANDOM FOREST V1
##
## Karl Cottenie
##
## 2024-09-29
##
##***************************

## _ Packages used -------
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(viridis)
# + scale_color/fill_viridis(discrete = T/F)
theme_set(theme_light())

#Script adapted from Sally Adamowicz, with consultation of below sources; last updated October 6, 2021. V1 poses challenge in Part 6. See V2 for example answer.

##SOURCES: During the preparation of this tutorial, I consulted the randomForest R package documentation as well as looked at the following tutorial by Saulo Pires de Oliveira, who uses the built-in iris dataset to show a simple example of using the randomForest package. I suggest you may also find it helpful to go through that tutorial.
#https://www.blopig.com/blog/2017/04/a-very-basic-introduction-to-random-forests-using-r/

#I suggest this posting is also helpful for understanding these concepts:
#https://towardsdatascience.com/understanding-random-forest-58381e0602d2

#In this script, we use the functions from the randomForest R package and apply them to DNA sequence classification, using publicly available data from BOLD.

####1 - OVERVIEW----

#Introduction: Recently, we have seen how we can extract simple sequence features using functions available through the Biostrings package. Today, we will explore one example of how such features can be used.

#Classification is a broad category of problem in diverse settings in machine learning and statistics. In biology, we often want to be able to predict what group a sample will belong to. An example would be wanting to assign individual human patients to disease risk groups on the basis of their biomarkers. Another major application of classification methods in bioinformatics is the assignment of biological sequences to taxonomic groups or genes.

#In the below example, we will use simple sequence features to classify sequences to one of two genes and to identify sequences to taxonomic groups.

#After this tutorial, try to build your own classifier for a problem of interest to you!

####2 - INSTALLING AND LOADING PACKAGES----

#Packages from CRAN:

#install tidyverse suite of packages, if needed, then load.
#install.packages("tidyverse")
library(tidyverse)

#install randomForest package, if needed, and load library.
#install.packages("randomForest")
library(randomForest)

#Package from Bioconductor:

#These are the lines to install Bioconductor packages for recent versions of R. Then, load library.
#install.packages("BiocManager")
#library(BiocManager)
#BiocManager::install("Biostrings")
library(Biostrings)

####3 - LOADING, CHECKING PARSING, AND FILTERING DATA----

#This is the code I used to obtain the dataset analyzed in this script. I ran the below line on October 3, 2021. We are expanding upon our "Cotesia" test dataset that we looked at last time (that genus is from the wasp family Braconidae) by also downloading a second genus within a different family of wasps (family Ichneumonidae). Both of the genera being analyzed belong to the order Hymenoptera (which consists of wasps, bees, ants, and allies).
#wasps <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Cotesia|Campoletis&format=tsv")

#This is the code I used to write to disk.
#write_tsv(wasps, "wasps.tsv")

#Now, we can read the provided tsv back in so that we all have exactly the same file.
dfBOLD <- read_tsv("../data/wasps.tsv")

#You should now see in your global environment a df with 6,542 observations of 80 variables.

#summary
summary(dfBOLD)

#check column names
names(dfBOLD)

#table will give us a count of records per marker code
table(dfBOLD$markercode)

#Or, we could use a tidyverse block of code as an alternative. This version will sort and assign the results to a tibble, which we could examine further.
dfBOLD.by.marker <- dfBOLD %>%
  filter(!is.na(nucleotides)) %>%
  count(markercode) %>%
  arrange(desc(n)) %>%
  print()

#remove - if not needed for downstream analysis.
rm(dfBOLD.by.marker)

#Downstream, we want to look at sequence features. These may be affected by sequence length and the portion of a given gene that is available. So, we will first have a look at sequence lengths and will consider creating a subset.

#histogram and summary of sequence lengths for COI. We can see a wide range of sequence lengths. Note that here we are using nchar() function to count the number of characters in each sequence (the sequences are strings). Here, we are summarizing sequence length only for the COI-5P marker.
hist(nchar(dfBOLD$nucleotides[dfBOLD$markercode == "COI-5P"]), xlab = "Sequence Length", ylab = "
     Frequency", main = "Frequency Histogram of COI Sequence Lengths")

summary(nchar(dfBOLD$nucleotides[dfBOLD$markercode == "COI-5P"]))

#Same for 28S
hist(nchar(dfBOLD$nucleotides[dfBOLD$markercode == "28S"]), xlab = "Sequence Length", ylab = "
     Frequency", main = "Frequency Histogram of 28S Sequence Lengths")

summary(nchar(dfBOLD$nucleotides[dfBOLD$markercode == "28S"]))

#We are filtering and trimming the data for further analysis. After filtering out records with no sequence, next we filter to retain for the COI-5P and 28S markers. Note the usage of the vertical pipe symbol "|" to indicate "or". We want all rows that are either COI or 28S.
dfCOI28S <- dfBOLD %>%
  filter(!is.na(nucleotides)) %>%
  filter(markercode == "COI-5P" | markercode == "28S")

#performing a couple of checks to verify that our filtering step above worked as intended. We are expecting just COI-5P and 28S data remaining. Also, remember that is.na() yields a logical vector (try out that part alone). TRUE = 1 and FALSE = 0. So, we can use sum to check if we have any records left having missing sequences; should be zero.
unique(dfCOI28S$markercode)
sum(is.na(dfCOI28S$nucleotides))

#We are next going to create another dataset with a filter for sequence length as well and to clean up the sequences. We will trim the ends, removing gaps and Ns from the ends. We'll also remove gaps throughout (as we want biological sequences, not alignments, here). And, we'll filter out sequences with a lot of missing data (> 5% Ns, after the end trimming for Ns and gaps). Variability in sequence length can be due to different primer sets being used, or to poor quality sequences. Long sequences would be due to the whole gene (rather than partial gene) being available, but most of these sequences are partial gene sequences. 28S, being a non-coding rRNA gene, can be expected to have some natural length variability due to indels (insertions and deletions during evolutionary history). Constraining the length variability would give us more consistent and comparable k-mer frequencies. Specific choices for such filtering steps would vary depending upon the goals of the analysis. And, exploring results across a range of parameters (i.e. sensitivity analysis) can also be helpful. I gave the df the name dfSeq, so that we have a short name to work with downstream. You may choose to use longer names if you wish to retain more information about the contents of a particular data frame. Note the usage of & (and) and | (or), which are used in a vectorized fashion.

#We're creating a new nucleotides column so that we can compare the old and new nucleotides columns.

dfSeq <- dfCOI28S %>%
    mutate(nucleotides2 = str_remove(nucleotides, "^[-N]+")) %>%
  mutate(nucleotides2 = str_remove(nucleotides2, "[-N]+$")) %>%
  mutate(nucleotides2 = str_remove_all(nucleotides2, "-+")) %>%
  filter(str_count(nucleotides2, "N") <= (0.05 * str_count(nucleotides)))

#Run the below, and have a look in the viewer to compare the old and new sequence columns.
dfSeqCompare <- cbind(dfSeq$nucleotides, dfSeq$nucleotides2)

#Next, we are assigning a variable name to hold the information about the first quartile sequence length for COI; i.e. at what value of sequence length do 25% of the observations fall below that sequence length and 75% of sequence lengths fall above that? (Note, here we don't need to specify na.rm = TRUE, as we already removed NAs, but I left that in as a reminder that often that argument is needed for this function.)
q1 <- quantile(nchar(dfSeq$nucleotides2[dfSeq$markercode == "COI-5P"]), probs = 0.25, na.rm = TRUE)
q1

q3 <- quantile(nchar(dfSeq$nucleotides2[dfSeq$markercode == "COI-5P"]), probs = 0.75, na.rm = TRUE)
q3

#We could have embedded this information in the filtering block above, but this can get a little hard to read. So, I suggest to consider breaking up each set of actions into reasonable and readable sized lines or blocks of code.

#Continuing with the filtering ... now, we can use our calculated q1 and q3 instead of hard coding in the sequence lengths. Below, we are filtering for sequence length, keeping COI sequences with length greater than or equal to the first quartile length value and also sequences with length equal to or less than the third quartile value. This is constraining length variability. In particular, very short sequences could have different properties compared to the whole gene region. We are keeping all 28S values here, as sequence length was more constrained, and there were fewer records for that marker.
dfSeq <- dfSeq %>%
  filter(((str_count(nucleotides2) >= q1 & str_count(nucleotides2) <= q3 & markercode == "COI-5P") | markercode == "28S"))

#some checks to make sure everything worked as expected. 
dim(dfSeq)
unique(dfSeq$markercode)
table(dfSeq$markercode)
sum(is.na(dfSeq$nucleotides2))
summary(str_count(dfSeq$nucleotides2[dfSeq$markercode == "COI-5P"]))
summary(str_count(dfSeq$nucleotides2[dfSeq$markercode == "28S"]))

#These checks pan out. We could add more (such as to check that the beginnings and ends of the sequences have been trimmed of Ns and gaps). Make the prediction about what should be the outcome before looking at the results.

#removing objects not needed downstream
rm(dfSeqCompare)

####4 - CALCULATING SEQUENCE FEATURES----

#Next, we are converting the nucleotides to a DNAStringSet (class) so that we can use functions from the Biostrings package. Note that this time we are editing in place (using same column name before and after assignment operator), as we have already checked out this function before.
dfSeq <- as.data.frame(dfSeq)
dfSeq$nucleotides2 <- DNAStringSet(dfSeq$nucleotides2)
class(dfSeq$nucleotides2)

#Calculating the nucleotide frequencies and appending onto our dataframe using cbind(). Here, we are using letterFrequency, rather than alphabetFrequency, as we don't want the Ns (remember we may still have up to 5% Ns in a given sequence, given our above filtering steps), and other uncertain characters included at this stage. (We could use a more stringent filter on Ns if we have a large sample size, such as in a sensitivity analysis or future study.)
dfSeq <- cbind(dfSeq, as.data.frame(letterFrequency(dfSeq$nucleotides2, letters = c("A", "C","G", "T"))))

#Look at dfSeq in viewer. Verify that the new columns are appended to the end of our data frame.

#Adding A, T, and G proportions in relation to total nucleotides (C would then be the remainder and so isn't an independent variable once we know A, T, and G) as these are helpful features less influenced by sequence length variation than are counts.
dfSeq$Aprop <- (dfSeq$A) / (dfSeq$A + dfSeq$T + dfSeq$C + dfSeq$G)

dfSeq$Tprop <- (dfSeq$T) / (dfSeq$A + dfSeq$T + dfSeq$C + dfSeq$G)

dfSeq$Gprop <- (dfSeq$G) / (dfSeq$A + dfSeq$T + dfSeq$C + dfSeq$G)

#After each batch of additions, check out dfSeq in the viewer. Look at the new columns appended on the end, and check to ensure that is what we want. Note that above we were doing basic arithmetic (vectorized fashion R style!).

#Adding dinucleotide frequency (k-mers of length 2), here using proportions as that let's us account for variability in sequence lengths.
dfSeq <- cbind(dfSeq, as.data.frame(dinucleotideFrequency(dfSeq$nucleotides2, as.prob = TRUE)))

#Let's check the viewer again.

#Adding trinucleotide frequency (k-mers of length 3)
dfSeq <- cbind(dfSeq, as.data.frame(trinucleotideFrequency(dfSeq$nucleotides2, as.prob = TRUE)))

#Note how the df size grows (number of variables), as we increase k. Have a look at dfSeq in viewer.

#Here is code to add k-mers of length 4 (as probabilities/proportions). For now, we won't run this as we will consider the simpler features first. But I'm sharing some example code here for your future reference, and note that we can continue to increase k further, such as if needed in the future to help with more difficult classification problems. Run this on your own and check it out.

#dfSeq <- cbind(dfSeq, as.data.frame(oligonucleotideFrequency(x = dfSeq$nucleotides2, width = 4, as.prob = TRUE)))

####5 - TRAINING CLASSIFICATION MODEL #1: MARKERS----

#First, we are going to give the random forest algorithm a relatively "easy" classification problem. We will ask whether the COI gene and 28S genes can be accurately classified from simple sequence features. As we saw in an earlier script, there was a large difference in nucleotide composition between these two genes in Cotesia, and therefore we SHOULD be able to classify sequences to one of these two genes on the basis of their sequence features. But let's see!

#We will confine this test to the Cotesia genus, which has both genes of interest (COI and 28S) represented.
dfCotesia <- dfSeq[dfSeq$genus_name == "Cotesia", ]
unique(dfCotesia$genus_name)

#converting string format back to character data, so we can more readily use tidyverse functions.
dfCotesia$nucleotides2 <- as.character(dfCotesia$nucleotides2)

#Seeing counts by markercode
table(dfCotesia$markercode)

#As our maximum sample size for these two genes is 82 for the 28S gene in dfCotesia, we will below sample 20 individuals (about 25% of the total for the 28S gene) from each gene to serve as the validation data set. These samples will be ENTIRELY set aside during the model training process. These samples should NOT be used until validating the model later. sample_n() is a very useful function for sampling a specific number of rows from each group. Here, we are sampling 20 rows from each group (i.e. each marker code). See the resulting object to confirm.

#setting seed, as we are using randomization in next step, and we want results to be reproducible.

set.seed(217)

dfValidation <- dfCotesia %>%
  group_by(markercode) %>%
  sample_n(20)

#Checking object. Make prediction about what we should see. checking sample size for each marker is the same.
table(dfValidation$markercode)

#Now, we are creating a training dataset that does NOT overlap with the validation set. To do this, we will first remove processids that were among the samples randomly selected for the validation dataset. We can do this by asking for processid's that are not (!) in the validation set using %in%. Second, among records remaining that are not in the validation set, we will pick 62 individuals for each gene to serve in the training set. Note that this is a small sample size; this script is an EXAMPLE. The sample sizes don't necessarily need to be exactly equal, but in general we need to aim for "class balance" (i.e. similar representation between classes, or groups). If there is a large imbalance in sample size, as would be the case here if we used all COI sequences, the classifier could achieve good performance simply by saying that everything is COI! That is NOT what we want... we want to train the classifier to recognize COI vs. 28S sequences.
set.seed(13)
dfTraining <- dfCotesia %>%
  filter(!processid %in% dfValidation$processid) %>%
  group_by(markercode) %>%
  sample_n(62)

#Checking we should have 62 records of each gene
table(dfTraining$markercode)

#Checking our variable names and the column numbers.
names(dfTraining)

#Next, we are building a classifier to separate the COI and 28S genes in these datasets, using first the A, T, and G proportion (columns 86-88) as predictors. Then, if needed, we can see if it is helpful to add more complex features (e.g. dinucleotide frequencies, 3-mers, etc.). The response variable is markercode; we are trying to predict which gene a sequence belongs to, on the basis of simple sequence features alone. Note I set ntree fairly small as this should be an easy classification problem and we have a small dataset (we want every row to be predicted several times).
gene_classifier <- randomForest::randomForest(x = dfTraining[, 86:88], y = as.factor(dfTraining$markercode), ntree = 50, importance = TRUE)

#Let's look at results.
gene_classifier
#Excellent performance!

#We can also dive into the specifics, such as by looking at the following. If we type the name of our classifier and $, then we can use arrows to scroll through the available options. See the documentation for the randomForest() function for more detailed information.

#relative importance of each feature
gene_classifier$importance

#number of times each row in the data set was left "out of bag" and predicted
gene_classifier$oob.times

#out of bag error rate based upon all prior trees built up to that iteration. (i.e. error rate for observations left out of that iteration and predicted based upon all decision trees built up to that point). Here, we can see we had no errors from the start!
gene_classifier$err.rate

#confusion matrix. (true positives, true negatives, false positives, and false negatives)
gene_classifier$confusion

#for each row, fraction of votes from the random forest for each class. In this example, there is agreement. (Contrast this with the outcome for a more difficult classification problem.)
gene_classifier$votes

#This worked so well, let's actually try with just 2 features
gene_classifier2 <- randomForest(x = dfTraining[, 86:87], y = as.factor(dfTraining$markercode), ntree = 50, importance = TRUE)
gene_classifier2
#still perfect performance.

#So, this was a very good classifier! Basically, these genes are very different in their nucleotide frequency, and this was easy for the classifier to pick up.

#But, will this work as a general tool, on sequences unseen during training the classifier? Here, we are using our first fitted classifier to predict the marker for our unseen data (set aside above as the validation data set) based upon nucleotide proportions.
predictValidation <- predict(gene_classifier, dfValidation[, c(70, 86:88)])

#Let's have a look at the data format here. There are the predictions.
predictValidation
class(predictValidation)
length(predictValidation)

#Create our own confusion matrix.
table(observed = dfValidation$markercode, predicted = predictValidation)

#Yes, this classifier also works on the unseen data.

#NOTE: This was a simple EXAMPLE, where we are only trying to distinguish TWO genes, which have very different sequence properties, as we saw in our earlier Biostrings script. As a next step, let's try to solve a more challenging classification problem! Then, move on to defining your own problem and try to build a classifier to solve it.

#Removing objects not needed further (note one of the objects below is generated by the example answers for the activity sheet, and you might not have that object name - that's OK ... you'll just get a notice that that wasn't found, not an error.)
rm(dfCOI28S, dfCotesia, dfSeq, dfTraining, dfValidation, gene_classifier, gene_classifier2, predictValidation, countsbygene)

####6 - TRAINING CLASSIFICATION MODEL #2: TAXONOMY----

#CHALLENGE. Using COI, and following the model in part 5, build a genus-level classifier for the two genera in the dfBOLD data set.

####7 - CODING CHALLENGE----

#RECOMMENDATION: Work in a study group for this challenge.

#CODING CHALLENGE: Build a similar random forest classifier but that tackles a more difficult problem. For example, can you build a classifier that separates COI and cytb sequences from the same genus? Those are both mitochondrial, protein-coding genes from the same taxonomic group, so may have similar features. So, that would likely be a more difficult classification problem. Or, choose a different example derived from your interests.

#TIP: You can obtain your data sets from public databases. Also, remember that R has a wide variety of built-in data sets that you can also explore.

#Type this to see list of available datasets
data()

#Have fun, and share your classifier with the class! You can post about what you did to the Discussion Board on CourseLink.
