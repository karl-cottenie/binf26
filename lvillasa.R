# lvillasa.R

## _ Packages used -------
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(viridis)
# + scale_color/fill_viridis(discrete = T/F)
theme_set(theme_light())
library(randomForest)
library(Biostrings)
# _ Code your own solution -----

# CREATE A NEW BRANCH - do this after every pull request has been merged!

usethis::pr_init(branch = "formidable") # change formidable w/ more informative branch name!

usethis::git_sitrep()
# check that Git local (project) Current local branch is formidable
# and this is differet from the git_sitrep output from line38!

# now do your R work
# for this assignment, create a script file named
# <your-university-user-name-from-your-uog-email-address>.R
# for me, that would be cottenie.R
# this file will consist of your commented script file for the randomforest code
# put it in the correct folder!

# COMMIT your new file with the changes in the "git" tab next to Environmnet

# _ Submit pull request -----

usethis::pr_push(),

# in the launched browser window, click "Create pull request"
# in the description, add - "This is my random forest script file :-)"
# click "Create pull request"

# then you wait for me to review and accept - this is your "quiz" for this week

# after I have merged the pull request

usethis::pr_pull() # make sure your branch is up-to-date

usethis::pr_finish() # cleaning everything up, deleting local 

usethis::git_sitrep()




# 1. Access the built-in Escherichia coli genome database object
genome <- BSgenome.Ecoli.NCBI.20080805
print(genome)

# 2. Extract a built-in sequence (e.g., Chromosome NC_010473)
raw_sequence <- genome[["NC_010473"]]

# 3. Convert it to a DNAStringSet (The standard R FASTA database structure)
# We can slice it (e.g., first 500 bases) to keep the example concise
fasta_db <- DNAStringSet(Views(raw_sequence, start = 1, end = 500))
fasta_db

# 4. Assign a FASTA header name to the sequence definition
names(fasta_db) <- "NC_010473.1 Escherichia coli str. K-12 substr. DH10B, complete sequence"

# 5. 3-mer 
oligonucleotideFrequency(fasta_db, width = 3)

# 6. 4-mer
oligonucleotideFrequency(fasta_db, width = 4)
