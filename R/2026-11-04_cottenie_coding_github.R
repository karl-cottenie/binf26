##***************************
## example of giusthub code
##
## Karl Cottenie
##
## 2026-11-04
##
##***************************

## _ Packages used -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(viridis)
library(Biostrings)
# Install the core sequence container and a built-in model genome (E. coli)
# BiocManager::install(c("BSgenome.Ecoli.NCBI.20080805"))
library(BSgenome.Ecoli.NCBI.20080805)
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

## example of a fasta file -----

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
