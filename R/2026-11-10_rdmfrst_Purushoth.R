#RANDOMFOREST MODEL-------
rm(list = ls())
##_ Packages used -------
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(viridis)
# + scale_color/fill_viridis(discrete = T/F)
theme_set(theme_light())
library(randomForest)
library(Biostrings)

##LOAD DATA -----
dfBOLD <- read_csv("../data/2026-10-26_t_spiralis.csv")

ls_non_te = c("simple_repeat", "Simple_repeat", "Low_complexity", "Satellite", "rRNA", "snRNA", "ARTEFACT", "Unknown")

ls_retro_te <- c("LINE/CR1", "LINE/RTE-BovB", "SINE/tRNA-RTE", "SINE/tRNA-Deu-L2", "LTR/Pao", "LTR/Gypsy", "LTR/Ngaro", "LTR/DIRS", "PLE", "PLE/Naiad")

ls_te <- c("DNA/TcMar-Tc1", "DNA/TcMar-Tc4", "DNA/TcMar-Mariner", "DNA/TcMar-Pogo", "DNA/hAT", "DNA/hAT-Ac", "DNA/PIF-Harbinger", "DNA/MULE-MuDR", "DNA/CMC-Chapaev", "RC/Helitron")

#k-mer
dfBOLD$nucl_seq <- as.character(dfBOLD$nucl_seq)

dfBOLD$A <- str_count(dfBOLD$nucl_seq, "A")
dfBOLD$T <- str_count(dfBOLD$nucl_seq, "T")
dfBOLD$C <- str_count(dfBOLD$nucl_seq, "C")
dfBOLD$G <- str_count(dfBOLD$nucl_seq, "G")

total_bases <- dfBOLD$A + dfBOLD$T + dfBOLD$C + dfBOLD$G
dfBOLD$Aprop <- dfBOLD$A / total_bases
dfBOLD$Tprop <- dfBOLD$T / total_bases
dfBOLD$Cprop <- dfBOLD$C / total_bases
dfBOLD$Gprop <- dfBOLD$G / total_bases

set.seed(217)

top_classes <- names(sort(table(dfBOLD$matching_class), decreasing = TRUE)) [1:2]
top_classes

df_task5 <- dfBOLD %>%
  filter(matching_class %in% top_classes)

rf_task5 <- randomForest(
  x = df_task5[, c("Aprop", "Tprop", "Cprop", "Gprop")], 
  y = as.factor(df_task5$matching_class), 
  ntree = 50,
  importance = TRUE)

rf_task5
saveRDS(rf_task5, file = "../output/rf_task5_model.rds")