# this script takes two input table: the original split unstratified feature table and the finalzed pheno table and 
# return two output tables that could later be output to files:
# 1. The reordered feature table
# 2. The horizontal pheno table.

library(tidyverse)

prepare_lefse_input_tables <- function(pheno, cts) {
  CTS = cts %>% 
    rename_all(funs(str_replace(., '_.+$',''))) 
  
  CTS <- bind_cols(CTS[,1],CTS[,pheno$Sampleid]) 
  return(cts=CTS)
}