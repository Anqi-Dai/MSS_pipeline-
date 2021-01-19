# this script will take 
# a counts table in the format of like "humann2_genefamilies_uniref90_go_cpm_unstratified_reordered.tsv",
# And a pheno table 
# And remove the rows that are zero in the counts table and do a Bray Curtis distance PCoA
# And return a ggplot object

library(tidyverse)
library(vegan)

plot_BC_PCoA <- function(cts, pheno, Var) {

  # remove the rows that are zero in the counts table
  idx_row_zero <- cts %>% 
    dplyr::rename(fea_name = names(.)[1]) %>% 
    column_to_rownames('fea_name') %>% 
    as.matrix() %>% 
    apply(., 1, sum)
  
  # calculate the bc distance
  CTS <- cts %>% 
    dplyr::rename(fea_name = names(.)[1]) %>% 
    column_to_rownames('fea_name') %>% 
    as.matrix() 
  
  bc_dist <-  vegdist(t(CTS), "bray")
  pc_bc <- cmdscale(bc_dist, k = 2) 

  
  g <- pc_bc %>% 
    as.data.frame() %>% 
    rownames_to_column('sampleid') %>% 
    full_join(pheno, by = 'sampleid') %>% 
    ggscatter(x = 'V1', y = 'V2', color =  Var, fill = Var, palette = 'lancet')+
    labs(x = 'PCo 1', y = 'PCo 2', title = str_glue('PCoA plot for response to {Var}')) 
  
  return(g)
}