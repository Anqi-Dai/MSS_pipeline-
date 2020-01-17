#! /usr/bin/env Rscript
'Write out csv file of per sample RPKM
Usage:
  produce_per_sample_shortbred_result.R <input_path> <output_fn>

Arguments:
  <input_path> a folder path that is the quantification result from shortbred with the path
  <output_fn> a csv file with two columns, the fid and the RPKM in a path
' -> doc

# install if the package is missing
list.of.packages <- c("tidyverse", "docopt")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos = "http://cran.us.r-project.org")

suppressMessages(library(tidyverse))
suppressMessages(library(docopt))

arguments <- docopt(doc)

inpath <- arguments$input_path

outfn <- arguments$output_fn

# the code

fns <- list.files(inpath, full.names = T)

all <- fns %>% 
  set_names(fns) %>% 
  map(~ read_tsv(.)) %>% 
  bind_rows(.id = 'sampleid') %>% 
  mutate(sampleid = str_replace(sampleid, '.+//',''),
         sampleid = str_replace(sampleid, '_short.+$','')) %>% 
  group_by(sampleid) %>% 
  summarise(RPKM = sum(Count))

all %>% 
  write_csv(outfn)