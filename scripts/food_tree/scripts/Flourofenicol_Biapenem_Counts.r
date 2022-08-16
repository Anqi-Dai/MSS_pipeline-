rm(list=ls())


library(googlesheets4)
# https://twdtradewinds.com/login/
library(ggplot2)
library(lubridate)
library(data.table)
library(ggplot2);
library(dplyr);
library(reshape2);
library(tidyr);
library(forcats)
library(stringr)
library(Hmisc)
library(gtsummary)
library(flextable)
library(fuzzyjoin)
library(stringi)
library(scales)
library(ggrepel)
library(tidyverse)

# Load in the CFU count data from google sheets
#alldat<- read_sheet("https://docs.google.com/spreadsheets/d/1Nt7jVu_NHtfk0ukZ7mRcMzQMeWlvXdU26f3BvvnWleA/edit?usp=sharing")
alldat<- read.csv("~/Desktop/Antibiotics/Exp9_Flourfenicol_Biapenem_Diets.csv")

### OR load from Excel on your desktop
# alldat <- read.csv("/Users/peledj/Downloads/VRE experimental design - AllCountData.csv")
twdkey <- readxl::read_excel("/Users/rangesam/Downloads/MemorialSloan-Fisher-T334-6SPR-PO#C21724859-06.02.21.xls")
colnames(twdkey) <- c("TubeID", "tubewt")
experiment <- 10
exp_title <- "Biapenem + Diet"


alldat1 <- alldat %>% 
  left_join(twdkey, by='TubeID') %>% 
  mutate(sweight = Total_weight_stool_and_tube_g - tubewt) %>%
  mutate(Log_CFUs_per_GramStool = (((((CFUs*(100/PlatedmicroL))/100)/(10^Dilution) )*1000)/sweight)+1 ) %>%
  filter(ProblemData==0) %>% 
  mutate(Day=factor(Day))
  #filter(Experiment==experiment) 


ggplot(alldat1, aes(x=Day, fill=Day, y=as.numeric(Log_CFUs_per_GramStool))) + 
  geom_dotplot(dotsize=3, binaxis='y', stackdir='center', position=position_dodge(0.5), binwidth = 0.05) +
  ylab("CFUs/gram of stool, Log10 Scale") +
  scale_y_log10(breaks = c(1 %o% 10^(-1:12)), 
                labels = trans_format("log10", math_format(10^.x))) +
  facet_grid(.~Treatment)+ ggtitle(exp_title) +

  theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust = 0.5)) +
  geom_text_repel(data=alldat1, aes(Day,as.numeric(Log_CFUs_per_GramStool), label=Mouse_identifier), fontface='bold')+
  theme(legend.position = "none")

