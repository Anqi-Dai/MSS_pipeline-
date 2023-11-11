
# This script will format a beautiful author list for manuscript submission

# install.packages("officer")
# browseVignettes(package = "officer")    # https://ardata-fr.github.io/officeverse/

library(tidyverse)
library(readxl)
library(officer)


########### If needed, use this to generate a template file to populate manually
########### in Excel and can save as an xls file to read back in

# template <- data.frame(Name_with_middle_initial = c("Jonathan U. Peled"),
#            Order.listed = c(1),
#            Equal.contribution = c(TRUE), 
#            Degree.with.periods. = c("M.D., Ph.D."),
#            email = c("peledj@mskcc.org"),
#            Affiliation1 = c("Adult Bone Marrow Transplantation Service, Department of Medicine, Memorial Sloan Kettering Cancer Center, New York, NY, USA"),
#            Affiliation2 = c("Department of Immunology, Sloan Kettering Institute, Memorial Sloan Kettering Cancer Center, New York, NY, USA"),
#            Affiliation3 = c(NA))
##write.csv(template, "author_list_template.csv",row.names=FALSE)


########## Pull in Excel file of author affiliations

# author.filename <- "/Volumes/vandenbrinklab/Susan/Human GVHD/authorship_COI/Author_list_TCR_manuscript_fromTsoni.xlsx"
author.filename <- "/Users/peledj/Desktop/desktop_stuff/Susan Dewolf TCR seqencing/authorship_COI/Author_list_TCR_manuscript_fromTsoni.xlsx"
author.filename <- "/Users/peledj/Desktop/desktop_stuff/Susan Dewolf TCR seqencing/authorship_COI/Author_list_TCR_2022-01-27.xlsx"
df <- read_excel(author.filename)


df <- df[,colSums(is.na(df))<nrow(df)] # remove na columns
df <- df[rowSums(is.na(df))<ncol(df),]  # remove na rows

# add placeholder dummies so that each author has at least one primary affiliation
df$Affiliation1[is.na(df$Affiliation1)] <- "placeholder.institution"

# melt so that a unique ID can be assigned to each affiliation
df.m <- df %>% 
  arrange(Order.listed) %>% 
  select(-Equal.contribution, -Order.listed, -Degree.with.periods., -email) %>%  #  remove extra columns
  reshape2::melt(id.vars = c("Name_with_middle_initial"),
                       variable.name = "affiliation.ordinal",
                       value.name = "affiliation")

# pull in author order
df.m <-df.m %>% 
  left_join(
    select(df, one_of(c("Name_with_middle_initial", "Order.listed"))),
    by = "Name_with_middle_initial") 

# Define order of affiliations according to author order
affiliations.list <- df.m %>% 
  filter(!is.na(affiliation)) %>% 
  arrange(Order.listed) %>% 
  filter(!duplicated(affiliation)) %>% 
  select(one_of(c("affiliation", "Order.listed")))

affiliations.list$affiliation.id <- seq(1:nrow(affiliations.list))

# Push affiliation ids back into the melted author list
df.m <- left_join(df.m, 
                  select(affiliations.list, -Order.listed),
                  by = "affiliation")

df.m <- df.m %>% filter(!is.na(affiliation)) #remove rows from melted list that have no secondary, tertiary affiliation

df.m <- df.m %>% 
  group_by(Name_with_middle_initial) %>% 
  mutate(superscripts = paste(affiliation.id[1], 
                               affiliation.id[2][!is.na(affiliation.id[2])], 
                               affiliation.id[3][!is.na(affiliation.id[3])],sep = ","))

# Remove extra commas
df.m$superscripts <- gsub(",$", "", df.m$superscripts)
df.m$superscripts <- gsub(",$", "", df.m$superscripts)



author.superscripts <- df.m %>% 
  filter(!duplicated(Name_with_middle_initial)) %>% 
  select(Name_with_middle_initial, superscripts)

df <- left_join(
  df,author.superscripts, by = "Name_with_middle_initial")


#OPTIONAL: Poor Man's Output to Excel
setwd("~/Desktop/desktop_stuff/"); getwd()
# List of authors
# paste0(df$Name_with_middle_initial, ", ", df$Degree.with.periods., df$superscripts) %>% 
  # write.csv("author.list_2019-10-16.csv")

# List of institutions
# paste0(affiliations.list$affiliation.id, affiliations.list$affiliation) %>% 
  # write.csv("institution.list_2019-10-16.csv")


############ FANCIER OPTION: generate word document

# Need to generate, on your own, a Word document into which R will write the author list
# So make a blank word document and save it the appropriate folder, 
# call it "author_list_1.docx"
setwd("/Users/peledj/Desktop/desktop_stuff/Susan Dewolf TCR seqencing/authorship_COI")
getwd()

#define the superscript format
sup.format <- fp_text_lite(vertical.align = "superscript") 

# If including degrees, need to add a comma & a space before each degree
df$Degree.preceeded_by_comma_space <- paste(", ", df$Degree.with.periods.)

# generate lists of ftext objects for each name chunk
names_ftext <- df$Name_with_middle_initial %>%
  map(ftext)

degree_ftext <- df$Degree.preceeded_by_comma_space %>%
  map(ftext)

affiliation_ftext <- df$superscripts %>%
  map(ftext, sup.format)

# concatenate name ftext objects in order with ","
fpar_args <- list()
for (i in seq_along(names_ftext)) {
  print(i)
  fpar_args <- c(
    fpar_args,
    list(names_ftext[[i]]),
    # list(degree_ftext[[i]]),       ####Toggle on and off for degrees
    list(affiliation_ftext[[i]])
  )
  if (i < length(names_ftext)) {
    fpar_args <- c(fpar_args, ", ")
  }
}

# calls fpar using fpar_args as arguments
formatted_names <- do.call(fpar, fpar_args)

# Generate the author list paragraph chunk
read_docx() %>%
  body_add_fpar(formatted_names) %>%
  print(target = "author_list_2.docx")


# Generate list of institutions: 
# This is output as a csv that you copy into Word & fix the superscripts in the affiliation list manually
paste0(affiliations.list$affiliation.id, affiliations.list$affiliation) %>%
write.csv("institution.list_2022-01-27.csv", row.names=FALSE)
