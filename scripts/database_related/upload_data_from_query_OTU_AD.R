#Jun/01/2017
#This script is going to update data to my table (sufix `_ag`) on Xavier's database (`plvglover1`,  `microbiome`). 
#
#To assure proper upload, table is going to perform some steps:
#     c1. matches between fields in input data frame and uploaded table;
#     c2. key value incremental to highest number in current table;
#
#Query numbers:
#1. Upload to samples_ag table.
#     TODO: I HAVE to check where the data comes from.
#     1.1 : samples_castory_ag
#     1.2 : samples_duke_ag
#     1.3 : samples_castory_ag
#2. Upload to otu_ag
#     Take data from 	`total.8.otu-tax.biom`. Includes: uploaded_date and run_date in file.
#3. Upload to counts_ag table.
#     Take data from 	`total.8.otu-tax.biom`. Includes: uploaded_date and run_date in file.

#4. Upload to fmt_ag
#     TODO: I HAVE to check where the data comes from.
#5. Clone patients_ag table.
#     TODO: I HAVE to check where the data comes from.
#
#6. antibiotics data.
#
library(tidyverse)
library(seqinr);
library(stringr); #For string replace.
#library(stringi); #For string replace.
#library(caroline); #for `dbWriteTable2`; #https://stackoverflow.com/questions/30276108/how-to-use-dbwritetable2-properly
library("biomformat");
library(tidyverse)
#library("xlsx"); #For read.xlsx function;

if(file.exists('~/projects/general/library/antoniostats/intervalcluster2.R')){
	source('~/projects/general/library/antoniostats/intervalcluster2.R')
}else{
  print("Intervalcluster2.R not loaded! Check path if you need it!")
}

source('/Users/daia1/pipeline/scripts/database_related/get_data_from_query_OTU.R')
source('/Users/daia1/pipeline/scripts/database_related/db_connect.R'); #Start a connection here; The connection variable is `con`.
getwd()
# temp for when I'm using Marissa's laptop
source('/Users/daia1/Downloads/MSK/MSS_pipeline-/scripts/database_related/get_data_from_query_OTU.R')
# if there is a sourcing problem of the above script then go to that script and run directly , the key is to
# have the get_data_from_query_otu function
source('/Users/daia1/Downloads/MSK/MSS_pipeline-/scripts/database_related/db_connect.R');
#warning("TO DO: ADD a step to get last `id` from database.")

my_dbWriteTable <- function(con, table_name, d_set_to_upload_reordered){
  #I created this function to upload data to table and skip duplicates.
  #I should add functionality to add `NULL` or `NA` values;
  
  
  stop("THIS HAS TO BE UPDATED! I AM SEARCHING FOR A GOOD SOLUTION TO AVOID DUPLICATE CONFLIT. I tried `ON CONFLICT`` option, but it did not work");
  
  q0p3_count_pre = get_data_from_query_OTU(0.3,table_name); #get count before upload;
  
  q0p1 = get_data_from_query_OTU(0.1,table_name);
  if( !all(q0p1$column_name==colnames(d_set_to_upload_reordered) ) ){
      stop();
  }
  
  colnames_str = paste(colnames(d_set_to_upload_reordered),collapse = ",");
  data_type_str = sprintf("(%s)", paste(q0p1$data_type,collapse = ",") );
  data_type_str = gsub(pattern = "integer", "%d",data_type_str);
  data_type_str = gsub(pattern = "text", "'%s'",data_type_str);
  data_type_str = gsub(pattern = "date", "'%s'",data_type_str);
  
  values = do.call("sprintf", c(data_type_str, d_set_to_upload_reordered));
  
  #values = apply( d_set_to_upload_reordered[ , colnames(d_set_to_upload_reordered) ] , 1 , paste , collapse = " , " );
  #values = sprintf("(%s)",values);
                  
  query= sprintf("insert into %s (%s) values %s", table_name, colnames_str,
                 paste(values, collapse = ","))
  
  dbSendQuery(con, query);
  
  q0p3_count_post = get_data_from_query_OTU(0.3, table_name); #get count after upload;
  
  total_insert = q0p3_count_post$count - q0p3_count_pret$count;
  print(sprintf("%d rows were inserted!", total_insert));
  
}

upload_data_from_query_OTU_check_and_submission <- function(table_name, d_set_to_upload ){
  #Check if d_set is in proper format to be uploaded in `table_name`.
  
  #Add uploaded_date to date_frame.
  uploaded_date = format(Sys.time(),"%m-%d-%Y");
  d_set_to_upload$uploaded_date = uploaded_date;
  
  #Assign incremental key values starting on `maximum` value in current table.
  
  q_key_max_cur = get_data_from_query_OTU(0.2,table_name);
  if(is.na(q_key_max_cur$max)){
    q_key_max_cur$max = 0;
  }
  d_set_to_upload$key = ( 1:length(d_set_to_upload$uploaded_date) ) + q_key_max_cur$max;
  
  #Get columns for `column names`
  q_column_names = get_data_from_query_OTU(0.1,table_name);
  
  #stop("Add step to compare columns, re-order data, add key and upload!");
  if( ! all(is.element(q_column_names$column_name, colnames(d_set_to_upload) )) ){
    stop( sprintf("Some fields in table %s are not defined in upload dataset", table_name));
  }
  
  d_set_to_upload_reordered = d_set_to_upload[,q_column_names$column_name]
  
  #Add the step to write the data into table!
  print("data getting uploaded!");
  
  dbWriteTable(con, table_name, value = d_set_to_upload_reordered, append = TRUE, row.names = FALSE);
  #my_dbWriteTable(con, table_name, d_set_to_upload_reordered);
  print( sprintf("upload: done! %d rows", length(d_set_to_upload_reordered$key)));
  
}

update_data_from_query_OTU_check_and_submission <- function(table_name, d_set_to_upload){
  #This version check duplicates in table (according to table keys) and add only new, unique data.
  #It intermediates this step with a tempory table.
  #Some interesting discussion on skip duplicates: https://stackoverflow.com/questions/1009584/how-to-emulate-insert-ignore-and-on-duplicate-key-update-sql-merge-with-po
  #Rationale:
  #     1. Create temp_updating table with same structure as target table
  #     2. Upload data to temporary table
  #     3. Retrieve data in temporary table that is not in target table
  #     4. Remove temporary table.
  #     5. Return only novel data to be uploaded;
  
  table_name <- table_name
  d_set_to_upload <- d_set
  #Clean temp_table
  temp_table="temp_updating";
  query_check_temp = sprintf("SELECT EXISTS (SELECT 1 FROM   information_schema.tables WHERE table_name = '%s')", 
                             temp_table);
  q_check = dbGetQuery(con,query_check_temp);
  if(q_check$exists){
    dbSendQuery(con, sprintf("drop table %s", temp_table));  
  }
  
  #Add uploaded_date to date_frame.
  uploaded_date = format(Sys.time(),"%m-%d-%Y");
  d_set_to_upload$uploaded_date = uploaded_date;
  
  #Assign incremental key values starting on `maximum` value in current table.
  q_key_max_cur = get_data_from_query_OTU(0.2,table_name);
  if(is.na(q_key_max_cur$max)){
    q_key_max_cur$max = 0;
  }
  d_set_to_upload$key = ( 1:length(d_set_to_upload$uploaded_date) ) + q_key_max_cur$max;
  
  #Get columns for `column names`
  q_column_names = get_data_from_query_OTU(0.1,table_name);
  #q_column_unique_names = get_data_from_query_OTU(0.1,sprintf("%s_unique",table_name));
  
  #stop("Add step to compare columns, re-order data, add key and upload!");
  if( ! all(is.element(q_column_names$column_name, colnames(d_set_to_upload) )) ){
    stop( sprintf("Some fields in table %s are not defined in upload dataset", table_name));
  }
  
  d_set_to_upload_reordered = d_set_to_upload[,q_column_names$column_name]
  
  print("    creating temporary table to check for duplicates");
  query_create_temp_table = sprintf("CREATE TABLE %s AS SELECT * FROM %s WHERE 1=2",temp_table, table_name);
  dbSendQuery(con,query_create_temp_table);
  
  dbWriteTable(con, temp_table, value = d_set_to_upload_reordered, append = TRUE, row.names = FALSE);
  #my_dbWriteTable(con, table_name, d_set_to_upload_reordered);
  
  ###getting data that does not violate unique constrain in target table.
  query_select_unique_key = paste("select ccu.column_name, ccu.constraint_name from information_schema.table_constraints as tc",
                                  " natural join information_schema.constraint_column_usage as ccu ",
                                  #" on tc.table_name=ccu.table_name and tc.constraint_name=ccu.constraint_name ",
                                  sprintf(" where tc.table_name='%s' and tc.constraint_type='UNIQUE'",table_name),
                                  sep = "");
  q_unique_key = dbGetQuery(con, query_select_unique_key);
  if(length(unique(q_unique_key$constraint_name))>1){
    stop("This table seems to have more than one rule for unique constraint. This script is not ready to work with it");
  }
  unique_constraint_str = paste(q_unique_key$column_name,collapse = ",");
#  query_select_novel_data = sprintf("select * from %s where %s as not in (select %s from %s)",
#                                    temp_table,
#                                    unique_constraint_str,
#                                    unique_constraint_str,
#                                    table_name);
  unique_constraint_join_str = paste(sprintf("t1.%s=t2.%s",q_unique_key$column_name,q_unique_key$column_name), collapse = " AND ")
  
  query_select_novel_data = paste(sprintf("select t1.* from %s as t1 left join %s as t2", temp_table, table_name),
                                  sprintf(" on %s ", unique_constraint_join_str),
                                  sprintf("where t2.%s is NULL", q_unique_key$column_name[1]),
                                  collapse = " ");
  
  d_set_unique = dbGetQuery(con,query_select_novel_data);
  
  print(sprintf("%d/%d novel rows to be updated!", 
                dim(d_set_unique)[1],
                dim(d_set_to_upload_reordered)[1]) );
  if(length(d_set_unique)==0){
    print("      No updates! No novel data.");
    return();
  }
  upload_data_from_query_OTU_check_and_submission(table_name,d_set_unique);
}


upload_data_from_query_OTU <- function(query_number, ...){
  
  args = list(...);
  if(length(args)){
    input_data_file = args[[1]];
  }
  if(query_number==1){
    table_name = "shotgun_lookup_ad";
    
    d_set_input = fread("~/Downloads/MSK/Catalog/data/new_20211014.csv")
    
    d_set=data.frame(directory=d_set_input$directory,
                     projectid=d_set_input$projectid,
                     sampleid=d_set_input$sampleid,
                     fid=d_set_input$fid);
    
    update_data_from_query_OTU_check_and_submission(table_name, d_set);
  }
  
  if(query_number==2){
    table_name = "picrust2_pathway_counts";
    
    d_set_input = read_csv('~/pipeline/scripts/picrust2/data/normalized_picrust2_pred_pathway_abundance_all.csv')

    
    d_set=data.frame(pwid=d_set_input$PWID,
                     sampleid=d_set_input$sampleid,
                     cpm=d_set_input$cpm);
    
    update_data_from_query_OTU_check_and_submission(table_name, d_set);
  }
  
  if(query_number==3){
    table_name = "metacyc_pathway_name";
    
    d_set_input = read_tsv('~/pipeline/scripts/shotgun_pipeline/data/metacyc_pathway_name_and_ID.tsv')
    
    
    d_set=data.frame(pwid=d_set_input$PWID,
                     pw_name=d_set_input$pw_name);
    
    update_data_from_query_OTU_check_and_submission(table_name, d_set);
  }
  
  if(query_number==4){
    table_name = "metacyc_pathway_ontology";
    
    d_set_input = read_csv('~/pipeline/scripts/shotgun_pipeline/data/metacyc_pathway_class_and_superclass_levels.csv', col_types = 'ccccccccccccccccc')
    
    
    d_set=data.frame(pwid=d_set_input$pwid,
                     l1=d_set_input$L1,
                     l2=d_set_input$L2,
                     l3=d_set_input$L3,
                     l4=d_set_input$L4,
                     l5=d_set_input$L5,
                     l6=d_set_input$L6,
                     l7=d_set_input$L7,
                     l8=d_set_input$L8,
                     l9=d_set_input$L9,
                     l10=d_set_input$L10,
                     l11=d_set_input$L11,
                     l12=d_set_input$L12,
                     l13=d_set_input$L13,
                     l14=d_set_input$L14,
                     l15=d_set_input$L15,
                     l16=d_set_input$L16);
    
    update_data_from_query_OTU_check_and_submission(table_name, d_set);
  }
  
  if(query_number==5){
    table_name = "qpcr_16s_ag";
    
    d_set_input = fread("~/Downloads/MSK/Catalog/data/qpcr_20211020.csv")
    
    d_set=data.frame(sample_id=d_set_input$sample_id,
                     copy_number_16s=d_set_input$copy_number_16s,
                     copies_16s_per_g=d_set_input$copies_16s_per_g,
                     comments=d_set_input$comments,
                     sample_id_unique = d_set_input$sample_id_unique);
    
    update_data_from_query_OTU_check_and_submission(table_name, d_set);
  }
  
  if(query_number==6){
    table_name = "samples_castori_ag";
    input_data_file = '~/Downloads/MSK/MSS_pipeline-/scripts/dada2_pipeline/data/tblSamples.csv'
    Samples_castori_center_file = input_data_file;
    
    castori_downloaded_date = '2021-10-07'
    #castori_data = read.table(Samples_castori_center_file,sep=",", quote = "", comment.char = "", header = T);
    castori_data = read_csv(Samples_castori_center_file)
    mrn_integer = as.character(castori_data$MRN);
    mrn_integer = suppressWarnings( as.numeric(mrn_integer) );
    mrn_integer[is.na(mrn_integer)]=-1;
    mrn_string = as.character(castori_data$MRN);
    mrn_string[mrn_integer!=-1]=NA;
    uploaded_date = format(Sys.time(),"%m-%d-%Y");
    d_set = data.frame(#id = 1:length(castori_data$Sample_ID),
      sampleid = castori_data$Sample_ID,
      mrn = mrn_integer,
      mrn_str = mrn_string,
      datecollection = lubridate::mdy_hms(castori_data$DateCollection),
      sampletype = castori_data$SampleType,
      consistency = castori_data$Consistency,
      datereceived = as.Date(castori_data$DateReceived,"%d-%b-%y"),
      datealiquot = as.Date(castori_data$DateAliquot,"%d-%b-%y"),
      boxrawstool = castori_data$BoxRawStool,
      numberaliquots = castori_data$NumberAliquots,
      numberstartingaliquots = castori_data$NumberStartingAliquots,
      comment = castori_data$Comment,
      timecollected = castori_data$TimeCollected,
      dropofftime = castori_data$DropOffTime,
      castori_downloaded_date =  lubridate::ymd(castori_downloaded_date),
      #uploaded_date = uploaded_date,
      stringsAsFactors = F);
    
    no_date_str = as.Date("5/5/5555", "%m/%d/%Y");
    col_names = colnames(d_set);
    is_date_type = is.element(col_names, c( "DateCollection", "DateReceived", "DateAliquot", "castori_downloaded_date"));
    ind_date_type = which(is_date_type);
    d_set[,is_date_type][is.na(d_set[,is_date_type])] = no_date_str;
    d_set[,ind_date_type]=format(d_set[,ind_date_type],"%m-%d-%Y")
    #d_set[,ind_date_type]= as.character(format(d_set[,ind_date_type],"%m-%d-%Y"));
    
    
    update_data_from_query_OTU_check_and_submission(table_name, d_set); 
  }
  
}


