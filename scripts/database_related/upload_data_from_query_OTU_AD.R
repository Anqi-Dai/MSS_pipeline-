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

library(seqinr);
library(stringr); #For string replace.
#library(stringi); #For string replace.
#library(caroline); #for `dbWriteTable2`; #https://stackoverflow.com/questions/30276108/how-to-use-dbwritetable2-properly
library("biomformat");
#library("xlsx"); #For read.xlsx function;

if(file.exists('~/projects/general/library/antoniostats/intervalcluster2.R')){
	source('~/projects/general/library/antoniostats/intervalcluster2.R')
}else{
  print("Intervalcluster2.R not loaded! Check path if you need it!")
}

source('~/MSK/work/microbiome_db/SQL/scripts/get_data_from_query_OTU.R')
source('~/MSK/work/microbiome_db/SQL/scripts/db_connect.R'); #Start a connection here; The connection variable is `con`.

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
    
    d_set_input = fread("/Volumes/vandenbrinklab/Angel_Dai/Full_human_shotgun_catalog/full_human_shotgun_catalog_updated.csv")
    # new samples from Italy and upenn
    d_set_input = fread("/Volumes/vandenbrinklab/Angel_Dai/Full_human_shotgun_catalog/shotgun_cart_earlier.csv")
    
    d_set_input = fread("/Volumes/vandenbrinklab/Angel_Dai/Full_human_shotgun_catalog/shotgun_update_1204.csv")
    d_set_input = fread("/Volumes/vandenbrinklab/Angel_Dai/Full_human_shotgun_catalog/shotgun_update_1207.csv")
    d_set_input = fread("/Volumes/vandenbrinklab/Angel_Dai/Full_human_shotgun_catalog/shotgun_update_1208.csv")
    
    d_set=data.frame(directory=d_set_input$directory,
                     projectid=d_set_input$projectID,
                     sampleid=d_set_input$sampleid,
                     fid=d_set_input$fid);
    
    update_data_from_query_OTU_check_and_submission(table_name, d_set);
  }
  
}


