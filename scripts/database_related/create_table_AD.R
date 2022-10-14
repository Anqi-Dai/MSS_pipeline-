#Jun/20/2017
#This script pre-defined a set of table. When possible, tables are named according to Xavier's database (`plvglover1`,  `microbiome`). 
#The new table initiates empty.
#Every table are created with field: `key` for primary key and `uploaded_data` to keep track when data was uploaded.
#otu_ag table has the additional field `running_day` to keep track of when dataset was created.
#
#The default annotation for fields that refer to foreign keys is the `FOREIGNTABLENAME` + `_key` `FOREIGNTABLENAME_key`.

library(RPostgres)
pFile <- read.table("~/dbConfig.txt", header = TRUE, sep = ",", colClasses = "character")
#drv <- dbDriver("PostgreSQL")
con <- dbConnect(RPostgres::Postgres(),
                 dbname = "microbiome",
                 host = "plvglover1.mskcc.org", port = 5432,
                 user = pFile$user, password = pFile$pass)



source('db_connect.R'); #Initialize connection.


check_column_exists <- function(connection, table, column) {
  
  
  query = paste( "SELECT column_name FROM information_schema.columns",
                 sprintf("WHERE table_name='%s' and column_name='%s'",
                         table,
                         column),
                 sep=" ");
  
  r = dbGetQuery(con, query);
  
  col_exists = length(r)>0;
  return(col_exists);
  
}

add_field <- function(connection, target_table, new_field, field_type) {
  #reference_table : table to be cloned.
  
  query = paste( sprintf("alter table %s ", target_table),
                 sprintf("ADD column %s %s", new_field, field_type),
                 sep=" ");
  
  col_exists = check_column_exists(con, target_table, new_field);
  
  if(col_exists){
    warning( sprintf("table %s already exists", target_table));
  }else{
    dbSendQuery(con, query);
  }
  
}

grant_access <- function(connection, grant_type="restricted", table_name){
  
  query_grant_read_access = sprintf("GRANT SELECT ON table %s TO littmane, taury, visitor, xavierj, gomesa, mooreg, peledj, guest",
                                    table_name);
  #When creating a new user, I need to grant `usage` to schema (see this: https://dba.stackexchange.com/questions/98892/minimal-grants-for-readonly-single-table-access-on-postgresql)
  #e.g. grat usage on schema public to GUEST
  
  dbSendQuery(con, query_grant_read_access);
  
  if(grant_type=="opened"){
    #query_grant_access = "GRANT SELECT, INSERT, UPDATE, DELETE ON ALL TABLES IN SCHEMA public TO littmane, taury, visitor, xavierj, gomesa";
    query_grant_access = sprintf("GRANT SELECT, INSERT, UPDATE, DELETE ON table %s TO littmane, taury, visitor, xavierj, gomesa",
                                 table_name);
  }
  
  if(grant_type=="restricted"){
    query_grant_access = sprintf("GRANT SELECT, INSERT, UPDATE, DELETE ON table %s TO littmane, taury, xavierj, gomesa",
                                 table_name);
    
    
  }
  
  dbSendQuery(con, query_grant_access);
}

add_unique <- function(connection, table_name, unique_set=NULL, not_null_constraint=T){
  
  if(is.null(unique_set)){
    return();
  }
  unique_str = sprintf("(%s)",
                       paste(unique_set, collapse = ","));
  
  query_unique = sprintf("ALTER TABLE %s ADD CONSTRAINT %s_unique UNIQUE %s;",
                         table_name,
                         table_name,
                         unique_str);
  
  dbSendQuery(con, query_unique);
  
  if(not_null_constraint){
    for(us in unique_set){
      not_null_str = sprintf( "alter table %s alter %s set not null",
                              table_name,
                              us);
      dbSendQuery(con, not_null_str);
    }
  }

}

create_table <- function(connection, table_name, table_fields, field_type, unique_set = NULL, 
                         add_key_field=T, create_uploaded_date_field=T, access_type="restricted") {
  #reference_table : table to be cloned.
  
  #TO DO: Create primary key: `ALTER TABLE samples_temp_ag ADD PRIMARY KEY (id);`
  already_exists = dbExistsTable(con, table_name); 
  
  if(add_key_field){
    #This is how `primary key` is created: alter table patients_ag add column key serial primary key ;
    if( all(table_fields!="key") ){
      table_fields = c("key", table_fields);
      field_type = c("serial primary key", field_type);
    }
  }
  
  table_and_field_str = paste(sprintf("%s %s", table_fields, field_type), collapse = " , ");
  query = paste( sprintf("create table %s ", table_name),
                 sprintf("( %s )", table_and_field_str),
                 sep=" ");
  
  
  if(already_exists){
    warning( sprintf("table %s already exists", table_name));
  }
  else{
    dbSendQuery(connection, query);
  }
  
  if(create_uploaded_date_field){
    new_field = "uploaded_date";
    field_type = "date"
    add_field(connection, table_name, new_field, field_type);
  }
  
  grant_access(connection, access_type, table_name);
  
  #if(!is.null(unique_set)){
  add_unique(connection, table_name, unique_set);    
  #}

  # 1
  return(table_name);
}


create_table_type <- function(query_type, ...){
  #1. Samples: I believe it is better to have three distinct tables and join them with proper id. 
  #I should add a field `source` on counts table to describe where the sample_id cames from. 
  #1.1 samples_castori_ag;
  #1.2 samples_duke_ag;
  #1.3 samples_regensburg_ag;
  
  if(query_type==1){
    table_name = "shotgun_lookup_ad";
    
    table_fields=c("directory","projectid","sampleid","fid");
    
    field_type=c("text","text","text","text");

    unique_set = c("projectid","sampleid");
    
    create_table(con, table_name, table_fields, field_type, unique_set = unique_set,access_type="restricted");
  }
  
  if(query_type==2){
    table_name = "picrust2_pathway_counts";
    
    table_fields=c("pwid","sampleid","cpm");
    
    field_type=c("text","text","real");
    
    
    unique_set = c("pwid","sampleid");
    
    create_table(con, table_name, table_fields, field_type, unique_set = unique_set,access_type="restricted");
  }

  
  if(query_type==3){
    table_name = "metacyc_pathway_name";
    
    table_fields=c("pwid","pw_name");
    
    field_type=c("text","text");
    
    
    unique_set = c("PWID");
    
    create_table(con, table_name, table_fields, field_type, unique_set = unique_set,access_type="restricted");
  }
  
  if(query_type==4){
    table_name = "metacyc_pathway_ontology";
    
    table_fields=c("pwid","l1","l2","l3","l4",
                   "l5","l6","l7","l8","l9",
                   "l10","l11","l12","l13","l14",
                   "l15","l16");
    
    field_type=c("text","text","text","text","text",
                 "text","text","text","text","text",
                 "text","text","text","text","text",
                 "text","text");
    
    
    unique_set = c("pwid");
    
    create_table(con, table_name, table_fields, field_type, unique_set = unique_set,access_type="restricted");
  }
  
  if(query_type==5){
    table_name = "metacyc_pathway_name";
    
    table_fields=c("pwid","pw_name");
    
    field_type=c("text","text");
    
    
    unique_set = c("PWID");
    
    create_table(con, table_name, table_fields, field_type, unique_set = unique_set,access_type="restricted");
  }
  
  
  if(query_type==6){
    table_name = "humann_shotgun_path_abundance_cpm_unstratified";
    
    table_fields=c("full_id","pathway", "pabun_cpm");
    
    field_type=c("text","text", "real");
    
    
    unique_set = c("full_id","pathway");
    
    create_table(con, table_name, table_fields, field_type, unique_set = unique_set,access_type="restricted");
  }
  
  if(query_type==7){
    table_name = "humann_shotgun_path_abundance_cpm_stratified";
    
    table_fields=c("full_id","pathway_spp", "pabun_cpm");
    
    field_type=c("text","text", "real");
    
    
    unique_set = c("full_id","pathway_spp");
    
    create_table(con, table_name, table_fields, field_type, unique_set = unique_set,access_type="restricted");
  }
  
  if(query_type==8){
    table_name = "humann_shotgun_ko_cpm_unstratified";
    
    table_fields=c("full_id","ko", "ko_cpm");
    
    field_type=c("text","text", "real");
    
    
    unique_set = c("full_id","KO");
    
    create_table(con, table_name, table_fields, field_type, unique_set = unique_set,access_type="restricted");
  }
  
  if(query_type==9){
    table_name = "humann_shotgun_ko_cpm_stratified";
    
    table_fields=c("full_id","ko_spp", "ko_cpm");
    
    field_type=c("text","text", "real");
    
    
    unique_set = c("full_id","ko_spp");
    
    create_table(con, table_name, table_fields, field_type, unique_set = unique_set,access_type="restricted");
  }
  
  if(query_type==10){
    table_name = "test9";
    
    table_fields=c("id","value");
    
    field_type=c("integer","text");
    
    
    unique_set = c("id");
    
    create_table(con, table_name, table_fields, field_type, unique_set = unique_set,access_type="restricted");
  }
  
}
#### !!!!! all characters will be converted to non captital letters 


create_table_prepare <- function(d_set_input){
  #This function prepares fields for creating a new data table.
  
  field_name_str = gsub("\\.","_",tolower(colnames(d_set_input)));
  n_fields = length(field_name_str); 
  cat("table_fields=c(");
  n_loops = ceiling( (n_fields-1)/5);
  for(i in 1:(n_loops)){
    i_start=5*(i-1) +1;
    i_end= min(5*(i), n_fields-1);
    if(i==1){
      cat("\"");
    }else{
      cat("\n\"");
    }
    cat(field_name_str[i_start:i_end],sep="\",\"");
    cat("\",");
  }
  cat(sprintf("\"%s\");",field_name_str[n_fields]));
  cat("\n");
  
  #Printing field_type_str
  field_type_str = rep("NA", n_fields);
  for(i in 1:n_fields){
    cur_class = class(d_set_input[[i]]);
    
    if(cur_class=="factor"){
      field_type_cur="text";
    }
    
    if(cur_class=="numeric" | cur_class=="integer"){
      if(all(floor(d_set_input[[i]])==d_set_input[[i]]) &
         !is.na(all(floor(d_set_input[[i]])==d_set_input[[i]]))){
        field_type_cur="integer";
      }
      else{
        field_type_cur="real";
      }
    }
    if(cur_class=="Date"){
      field_type_cur="date";
    }
    
    if(cur_class=="character"){
      field_type_cur="text";
    }
    
    if(cur_class=="logical"){
      field_type_cur="boolean";
    }
    
    field_type_str[i] = field_type_cur;
  }
  
  cat("\n")
  
  cat("field_type=c(");
  
  for(i in 1:n_loops){
    i_start=5*(i-1) +1;
    i_end= min(5*(i), n_fields-1);
    if(i==1){
      cat("\"");
    }else{
      cat("\n\"");
    }
    cat(field_type_str[i_start:i_end],sep="\",\"");
    cat("\",");
  }
  cat(sprintf("\"%s\");",field_type_str[n_fields]));
  cat("\n");
  
  cat("\n");
  #Printing d_set part in upload_data...
  cat("d_set=data.frame(");
  field_name_raw = colnames(d_set_input);
  for(i in 1:(n_fields-1)){
    cat(sprintf("%s=d_set_input$%s,\n",field_name_str[i],field_name_raw[i]))
  }
  cat(sprintf("%s=d_set_input$%s);",field_name_str[n_fields],field_name_raw[n_fields]))
  
}

library(tidyverse)
tibble(id = seq(1,3), value = c('a','b','c')) %>% write_csv('test.csv')
d_set_input <- read_csv('test.csv')
create_table_prepare(d_set_input)
create_table_type(10)
