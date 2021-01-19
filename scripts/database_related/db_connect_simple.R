#This file is created to connect to database.
# read username and password from a comma-separated file with headers "user, "
#
#

#Aug/09/2019
#We had the `maximum connections` capacity reached and couldn't connect anymore. Many of the connections of from user `postgres`. 
#I solved it by forcing to terminate all connection:
#SELECT *, pg_terminate_backend(pid) FROM pg_stat_activity;
#
#After connection, use this to see list of tables:
#     dbListTables(con)
#Use this to load a table
#My_table = dbReadTable(con,TABLE_NAME)
#
#
# Use the this function in the attached script to connect:
#   connect_database(username, password);
# At your first usage, change the password:
#   dbSendQuery(conn = con,"alter role YOUR_USERNAME password 'YOUR_PASSWORD'”);
# Use this function to list table from database:
#   list_table_from_database()
# You can use this function to load a table:
#   get_table_from_database()
#


library(RPostgreSQL);
library(data.table);

connect_database <- function(config_file="~/dbConfig.txt"){
  
  config = fread(config_file);
  
  drv <- dbDriver("PostgreSQL");
  #print("if you see `expired PostgreSQLDriver error, try to remove object `drv`");
  if(exists("drv")) {
  	for (con2 in dbListConnections(drv)) { dbDisconnect(con2)};
  }
  

  # loads the PostgreSQL driver
  # creates a connection to the postgres database
  # note that "con" will be used later in each connection to the database
  con <<- dbConnect(drv, dbname = "microbiome",
                   host = "plvglover1.mskcc.org", port = 5432,
                   user = config$user, password = config$pass);
  
  if(config$pass=="test123456"){
    print("You are using the default password `test123456`")
    print(prompt="Please, update your password: DO NOT USER YOUR MSK PASSWORD!")
    new_pass = readline(prompt = "type your new password:")
    pw_update_line = sprintf("alter role %s password '%s'",
                             config$user,
                             new_pass);
    dbSendQuery(conn = con,pw_update_line)
    config$pass=new_pass;
    write.table(config,config_file,row.names = F,quote = F,sep=",");
    print("Your password was updated! Config file was changed accordingly")
  }
}

get_table_from_database <- function(table_name){
  
  tb = dbReadTable(con,table_name);
  tb = data.table(tb)
  assign(table_name, tb, envir = .GlobalEnv );
}

list_table_from_database <- function(pattern=NULL){
  
  tb = dbListTables(con);
  
  if(!is.null(pattern)){
    tb = tb[grep(pattern,tb,ignore.case = T)];
  }
  return(tb);
  
}

get_table_from_database_predefined_filter <- function(table=NA, reload=T){
  #This function filters only the sample with most counts for asv_counts_ag and asv_alpha_diversity_ag 
  #for cases the same sample was sequenced more than once. 
   
  query_set = c("asv_counts_ag", "asv_alpha_diversity_ag");
  if(is.na(table) | !(table %in% query_set) ){
    cat(sprintf("\nA predefined filtered table is available for:\n"))
    cat(query_set,sep = "\n");
  }
  
  if(exists(table) & reload==F){
    return();
  }
  
  if(table=="asv_counts_ag"){
    get_table_from_database("asv_counts_ag");
    oligos_id_filtered = unique(asv_counts_ag[count_total>200][order(sampleid,-count_total)][!duplicated(sampleid)]$oligos_id);
    asv_counts_ag = asv_counts_ag[count_total>200][order(sampleid,-count_total)][oligos_id %in% oligos_id_filtered]
    asv_counts_ag[, count_relative := count/count_total];
    asv_counts_ag$filtered_for_highest_coverage_run=T;
    asv_counts_ag <<- asv_counts_ag;
    
  }else if(table=="asv_alpha_diversity_ag"){
    get_table_from_database("asv_alpha_diversity_ag");
    oligos_id_filtered = unique(asv_alpha_diversity_ag[count_total>200][order(sampleid,-count_total)][!duplicated(sampleid)]$oligos_id);
    asv_alpha_diversity_ag = asv_alpha_diversity_ag[count_total>200][order(sampleid,-count_total)][oligos_id %in% oligos_id_filtered]
    asv_alpha_diversity_ag$filtered_for_highest_coverage_run=T;
    asv_alpha_diversity_ag <<- asv_alpha_diversity_ag;
  }
  
  print(sprintf("table %s is loaded and filtered for duplicates. Only the replicate of highest coverage is retained.",
                table));
  
  
}


#This disconnect all connections: for (con in dbListConnections(drv)) { dbDisconnect(con)}
