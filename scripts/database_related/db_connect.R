#This file is created to connect to database.
# read username and password from a comma-separated file with headers "user, "
#
#

#Aug/09/2019
#We had the `maximum connections` capacity reached and couldn't connect anymore. Many of the connections of from user `postgres`. 
#I solved it by forcing to terminate all connection:
#SELECT *, pg_terminate_backend(pid) FROM pg_stat_activity;
#
setwd('~/MSK/work/microbiome_db/SQL/scripts/')
library(RPostgreSQL);

if(!exists("IS_DBCONNECTED")) {
  
  drv <- dbDriver("PostgreSQL");
  #print("if you see `expired PostgreSQLDriver error, try to remove object `drv`");
  if(exists("drv")) {
  	for (con in dbListConnections(drv)) { dbDisconnect(con)};
  }
  pFile <- read.table("../config/dbConfig.txt", header = TRUE, sep = ",", colClasses = "character")
  
  # loads the PostgreSQL driver
  # creates a connection to the postgres database
  # note that "con" will be used later in each connection to the database
  con <- dbConnect(drv, dbname = "microbiome",
                   host = "plvglover1.mskcc.org", port = 5432,
                   user = pFile$user, password = pFile$pass)
  rm(pFile) # removes the password
  
  IS_DBCONNECTED=T;
}


#This disconnect all connections: for (con in dbListConnections(drv)) { dbDisconnect(con)}


