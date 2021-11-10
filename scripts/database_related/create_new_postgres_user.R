#Oct/22/2019
#
#
#Create new user (read only) accounts for postgres database.
#https://tableplus.com/blog/2018/04/postgresql-how-to-create-read-only-user.html
#
#
#Each user account will be created with a temporary password that will expire in one week from its creation.
#
#
#In order to set a new password, users need to run the following command:
#dbSendQuery(conn = con,"alter role YOUR_USERNAME password 'YOUR_PASSWORD'");
#
#
#New users will be granted View permission for all tables.
#
#
#Dropping user:
#  drop owned by test; drop role test; (https://stackoverflow.com/questions/9840955/postgresql-drop-role-fails-because-of-default-privileges)
#
#
library(DBI)

source('get_data_from_query_OTU.R');

create_new_postgres_user <- function(user_name, temp_pass){
  
  user_exist = dbGetQuery(con, sprintf("SELECT 1 FROM pg_roles WHERE rolname='%s'",user_name))
  if(length(user_exist$`?column?`) == 1){
    # if(length(user_exist)==1){
    print(sprintf("user %s already exists", user_name));
    return();
  }
  expiration_time = Sys.Date()+3;
  create_user_str = sprintf("create user %s password '%s' valid until '%s'",
                            user_name,
                            temp_pass,
                            expiration_time);
  
  dbSendQuery(con, create_user_str);
  # 1
  
  grant_connect_str = sprintf("GRANT CONNECT ON DATABASE microbiome TO %s", user_name);
  
  dbSendQuery(con, grant_connect_str);
  

  grant_schema_acess = sprintf("GRANT USAGE ON SCHEMA public TO %s;",user_name);
  dbSendQuery(con, grant_schema_acess);
  
  grant_table_read_access = sprintf("GRANT SELECT ON ALL TABLES IN SCHEMA public TO %s\n",user_name);
  dbSendQuery(con, grant_table_read_access);
  
  #Make default privilege to user:
  grant_table_read_access = sprintf("ALTER DEFAULT PRIVILEGES IN SCHEMA public GRANT SELECT ON TABLES TO %s;",
                                      user_name);
  
  dbSendQuery(con, grant_table_read_access);
  # 
  cat(sprintf("user: %s\ntemporary password: %s\npass expiration: %s\n",
              user_name,
              temp_pass,
              expiration_time));
  
}

block_user_who_did_not_change_password <- function(user, temp_password){
  #If user uses temp_password, the user is assigned a new random password (which will lock their account);

  temp_password = 'test123456'
  a=as.numeric(Sys.time())
  set.seed(a)
  
  x=runif(1)
  random_password = sprintf("%2.32f",x);
  
  con_temp = tryCatch(dbConnect(drv, dbname = "microbiome",
                                host = "plvglover1.mskcc.org", port = 5432,
                                user = user, 
                                password = temp_password),
                      error = function(e) e)
  
  if(class(con_temp)[1]=="simpleError"){
    print(sprintf("User %s has changed her/his password!",user));
    remove_user_expiration_query = sprintf("alter user %s VALID UNTIL 'infinity'",
                                           user);
    dbSendQuery(con,remove_user_expiration_query);
    return(0);
  }else{
    print(sprintf("User %s did not change her/his password!",user));
    return(0);
    
    ##I am not changing user password to something arbitrary. I am just keeping it expired!
    block_users_query = sprintf("alter user %s password '%s'",
                                user,
                                random_password);
    dbSendQuery(con_temp,block_users_query);
    dbDisconnect(con_temp);
    
    warning(sprintf("Password for user %s is shuffled! Make sure they change password!",user));
  }
  
}

