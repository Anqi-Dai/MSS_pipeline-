#Oct/22/2019
#
#Create username specific access to database in VDB lab. 
#It removes the general `guest` account.
#
#
#Dropping user:
#  drop owned by test; drop role test; (https://stackoverflow.com/questions/9840955/postgresql-drop-role-fails-because-of-default-privileges)
#
#In order to set a new password, users need to run the following command:
#dbSendQuery(conn = con,"alter role YOUR_USERNAME password 'YOUR_PASSWORD'");
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# when there is request for new connection, add the new username to the list downbelow 
# after the person changed the password, run this script again to unlock this user. then he should be good to go.
# (if the person doesn't change pw within 3 days), run : 
#dbSendQuery(conn = con, "drop owned by <user>; drop role <user>") 
#to drop the user 

# (run this script and look down below to check who changed and who didn't  )
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(DBI)
source('create_new_postgres_user.R');
 
vdb_users <- c("peledj",
               "burgosdm",
               "docampom",
               "slingerj",
               "nguyenc1",
               "daia1",
               "markeyk",
               "armijog",
               "miltiado",
               "andrlovh",
               "lindners",
               "breretod",
               "giardinp",
               "khann2",
               "zuanellc",
               "hatfielp",
               "clurmana",
               "smithm4",
               #"test_user"
               "adintorp", 
               "ngb",
               "chens8",
               "elkriefa",
               "mab4025",
               "harfordl",
               "jesusfar",
               "crossj",
               'funnellt',
               "watersn",
               "maclachlan",
               "feit1",
               "rajs");
#               "gomesa_cluster");

collaboration_users <- c("shouvalr",
                        "lesokhia",
                        "rollingt", 
                        "chowelld",
                        "krishnac");

#temporary_pw=my.name <- readline(prompt="Enter temporary password: ")
temporary_pw='test123456';

vdb_and_collaborator_users <- c(vdb_users,
                                collaboration_users);
#stop("Stopping before creating user");

for(i in 1:length(vdb_and_collaborator_users)){
  create_new_postgres_user(vdb_and_collaborator_users[i],temporary_pw);
}
cat("\n\n")
#create_new_postgres_user(vdb_users[length(vdb_users)],temporary_pw);
#dbDisconnect(con)
if(1){
  #Run this part to remove passwords!
  for(i in 1:length(vdb_and_collaborator_users)){
  block_user_who_did_not_change_password(vdb_and_collaborator_users[i],
                                         temporary_pw)
  }
}

#block_user_who_did_not_change_password(vdb_users[length(vdb_users)],temporary_pw);