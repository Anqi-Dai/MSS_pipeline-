library(DBI)
library(RPostgreSQL)
drv <- dbDriver("PostgreSQL");
sessionInfo()
pFile <- read.table("~/dbConfig.txt", header = TRUE, sep = ",", colClasses = "character")

print("Connecting with good credentials")

con <- dbConnect(drv, dbname = "microbiome",
                 host = "plvglover1.mskcc.org", port = 5432,
                 user = pFile$user, password = pFile$pass)
print("connecting with bad credentials")
con <- dbConnect(drv, dbname = "microbiome",
                 host = "plvglover1.mskcc.org", port = 5432,
                 user = pFile$user, password = "not the password")
