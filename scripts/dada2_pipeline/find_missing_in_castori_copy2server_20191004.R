#Jan/23/2020
#
#I am changing this script to look use `asv_alpha_diversity_ag` to check if files were sequenced.
#
#nodate
#
#Check samples id in castori that are not in our database

source('~/projects/human_microbiota/library/sql/get_data_from_query_OTU.R');

#path_castori_mounted="/Volumes/castoricenter";
path_castori_mounted="/tmp/castoricenter";

#require(xlsx);
# sequenced_castori = read.xlsx("/Volumes/castoricenter/Human.Sequencing.Data/Sequenced.BMT.xlsx", sheetIndex = 1 ); #Information for sequenced samples
if(!exists("sequenced_castori")){
  #sequenced_castori = read.xlsx("/Volumes/castoricenter/Human.Sequencing.Data/Sequenced.BMT.xlsx", sheetIndex = 1 ); #Information for sequenced samples
  sequenced_castori = openxlsx::read.xlsx(sprintf("%s/Human.Sequencing.Data/Sequenced.BMT.xlsx",path_castori_mounted),
                                          sheet=1);#sheetIndex = 1 ); #Information for sequenced samples
  sequenced_castori = sequenced_castori[,1:8];
  sequenced_castori = data.table(sequenced_castori);
}
# q11 = get_data_from_query_OTU(11); #Information for Castori samples;
# q11$sampleid = q11$Sample_ID
q11 = get_table_psql("samples_castori_ag");
get_table_psql("asv_alpha_diversity_ag");#q12 = get_data_from_query_OTU(12); #Alpha diversity table contains the samples that are sequenced.

sample_id_not_analyzed_yet = q11[!sampleid %in% asv_alpha_diversity_ag$sampleid]$sampleid; #Define sample_ids in castori that are not in postgreSQl database.

#sample_id_not_analyzed_yet = q11$Sample_ID; #Selecting all samples;

sequenced_castori$sample_in_DB = F;
sequenced_castori$sample_in_DB[!sequenced_castori$Sample_ID %in% sample_id_not_analyzed_yet ] = T;
         
sequenced_castori$Pool.Number = as.character(sequenced_castori$Pool.Number);
#Check pools that may have not been submitted to OTU call pipeline
pending_pools = unique(sequenced_castori[Sample_ID %in% sample_id_not_analyzed_yet &
                                           sequenced_castori$Miseq.V4.V5.sequenced=="yes"]$Pool.Number);


#Generate a table that counts number of sample ids in `pools` that are not (FALSE) or in (TRUE) postgreSQl database.
t = table(sequenced_castori[sequenced_castori$Pool.Number %in% pending_pools,
                            c("Pool.Number","sample_in_DB")]);

list_of_missing_pools = sort(as.numeric(names(t[t[,2] < t[,1],1]))); #I assume the pools with more samples not in database than in database to be pending pools. I make this assumption because some samples may have been sequenced more than once and could be part of more than a single pool. In fact, 936 samples were sequenced more than once.

#stop();

copy_pool_from_castori <- function(pool_number, to_cluster=F){
  
  castori_path=sprintf("%s/Human.Sequencing.Data/Miseq/./Sample_pool%d_complete",
                       path_castori_mounted,
                       pool_number);
  if(!dir.exists(castori_path)){
    castori_path=sprintf("%s/Human.Sequencing.Data/Miseq/Sample_pool%d_complete_HiSeq",
                         path_castori_mounted,
                         pool_number);
    if(!dir.exists(castori_path)){
      print(pool_number);
      stop(sprintf("path %s does not exist in castori center",castori_path));
    }
  }
  
  #output_path="/Volumes/vandenbrinklab/deep_sequencing/Human_data_Mar142018_update/";
  if(to_cluster==F){
    output_path="/Volumes/vandenbrinklab/deep_sequencing/Human_data_castori_update/";
  }else{
    output_path="lilac.mskcc.org:/data/brinkvd/gomesa/e63data/pipeline_16S_call/Human_data_castori_update/";
  }
  sync_str = sprintf("rsync -vrtgoD -R %s/ %s",castori_path, output_path);
  if(0){
    #I made this comment to force copying oligos file.
    castori_last_path = tail(strsplit(castori_path,"/")[[1]],1)
    output_path=sprintf("%s/%s",output_path,castori_last_path);
    sync_str = sprintf("scp %s/*.oligos %s",castori_path, output_path);
  }
  print(sync_str);
  system(sync_str);
  return(sync_str);
}

copy_straight_to_cluster=T;
for( i in 1:length(list_of_missing_pools)){
  #Missing pool 585, 652 (it is hiseq)
  copy_pool_from_castori(list_of_missing_pools[i],to_cluster = copy_straight_to_cluster);
}