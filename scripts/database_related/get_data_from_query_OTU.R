#May/17/2017
#
#This script creates a pre-defined set of queries that I will use to analyze 16S data.
#
#0. Get a table according to arg[1] input;.
#   arg[1] : table_name;
#   0.1: Get column names from a table
#   0.2: Get maximum `key` value from a table.
#1. Get samples table
#2. Link FMT, OTU and counts.
#3. Get data from `counts_ag`;
#    3.1 get counts joined with OTU.
#    3.11 Get castori + count + OTU information (count_ag + samples_castori_ag + otu_ag + patients_ag);
#    3.12 Get duke + count + OTU information (count_ag + samples_castori_ag + otu_ag + patients_ag);
#    3.13 Get Regensburg + count + OTU information (count_ag + samples_castori_ag + otu_ag + patients_ag);
#    3.3: Get sample_id, count_group and count_total. Count_group counts total for selected input key.
#           input: key_set_group.
#4. Patient + sample table.
#
#8. UGI vs LGI file (from Doris and Djamila);
#
#11. Get castori data from `access database` (`(...)/castori_access.R`) ;
#   11.1 Get samples with day_relative_to_HCT (first HCT).
#13: Get data prepared for survival analysis (requires 13.1 and 13.2).
#   13.1 : (PREPARING DATA)
#    13.11 : Prepare data for MSK only
#    13.12 : Regensburg only
#    13.13 : Duke only
#   13.2 : (PREPARING DATA)

source('~/MSK/work/microbiome_db/SQL/scripts/db_connect.R'); #Initialize connection.
library("data.table"); #I downloaded this to compute cummulative count; 
#https://stackoverflow.com/questions/18925600/r-cumulative-sum-by-condition

# dbListTables(con) %>% data.frame() %>% filter(grepl("_ag",.))


get_data_from_query_OTU_resource <- function() {
  if(exists("IS_DBCONNECTED")){
    rm(IS_DBCONNECTED, pos = ".GlobalEnv");
  }
  source('~/MSK/work/microbiome_db/SQL/scripts/db_connect.R'); #Initialize connection.
}


get_data_from_query_OTU = function(query_id, ...) {
  
  query = '';
  
  args = list(...);
  
  if(query_id==-1){
    #query is input in args[1];
    
    query = args[[1]];
  }
  
  if(query_id==0) {
    
    table_name = args[1];
    q = postgresqlReadTable(con, name=table_name); #Get `table` directly, without need of a query.
    
  }
  
  if(query_id==0.1){
    #get column names from a table;
    table_name = args[1];
    
    query= sprintf("select column_name, data_type  from information_schema.columns where table_name='%s'",table_name);
  }

  if(query_id==0.2){
    #get maximum `key` value from a table;
    table_name = args[1];
    
    query= sprintf("select max(key) from %s",table_name);
  }
  
  if(query_id==0.3){
    #get counts from a table;
    table_name = args[1];
    
    query= sprintf("select count(*) from %s",table_name);
  }
  
  if(query_id==0.5 | query_id=="showtables"){
    query="SELECT tablename FROM pg_catalog.pg_tables where tablename like '%_ag' order by tablename";
  }
  
  if(query_id==1) {
    #Get all distinct `(sample_id, source) and total counts` from count table.
    
    query = "select cag.sampleid, cag.source_table, sum(cag.count) as total_coverage, count(cag.otu_key) as n_otus from counts_ag as cag group by cag.sampleid, cag.source_table";

    #table_name = "samples";
    #q = get_data_from_query_OTU(0,table_name);
    
  }
  
  if(query_id==2){
    
    query_cur = "select o.*, ob.name, ob.unique_name from (otu_ag as o left join otu_blast_ag as ob on o.otustring=ob.otustring ) order by o.key";
    q = dbGetQuery(con, query_cur);
    q$species_blast = q$unique_name;
    q$genus_blast =sapply(q$unique_name, function(x) {strsplit(x, "_")[[1]][[1]]});
    q$species_gg_and_blast = q$species;
    ind_not_gg_and_with_blast = q$species_gg_and_blast=="" & !is.na(q$species_blast);
    q$species_gg_and_blast[ind_not_gg_and_with_blast]=q$species_blast[ind_not_gg_and_with_blast];
    
    q$genus_gg_and_blast = q$genus;
    ind_not_gg_and_with_blast = q$genus=="" & !is.na(q$species_blast);
    q$genus_gg_and_blast[ind_not_gg_and_with_blast]=q$genus_blast[ind_not_gg_and_with_blast];
    
  }
  
  if(query_id=="2_str"){
    
    
    q2 = get_data_from_query_OTU(2); #Get `OTU_key` information;
    q2 = q2[order(q2$key),];
    
    tax_level_keys = c("k", "p", "c", "o", "f", "g", "s", "otu"); #This indicates the `level` in which data is being grouped;
    tax_level_fields = c("kingdom", "phylum", "class", "ordr", "family", "genus_gg_and_blast", "species_gg_and_blast", "otustring");
    q = q2#otu_key ;#q2;
    #Convert value from q2 to `tax_level` str.
    for(i in 1:length(tax_level_fields)){
      
      for(j in 1:i){
        value = q[[tax_level_fields[ j ] ]];
        value[value==""]="NA";
        tax_level_str_cur = paste(tax_level_keys[j],
                                  value, #q[[tax_level_fields[ j ] ]],
                                  sep = "__");
        if(j==1){
          tax_level_str = tax_level_str_cur;
        }else{
          tax_level_str = paste(tax_level_str, tax_level_str_cur,sep =  "|");
        }
      }
      tax_str_label = sprintf("%s_full",tax_level_keys[j])
      q[[tax_str_label]] = tax_level_str;
      
    }
    
  }
  
  if(query_id==3){
    #select `counts_ag` and includes a column for counts total per sampleid.
    query = "select * from counts_ag as c join (select sampleid, sum(count) as count_total from counts_ag group by sampleid) as ct on c.sampleid=ct.sampleid;";
    
  }
  
  if(query_id==3.01 | query_id=="3_subset") {
    #3.01. Select data from count in which sample_ids as listed in where_query;
    
    sampleid_subset = c(as.character(args[[1]]));# "THISISARANDOMSTRINGTHATSHOULDNOTEXISTsas2323342");
    sampleid_subset_str = paste(sprintf("'%s'",sampleid_subset),collapse = ",");
    where_query= sprintf(" where  c.sampleid in (%s)",sampleid_subset_str);
    counts_total_query =sprintf("select sampleid, sum(count) as count_total from counts_ag as c %s group by sampleid",
                                where_query); 
    query_3p01 = sprintf("select * from (counts_ag as c join otu_ag as o on c.otu_key = o.key) join (%s) as ct on c.sampleid=ct.sampleid %s ",
                    counts_total_query,
                    where_query);

    q = dbGetQuery(con, query_3p01);
    q$count_relative = q$count/q$count_total;
   # print (query);
  }
  
  if(query_id==3.02 | query_id=="count_total"){
    #select total count per sample_id
    query = "select sampleid, sum(count) as count_total from counts_ag group by sampleid";
    
  }
  
  if(query_id==3.1) {
    
    query = sprintf("select * from counts_ag as c join otu_ag as o on c.otu_key = o.key");
  }
  
  if(query_id==3.11) {
    #Select from sample_castori_ag;
    query=paste(" select ",
                " c.count, c.source_table, c.key as \"count_key\", c.samples_key, c.otu_key,",
                #"s.sampleid, s.patient_uid, s.timepoint,",
                "s.mrn,  s.sampleid, s.datecollection, ",
                "p.dateofbmt",
                #"o.kingdom, o.phylum, o.class, o.ordr, o.family, o.genus, o.species ", # I only need key to get color schema from an otu table.
                " from  (counts_ag as c join (samples_castori_ag as s join patients as p on s.mrn=p.mrn) on c.samples_key=s.key)",
                "join otu_ag as o on c.otu_key = o.key",
                "where c.source_table='samples_castori_ag'",
                sep=" ");
  }
  
  if(query_id==3.12) {
    #Select from sample_duke_ag;
    query=paste(" select ",
                " c.count, c.source_table, c.key as \"count_key\", c.samples_key, c.otu_key,",
                "s.sampleid, s.patient_uid, s.timepoint",
                #"o.kingdom, o.phylum, o.class, o.ordr, o.family, o.genus, o.species ", # I only need key to get color schema from an otu table.
                " from  (counts_ag as c join samples_duke_ag s on c.samples_key=s.key)",
                "join otu_ag as o on c.otu_key = o.key",
                "where c.source_table='samples_duke_ag'",
                sep=" ");
  }
  
  if(query_id==3.13) {
    #Select from sample_regensburg_ag;
    query=paste(" select ",
                " c.count, c.source_table, c.key as \"count_key\", c.samples_key, c.otu_key,",
                "s.sampleid, s.patient_uid, s.timepoint",
                #"o.kingdom, o.phylum, o.class, o.ordr, o.family, o.genus, o.species ", # I only need key to get color schema from an otu table.
                " from  (counts_ag as c join samples_regensburg_ag s on c.samples_key=s.key)",
                "join otu_ag as o on c.otu_key = o.key",
                "where c.source_table='samples_regensburg_ag'",
                sep=" ");
  }
  
  if(query_id==3.2){
    #Aggregates counts per taxonomy.
    #select `counts_ag` and includes a column for counts total per sampleid (as in query `3`) 
    
    q0_otu = get_data_from_query_OTU(0,"otu_ag");
    u = data.frame(unique(q0_otu[,3:9]));
    u$unique_index = 1:length(u$kingdom);
    m = merge(q0_otu, u);
    m = m[order(m$key),]
    q3 = get_data_from_query_OTU(3);
    a=aggregate( data.frame(count_tax=q3$count),
                 by=data.frame(count_total=q3$count_total,
                               m[q3$otu_key,c("kingdom",
                                              "phylum",
                                              "class",
                                              "ordr",
                                              "family",
                                              "genus",
                                              "species",
                                              "unique_index")],
                               sampleid=q3$sampleid),
                 sum)
    q=a;

  }
  
  if(query_id==3.3){
    #Select samples from a pre-defined subset and returns it.
    #(Created originally to get count of a taxonomic subgroup (e.g. Blautia,...))
    #
    otu_group_list = args[[1]]; #args[1] is a vector of otu keys.
    otu_group_list_srt = paste(otu_group_list, collapse = ",");
    
    
    query = paste("select c.sampleid, sum(c.count) as count_total, COALESCE(cg.count_group, 0) as count_group from counts_ag as c left join ",
                  "(select sampleid, sum(count) as count_group from counts_ag ",
                  sprintf(" where otu_key in (%s) ",otu_group_list_srt),
                  "group by sampleid ) as cg",
                  " on c.sampleid=cg.sampleid",
                  " group by c.sampleid,  cg.count_group",
                  sep=" ");
    
    
  }
  
  if(query_id==3.4 | query_id=="3_monodominance"){
    if(!exists("q3")){
      q3 = get_data_from_query_OTU(3);
    }
    q3$count_relative = q3$count/q3$count_total;
    
    q = q3[,.(count_tax = sum(count_relative)),
                         .(sampleid=sampleid,
                           otu_key)]
    
    q = q[q[, .I[count_tax == max(count_tax)],
            .(sampleid)]$V1,];

  }
  
  if(query_id==4){
    
    query = sprintf("select * from patients as p join samples_castori_ag as sc on p.mrn=sc.mrn");
  }
  
  if(query_id==5){
    #5. nutrition_demographics + castory:
    query = sprintf("select * from nutrition_demographics_ag as n join samples_castori_ag as sc on n.mrn=sc.mrn");
  
  }
  
  if(query_id==6){
    #6. antibiotics_ag + castory:
    query = "select * from ( antibiotics_ag as a left join  antibiotics_category_map_ag as acm on a.order_name=acm.order_name)";
  }
  
  if(query_id==6.02){
    #6.02 Medication + category. #At this moment it is querying exclusively `antibacterial` drugs.
    category = 'antibacterial';
    query = sprintf("select * from antibiotics_medication_ag as am join antibiotics_category_map_ag as ac on lower(am.order_name) = lower(ac.order_name) where category='%s'",
                    category);
    
  }
  
  if(query_id==8){
    #Read from table `patient_allo_ag`, remove 2nd and further transplants. Old query used to read from xlsx file (see query_id==8.01 for reference). 
    
    #warning("Fields in query 8 are now low case, but they used to not be. You may need to fix proper reference in case!");
    q0_patient_allo = get_data_from_query_OTU(0,"patient_allo_ag");
    
    q0_patient_allo = q0_patient_allo[order(q0_patient_allo$mrn, q0_patient_allo$hct),];
    q0_patient_allo = q0_patient_allo[grepl("Initial",q0_patient_allo$indication),];
    
    q0_patient_allo$source_simple = ifelse((grepl("Cord", q0_patient_allo$source, ignore.case = T)),"cord",
                                           ifelse((grepl("SBA|CD34", q0_patient_allo$source, ignore.case = T)),
                                                  "TCD", "Adult.unmodified"));
    q=q0_patient_allo;
  }
  
  if(query_id==8.01){
    stop("I use this to back up ")
    #8. Reads `UGI vs LGI file` (from Doris and Djamila).
    require(xlsx);
    
    #patient_UGI_vs_LGI_file="/Volumes/AdultBMT/Attending queries/Peled/Microbiota (16-834)/Microbiota & UGI vs LGI GVHD/16-834 Microbiota & UGI vs LGI\ GVHD_final.xlsx"; #/Volumes/MedShared/AdultBMT/Attending queries/Peled/Microbiota (16-834)/Microbiota & UGI vs LGI GVHD/16-834 Microbiota & UGI vs LGI GVHD_final.xlsx";
    patient_allo_file="/Volumes/vandenbrinklab/deep_sequencing/Clinical Annotation/Chart Pulls from Molly and Tsoni manual chart reviews/pt list 2017Jul17.xlsx";
    patient_allo= read.xlsx(patient_allo_file,sheetIndex = 1);

    patient_allo$MRN = as.numeric(as.character(patient_allo$MRN));
    patient_allo = patient_allo[order(patient_allo$MRN, patient_allo$HCT),];
    patient_allo = patient_allo[!duplicated(patient_allo$MRN),];
    #Fixing the `reading` of dates. It initially read dates as integer and I need to convert them back to proper number.
    #patient_UGI_vs_LGI$GI.Onset = as.Date( as.numeric( as.character(patient_UGI_vs_LGI$GI.Onset )),
    #                                       origin="1899-12-30" );
    #patient_UGI_vs_LGI$Onset = as.Date( as.numeric( as.character(patient_UGI_vs_LGI$Onset )),
    #                                    origin="1899-12-30" );
    
    q=patient_allo;
  }
  
  if(query_id==11){
    #Get info directly from access database server. 
    source('~/projects/human_microbiota/library/sql/castori_access.R');
    q=df_tblsamples;
    q$DateCollection = as.Date(q$DateCollection);
    q$MRN = as.numeric(as.character(q$MRN));
    
    #VV Those samples ids belong to a patient whose consent was revogated VV
    sample_id_withdrawn = c("FMT.0015A","FMT.0015B","FMT.0015C","FMT.0015D","FMT.0015E","FMT.0015F","FMT.0015G","FMT.0015H","FMT.0015I");
    mrn_to_be_removed = unique(q$MRN[q$Sample_ID %in% sample_id_withdrawn]);
    if(!is.null(mrn_to_be_removed)){
      mrn_removed = paste(mrn_to_be_removed,collapse = "; ");
      warning( sprintf("MRN ids found in tblSamples that must be removed: %s",mrn_to_be_removed));
    } 
    q = q[!q$MRN %in% mrn_to_be_removed,];
    
  }
  
  if(query_id==11.1){
    #11.1 Get patient and samples info together.
    #stop("Query 11.1 was made to get 2nd transplant information. I don't need it is necessary anymore, but may remove this stop if it seems important.")
    q_castori = get_data_from_query_OTU(11); #Get information from access database.
    
    q0_patient_allo = get_data_from_query_OTU(0,"patient_allo_ag");
    #q0_patient_allo$MRN = q0_patient_allo$mrn;
    #q0_patient_allo$HCT_date = q0_patient_allo$hct;
    #Adding column for `second transplant day`.
    #We censor date if patient goes to second transplant or in last contact day;
    censor_2nd_transplant_day = data.table(data.frame(mrn=q0_patient_allo$mrn,
                                                      hct=q0_patient_allo$hct));
    censor_2nd_transplant_day = censor_2nd_transplant_day[order(censor_2nd_transplant_day$hct),];
    censor_2nd_transplant_day$transplant_day_count=1;
    censor_2nd_transplant_day[,transplant_day_count :=cumsum(transplant_day_count), , by=list(mrn)];
    m_censor = merge( censor_2nd_transplant_day[censor_2nd_transplant_day$transplant_day_count==1,],
                      censor_2nd_transplant_day[censor_2nd_transplant_day$transplant_day_count==2,],
                      by="mrn",
                      all.x = T);
    if( !( all(colnames(m_censor)[c(2,4)] == c( "hct.x", "hct.y")) ) ){
      stop();
    }
    colnames(m_censor)[c(2,4)] = c("hct", "censor_second_transplant_day");
    
    #Recover `2nd transplant day` from original transplant table. Keeps only first transplant data.
    m_q0_patient_allo_and_2nd_transplant_day = merge(q0_patient_allo, m_censor, by=c("hct","mrn"));
    
    m_castori_and_patient = merge(q_castori, m_q0_patient_allo_and_2nd_transplant_day, by.x ="MRN", by.y="mrn");
    m_castori_and_patient$day_relative_to_hsct = as.numeric(m_castori_and_patient$DateCollection-m_castori_and_patient$hct);
    m_castori_and_patient$day_relative_to_second_transplant = as.numeric(m_castori_and_patient$DateCollection-m_castori_and_patient$censor_second_transplant_day);
    m_castori_and_patient$day_relative_to_second_transplant[is.na(m_castori_and_patient$day_relative_to_second_transplant)] = Inf;
    m_castori_and_patient$day_relative_to_last_contact = as.numeric(m_castori_and_patient$DateCollection-m_castori_and_patient$last_contact_date);
    
    #m_castori_and_patient_a = aggregate(data.frame(day_relative_to_hsct=m_castori_and_patient$day_relative_to_hsct,
    #                                               day_relative_to_second_transplant=m_castori_and_patient$day_relative_to_second_transplant,
    #                                               m_castori_and_patient=m_castori_and_patient$GvHD,
    #                                               day_relative_to_last_contact=m_castori_and_patient$day_relative_to_last_contact),
    #                                    by=data.frame(mrn=m_castori_and_patient$MRN,
    #                                                  sampleid=m_castori_and_patient$Sample_ID),
    #                                    function(x) {x[which.min(abs(x))]});
    
    m_castori_and_patient_a = data.frame(mrn=m_castori_and_patient$MRN,
                                         sampleid=m_castori_and_patient$Sample_ID,
                                         day_relative_to_hsct=m_castori_and_patient$day_relative_to_hsct,
                                         day_relative_to_second_transplant=m_castori_and_patient$day_relative_to_second_transplant,
                                         day_relative_to_last_contact=m_castori_and_patient$day_relative_to_last_contact,
                                         gvhd = m_castori_and_patient$gvhd,
                                         hct=m_castori_and_patient$hct,
                                         DateCollection=m_castori_and_patient$DateCollection,
                                         Grade = m_castori_and_patient$grade);
    
    warning("TO DO: There are multiple sampleids for same collection day. You need to filter it!"); 
    #/Volumes/vandenBrinkLab/deep_sequencing/R\ scripts-vdbserver/Sample_Map_10.R 
    
    m_castori_and_patient_a$pre_or_post=NA;
    m_castori_and_patient_a$pre_or_post[m_castori_and_patient_a$day_relative_to_hsct<(-30)]="Pre_d<-30";
    m_castori_and_patient_a$pre_or_post[m_castori_and_patient_a$day_relative_to_hsct>=(-30) & m_castori_and_patient_a$day_relative_to_hsct<0]="Pre";
    m_castori_and_patient_a$pre_or_post[m_castori_and_patient_a$day_relative_to_hsct==0]="HSCT_day";
    m_castori_and_patient_a$pre_or_post[m_castori_and_patient_a$day_relative_to_hsct<=30 & m_castori_and_patient_a$day_relative_to_hsct>0]="Post";
    m_castori_and_patient_a$pre_or_post[m_castori_and_patient_a$day_relative_to_hsct>30]="Post_d>30";
    
    q = m_castori_and_patient_a;
  }
  
  
  if(query_id==11.11){
    stop("This is back up of previous 11.1 version. 11.1 was updated!")
    #Add to `access database` patient information. This is important to compute `day_relative to HSCT`.
    
    q_castori = get_data_from_query_OTU(11); #Get information from access database.
    
    patient_info_file="/Users/gomesa/projects/human_microbiota/data/patient/pt list 2017Jul17_allo_pt_list.txt";
    patient_info=read.table(patient_info_file,header = T,quote ="", sep="\t");
    warning("This procedure considers only pateints with AlloHSCT. Last update in patient info: 2017/Jul/17");
    patient_info$HCT_date = as.Date(patient_info$HCT, "%d-%b-%y");
    patient_info$Last.Contact_date = as.Date(patient_info$Last.Contact, "%d-%b-%y");

    #Adding column for `second transplant day`.
    #We censor date if patient goes to second transplant or in last contact day;
    censor_2nd_transplant_day = data.table(data.frame(MRN=patient_info$MRN,
                                                      HCT_date=patient_info$HCT_date));
    censor_2nd_transplant_day = censor_2nd_transplant_day[order(censor_2nd_transplant_day$HCT_date),];
    censor_2nd_transplant_day$transplant_day_count=1;
    censor_2nd_transplant_day[,transplant_day_count :=cumsum(transplant_day_count), , by=list(MRN)];
    m_censor = merge( censor_2nd_transplant_day[censor_2nd_transplant_day$transplant_day_count==1,],
                     censor_2nd_transplant_day[censor_2nd_transplant_day$transplant_day_count==2,],
                     by="MRN",
                     all.x = T);
    if( !( all(colnames(m_censor)[c(2,4)] == c( "HCT_date.x", "HCT_date.y")) ) ){
      stop();
    }
    colnames(m_censor)[c(2,4)] = c("HCT_date", "censor_second_transplant_day");
    
    #Recover `2nd transplant day` from original transplant table. Keeps only first transplant data.
    m_patient_info_and_2nd_transplant_day = merge(patient_info, m_censor);

    m_castori_and_patient = merge(q_castori, m_patient_info_and_2nd_transplant_day, by="MRN");
    m_castori_and_patient$day_relative_to_hsct = as.numeric(m_castori_and_patient$DateCollection-m_castori_and_patient$HCT_date);
    m_castori_and_patient$day_relative_to_second_transplant = as.numeric(m_castori_and_patient$DateCollection-m_castori_and_patient$censor_second_transplant_day);
    m_castori_and_patient$day_relative_to_second_transplant[is.na(m_castori_and_patient$day_relative_to_second_transplant)] = Inf;
    m_castori_and_patient$day_relative_to_last_contact = as.numeric(m_castori_and_patient$DateCollection-m_castori_and_patient$Last.Contact_date);
    
    #m_castori_and_patient_a = aggregate(data.frame(day_relative_to_hsct=m_castori_and_patient$day_relative_to_hsct,
    #                                               day_relative_to_second_transplant=m_castori_and_patient$day_relative_to_second_transplant,
    #                                               m_castori_and_patient=m_castori_and_patient$GvHD,
    #                                               day_relative_to_last_contact=m_castori_and_patient$day_relative_to_last_contact),
    #                                    by=data.frame(mrn=m_castori_and_patient$MRN,
    #                                                  sampleid=m_castori_and_patient$Sample_ID),
    #                                    function(x) {x[which.min(abs(x))]});
    
    m_castori_and_patient_a = data.frame(mrn=m_castori_and_patient$MRN,
                                         sampleid=m_castori_and_patient$Sample_ID,
                                         day_relative_to_hsct=m_castori_and_patient$day_relative_to_hsct,
                                         day_relative_to_second_transplant=m_castori_and_patient$day_relative_to_second_transplant,
                                         day_relative_to_last_contact=m_castori_and_patient$day_relative_to_last_contact,
                                         gvhd = m_castori_and_patient$GvHD,
                                         HCT_date=m_castori_and_patient$HCT_date,
                                         DateCollection=m_castori_and_patient$DateCollection,
                                         Grade = m_castori_and_patient$Grade);
    
    warning("TO DO: There are multiple sampleids for same collection day. You need to filter it!"); 
    #/Volumes/vandenBrinkLab/deep_sequencing/R\ scripts-vdbserver/Sample_Map_10.R 
    
    m_castori_and_patient_a$pre_or_post=NA;
    m_castori_and_patient_a$pre_or_post[m_castori_and_patient_a$day_relative_to_hsct<(-30)]="Pre_d<-30";
    m_castori_and_patient_a$pre_or_post[m_castori_and_patient_a$day_relative_to_hsct>=(-30) & m_castori_and_patient_a$day_relative_to_hsct<0]="Pre";
    m_castori_and_patient_a$pre_or_post[m_castori_and_patient_a$day_relative_to_hsct==0]="HSCT_day";
    m_castori_and_patient_a$pre_or_post[m_castori_and_patient_a$day_relative_to_hsct<=30 & m_castori_and_patient_a$day_relative_to_hsct>0]="Post";
    m_castori_and_patient_a$pre_or_post[m_castori_and_patient_a$day_relative_to_hsct>30]="Post_d>30";
    
    q = m_castori_and_patient_a;
  }
  
  
  if(query_id==11.3){
    #Get Duke samples `clean` and use patient information to compute day_relative_to_HCT and also day_relative_to_engrafment.

    stop("I don't see a use for this. It seems a copy of 13.23")
    qm1_dukecleaning = get_data_from_query_OTU(-1,"select sd.sampleid, sd.patient_uid, collection_date, sd.timepoint,  transplant_date,  tp_type, collection_date-transplant_date as day_relative_to_hct from samples_duke_clean_ag as sd join patient_duke_ag pd on sd.patient_uid=pd.pid  order by sd.sampleid, pd.pid, transplant_date");
    qm1_dukecleaning = qm1_dukecleaning[qm1_dukecleaning$tp_type=="Allo",];
    warning("considering allo only");
    qm1_dukecleaning= qm1_dukecleaning[order(qm1_dukecleaning$patient_uid,qm1_dukecleaning$transplant_date),];
    qm1_dukecleaning = qm1_dukecleaning[!duplicated(qm1_dukecleaning$sampleid),];
    qm1_dukecleaning = qm1_dukecleaning[!qm1_dukecleaning$transplant_date=="5555-05-05",]; #Remove the cases without transplant information.
    
    #q11p31 = get_data_from_query_OTU(11.31);
    q = qm1_dukecleaning;
    #query = paste("select sd.sampleid, sd.patient_uid, collection_date, sd.timepoint,  transplant_date, collection_date-transplant_date as day_relative_to_hct from samples_duke_clean_ag as sd join patient_duke_ag pd on sd.patient_uid=pd.pid  order by sd.sampleid, pd.pid, transplant_date");#
    
  }
  
  if(query_id==11.31){
    #get engraftment_day_relative_to_hct to Duke patients;
    query_patient_and_engraftment_day = "select pd.pid, pd.transplant_date, pd.tp_type, pded.transplant_date as transplant_date_as_ed, pded.anc___500 from patient_duke_ag as pd left join patient_duke_engraftment_day_ag as pded on pd.pid=pded.pid and pd.transplant_date=pded.transplant_date";
    q_engraftment_day_duke = get_data_from_query_OTU(-1, query_patient_and_engraftment_day); 
    #if( any(q_engraftment_day_duke$transplant_date[is.na(q_engraftment_day_duke$transplant_date_as_ed)] != "5555-05-05")){
    #  stop("There is some unexpected mismatch between transplant_date in patients_duke and patients_duke_engrafment");
    #} #I changed all `5555-05-05` to NA in database.
    if( any(!is.na(q_engraftment_day_duke$transplant_date[is.na(q_engraftment_day_duke$transplant_date_as_ed)])) ){
        stop("There is some unexpected mismatch between transplant_date in patients_duke and patients_duke_engrafment");
    }
    
    if( any(q_engraftment_day_duke$transplant_date[!is.na(q_engraftment_day_duke$transplant_date_as_ed)] == "5555-05-05")){
      stop("I was not expecting any match case with transplant date `5555-05-05`");
    }
    q_engraftment_day_duke=q_engraftment_day_duke[!is.na(q_engraftment_day_duke$transplant_date_as_ed),]; #filter only cases with `transplant date`.
    q_engraftment_day_duke=q_engraftment_day_duke[q_engraftment_day_duke$tp_type=="Allo",];
    q_engraftment_day_duke = q_engraftment_day_duke[order(q_engraftment_day_duke$pid,q_engraftment_day_duke$transplant_date),]
    q_engraftment_day_duke = q_engraftment_day_duke[!duplicated(q_engraftment_day_duke$pid),];
    q_engraftment_day_duke$engraftment_day_relative_to_hct = q_engraftment_day_duke$anc___500 - q_engraftment_day_duke$transplant_date;
    
    q = q_engraftment_day_duke;
  }
  
  if(query_id==12){
    #Load alpha diverstiy!
    if(0){
      alpha_diversity1_file = "/Volumes/vandenBrinkLab/deep_sequencing/analysis_updated_May252017/alpha_diversity.txt";
      alpha_diversity1 = read.table(alpha_diversity1_file,sep = "\t",header = T,quote = "");
      alpha_diversity2_file = "/Volumes/vandenBrinkLab/deep_sequencing/nutrition_related/alpha_diversity_nutrition_Jul242017.txt"; #This is alpha_diversity for case we have nutrition data
      alpha_diversity2 = read.table(alpha_diversity2_file,sep = "\t",header = T);
      alpha_diversity = rbind(alpha_diversity1[,c("X","simpson_reciprocal", "shannon")],alpha_diversity2[,c("X","simpson_reciprocal", "shannon")]);
      alpha_diversity = aggregate(alpha_diversity[,c("simpson_reciprocal", "shannon")],by=list(X=alpha_diversity$X),mean)
      colnames(alpha_diversity)[1] = "sampleid";
      q = alpha_diversity;
    }
    query = "select * from alpha_diversity_ag";
    
  }
  
  if(query_id==13){
    #13 Collecting data for survival analysis
    #Joint data with patient_id (mrn/pid), institution source, sample_id, day relative to HSCT.
    q13p1 = get_data_from_query_OTU(13.1, args); #get patient_id and events;
    q13p2 = get_data_from_query_OTU(13.2); #get patient_id, sample_id, day_relative_to_hct.
    
    q = merge(q13p1, q13p2);
  }
  
  if(query_id==13.01){
    #13 Collecting data for survival analysis
    #Joint data with patient_id (mrn/pid), institution source, sample_id, day relative to HSCT.
    q13p11 = get_data_from_query_OTU(13.11, args); #get patient_id and events;
    q13p21 = get_data_from_query_OTU(13.21); #get patient_id, sample_id, day_relative_to_hct.
    
    q = merge(q13p21,q13p11, 
              by.x="patient_id",
              by.y="mrn");
  }
  
  if(query_id==13.02){
    #13 Collecting data for survival analysis
    #Joint data with patient_id (mrn/pid), institution source, sample_id, day relative to HSCT.
    q13p12 = get_data_from_query_OTU(13.12, args); #get patient_id and events;
    q13p22 = get_data_from_query_OTU(13.22); #get patient_id, sample_id, day_relative_to_hct.
    
    q = merge(q13p22,q13p12, 
              by=c("patient_id","institution"));
  }
  
  if(query_id==13.1){
    #13.1 Collecting data for overall survival analysis (I may update it to consider competing risks).
    # d_set_event (joins patient_id with event (death/cod/death_day) ).
    
    filter_index = args[1][[1]]; #Filter 0 (or null, nothing, length 0): Does not filter anything; `Allo` only is filtered by default;
    q13p11 = get_data_from_query_OTU(13.11,filter_index);
    q13p12 = get_data_from_query_OTU(13.12);
    q13p13 = get_data_from_query_OTU(13.13);
    q13p14 = get_data_from_query_OTU(13.14);
    
    q = rbind( data.frame(patient_id=q13p11$mrn,
                          event_day=q13p11$event_day,
                          event_type=q13p11$event_type,
                          institution="MSK_allo"),
               q13p12[,c("patient_id","event_day","event_type","institution")],
               q13p13[,c("patient_id","event_day","event_type","institution")],
               q13p14[,c("patient_id","event_day","event_type","institution")]);
    
  }
  
  if(query_id== -13.1){
    stop("obsolete! I am using subqueries 13.11, ..., 13.14 to get event_day and event_type per institution.")
    #13.1 Collecting data for survival analysis
    # d_set_event (joins patient_id with event (death/cod/death_day) ).
    #fields:
    #   - patient_id,
    #   - (mrn or uid),
    #   - event_day, 
    #   - event_type, 
    #   - event_competing_day
    #       Same as event_day; but update `relapse` with `relapse day`;
    #   - event_competing_type,
    #       Competing event types:
    #         0 = survival
    #         1 = relapse (Update day for day of relapse)
    #         2 = no engraftment
    #         3 - GVHD related death
    #         4- other
    #
    #   - institution
    
    #filter types: 
        #1. Use only `Allo` cases; Also 
    filter_index = args[1][[1]]; #Filter 0 (or null, nothing, length 0): Does not filter anything; `Allo` only is filtered by default;
                           #bit1: Remove TCD from MSK data.
                           #bit2: 
    
    if(is.null(filter_index)){
      filter_index=3;
      warning("using default filter_index=3. Allo only and remove TCD from MSK");
    }
    filter_index  = as.integer(filter_index);
    
    #Loading all data a beginning
    q0_patient_MSK_allo = get_data_from_query_OTU(8); #get_data_from_query_OTU(0,"patient_allo_ag");
    q0_patient_Regensburg = get_data_from_query_OTU(0,"patient_regensburg_ag");
    q0_patient_Duke = get_data_from_query_OTU(0,"patient_duke_ag");
    
    if(intToBits(filter_index)[1]==1){
      #Allo only is filtered by default. I have this as a check.
      # - Filter for "Allo" only; 
      # - Filter msk for non-TCD cases.
      tp_type_set <- c("Allo", "Auto", "N/A", "not transplanted", "Not transplanted", NA);
      if( !all(q0_patient_Duke$tp_type %in% tp_type_set) | 
          !all(q0_patient_Regensburg$kind_of_transplantation %in% tp_type_set) ){
        stop("Unexpected Tptype. I was expecting Allo or Auto only");
      }
      q0_patient_Regensburg = q0_patient_Regensburg[q0_patient_Regensburg$kind_of_transplantation=="Allo",];
      q0_patient_Duke = q0_patient_Duke[q0_patient_Duke$tp_type=="Allo",];
      cat("select only Allo");
    }
    
    if(intToBits(filter_index)[2]==1){
      #Filter TCD cases out of MSK data;
      q0_patient_MSK_allo$source_simple <- ifelse((grepl("Cord", q0_patient_MSK_allo$source, ignore.case = T)),"cord",
                                                  ifelse((grepl("SBA|CD34", q0_patient_MSK_allo$source, ignore.case = T)),"TCD", "Adult.unmodified"));
      q0_patient_MSK_allo = q0_patient_MSK_allo[q0_patient_MSK_allo$source_simple!="TCD",];
      
      cat("Filter out TCD cases");
      
    }
    
    #q0_patient_MSK_allo = get_data_from_query_OTU(0,"patient_allo_ag");
    
    #Some patients have duplicated field. Select only the first HCT case.
    q0_patient_MSK_allo = q0_patient_MSK_allo[order(q0_patient_MSK_allo$mrn, q0_patient_MSK_allo$hct),];
    q0_patient_MSK_allo = q0_patient_MSK_allo[!duplicated(q0_patient_MSK_allo$mrn),];
    
    q0_patient_MSK_allo$event_day = q0_patient_MSK_allo$last_contact - q0_patient_MSK_allo$hct;
    q0_patient_MSK_allo$event_type = q0_patient_MSK_allo$vital_status=="Dead";

    q0_patient_MSK_allo$event_competing_type = 4; #Initialize with all values equal to 4;
    ind_alive = q0_patient_MSK_allo$vital_status == "Alive";
    q0_patient_MSK_allo$event_competing_type[ind_alive]=0;
    if( ! (sum( is.na(q0_patient_MSK_allo$cod) & !ind_alive)==2 | !any( is.na(q0_patient_MSK_allo$cod) & !ind_alive) )){
      stop("I was expecting two NA for cod field or none (when filtering out TCD cases). If they was updated, this should be fixed.");
    }
    q0_patient_MSK_allo$event_competing_type[grepl("graft",q0_patient_MSK_allo$cod,ignore.case = T)] = 2;
    q0_patient_MSK_allo$event_competing_type[grepl("gvhd",q0_patient_MSK_allo$cod,ignore.case = T)] = 3;
    ind_relapse = q0_patient_MSK_allo$relapse=="Y" | q0_patient_MSK_allo$pod=="Y";
    q0_patient_MSK_allo$event_competing_type[ind_relapse ] = 1; #There are two cases without a labeled `cod` that is relapse;
    
    #Update `event_competing_day` for cases of relapse
    q0_patient_MSK_allo$event_competing_day = q0_patient_MSK_allo$event_day;
    ind_relapse=q0_patient_MSK_allo$event_competing_type==1; 
    pre_competing_day_relapse1 = q0_patient_MSK_allo$relapse_date[ind_relapse] - q0_patient_MSK_allo$hct[ind_relapse];
    pre_competing_day_relapse2 = q0_patient_MSK_allo$pod_date[ind_relapse] - q0_patient_MSK_allo$hct[ind_relapse];
    competing_day_relapse = apply(cbind(pre_competing_day_relapse1, pre_competing_day_relapse2),1, function(x) {min(x, na.rm = T)})
    q0_patient_MSK_allo$event_competing_day[ind_relapse] = competing_day_relapse;
    
    #Regensburg
    #q0_patient_Regensburg = get_data_from_query_OTU(0,"patient_regensburg_ag");
    warning("clean patient_Regensburg_ag table. Remove rows that are NAs in patient_id or hct`")
    q0_patient_Regensburg = q0_patient_Regensburg[!is.na(q0_patient_Regensburg$patient_id),];
    #Some patients have duplicated field. Select only the first HCT case.
    q0_patient_Regensburg = q0_patient_Regensburg[order(q0_patient_Regensburg$patient_id, q0_patient_Regensburg$hct),];
    q0_patient_Regensburg = q0_patient_Regensburg[!duplicated(q0_patient_Regensburg$patient_id),];
    
    ind_is_dead =!is.na(q0_patient_Regensburg$date_of_death__updated_july_2017);
    if(any (!is.na(q0_patient_Regensburg$last_contact_clean[ind_is_dead]) ) ){
      warning("last_contact_clean incompatible with date_of_death_updated_july_2017");
    }
    if(any(is.na(q0_patient_Regensburg$cause_of_death__updated_july_2017[ind_is_dead])) |
       any(!is.na(q0_patient_Regensburg$cause_of_death__updated_july_2017[!ind_is_dead])) ){
      stop("Incompatible alive/death cases have/don't have a value in cod");
    }
    q0_patient_Regensburg$event_day =  q0_patient_Regensburg$date_of_death__updated_july_2017 -
      q0_patient_Regensburg$hct;
    q0_patient_Regensburg$event_day[!ind_is_dead] = q0_patient_Regensburg$last_contact_clean[!ind_is_dead] -
      q0_patient_Regensburg$hct[!ind_is_dead];
    q0_patient_Regensburg$event_type = ind_is_dead; #1 means death, 0 means alive.
    
    q0_patient_Regensburg$event_competing_type = 4; #Initialize with all values equal to 4;
    q0_patient_Regensburg$event_competing_type[!ind_is_dead]=0;
    q0_patient_Regensburg$event_competing_type[grepl("gvhd",q0_patient_Regensburg$cause_of_death__updated_july_2017,ignore.case = T)] = 3;
    q0_patient_Regensburg$event_competing_type[grepl("graft",q0_patient_Regensburg$cause_of_death__updated_july_2017,ignore.case = T)] = 2;
    q0_patient_Regensburg$event_competing_type[q0_patient_Regensburg$relapse_pod=="Yes"] = 1;

    #Update `event_competing_day` for cases of relapse
    q0_patient_Regensburg$event_competing_day = q0_patient_Regensburg$event_day;
    ind_relapse=q0_patient_Regensburg$event_competing_type==1; 
    warning("Relapse day has to be fixed! Some is not there, others are MONTH/YEAR");
    q0_patient_Regensburg$relapse__pod_date[q0_patient_Regensburg$relapse__pod_date=="May 2015"] = "2015-05-15";
    q0_patient_Regensburg$relapse__pod_date[q0_patient_Regensburg$relapse__pod_date=="June 2015"] = "2015-06-15";
    q0_patient_Regensburg$relapse__pod_date[q0_patient_Regensburg$relapse__pod_date=="March 2015"] = "2015-03-15";
    
    q0_patient_Regensburg$event_competing_day = q0_patient_Regensburg$event_day;
    competing_day_relapse = as.Date(q0_patient_Regensburg$relapse__pod_date[ind_relapse], "%Y-%m-%d") - q0_patient_Regensburg$hct[ind_relapse];
    q0_patient_Regensburg$event_competing_day[ind_relapse] = competing_day_relapse;
    
    #Duke
    #q0_patient_Duke = get_data_from_query_OTU(0,"patient_duke_ag");
    #q0_patient_Duke = q0_patient_Duke[q0_patient_Duke$tp_type=="Allo",];
    
    #Some patients have duplicated field. Select only the first HCT case.
    q0_patient_Duke = q0_patient_Duke[order(q0_patient_Duke$pid, q0_patient_Duke$transplant_date),]
    q0_patient_Duke = q0_patient_Duke[!duplicated(q0_patient_Duke$pid),];
    
    q0_patient_Duke$last_contact_inferred = q0_patient_Duke$last_contact;
    q0_patient_Duke$last_contact_inferred[is.na(q0_patient_Duke$last_contact)] = as.Date("08/29/2017","%m/%d/%Y");
    q0_patient_Duke$last_contact_inferred[q0_patient_Duke$last_contact_inferred < q0_patient_Duke$transplant_date] = as.Date("08/29/2017","%m/%d/%Y");
    q0_patient_Duke$last_contact[is.na(q0_patient_Duke$last_contact)] = as.Date("08/29/2017","%m/%d/%Y");
    q0_patient_Duke$event_day = q0_patient_Duke$last_contact_inferred - q0_patient_Duke$transplant_date;
    ind_is_dead = q0_patient_Duke$vital_status=="Dead";
    q0_patient_Duke$event_type = ind_is_dead+0;
    #    q0_patient_Duke = q0_patient_Duke[q0_patient_Duke$tp_type=="Allo",]
    if(any(q0_patient_Duke$cause_of_death=="" & ind_is_dead)  | 
       any(q0_patient_Duke$cause_of_death!="" & !ind_is_dead) |
       any(is.na(q0_patient_Regensburg$cause_of_death__updated_july_2017) ) ){
      warning("Fix cod in DUKE");
      warning("Incompatible alive/death cases have/don't have a value in cod");
      warning("cod will be grouped as `others`");
      warning("Fix cod in DUKE");
    }
    #q0_patient_Duke = q0_patient_Duke[!(q0_patient_Duke$cause_of_death=="" & ind_is_dead),];
    q0_patient_Duke$event_competing_type = 4; #Initialize with all values equal to 4;
    ind_is_dead = q0_patient_Duke$vital_status=="Dead";
    q0_patient_Duke$event_competing_type[!ind_is_dead]=0; #Alive cases
    q0_patient_Duke$event_competing_type[grepl("gvhd",q0_patient_Duke$cause_of_death,ignore.case = T)] = 3;
    q0_patient_Duke$event_competing_type[grepl("graft",q0_patient_Duke$cause_of_death,ignore.case = T)] = 2;
    ind_relapse = !is.na(q0_patient_Duke$transplant_relapse) |
      !is.na(q0_patient_Duke$progressive_disease_date) |
      !is.na(q0_patient_Duke$dbtc_relapse_date) |
      grepl("Recurrent",q0_patient_Duke$cause_of_death,ignore.case = T);
    q0_patient_Duke$event_competing_type[ind_relapse] = 1;
    
    q0_patient_Duke$event_competing_day = q0_patient_Duke$event_day;
    
 
    pre_competing_day_relapse1 = q0_patient_Duke$transplant_relapse[ind_relapse] - q0_patient_Duke$transplant_date[ind_relapse];
    pre_competing_day_relapse2 = q0_patient_Duke$dbtc_relapse_date[ind_relapse] - q0_patient_Duke$transplant_date[ind_relapse];
    pre_competing_day_relapse3 = q0_patient_Duke$progressive_disease_date[ind_relapse] -q0_patient_Duke$transplant_date[ind_relapse];
    competing_day_relapse = apply(cbind(pre_competing_day_relapse1,
                                        pre_competing_day_relapse2,
                                        pre_competing_day_relapse3),
                                  1, 
                                  function(x) {min(x, na.rm = T)})
    
    q0_patient_Duke$event_competing_day[ind_relapse] = competing_day_relapse;
    
    d_set_event = rbind(data.frame(patient_id=q0_patient_MSK_allo$mrn,
                                   event_day=q0_patient_MSK_allo$event_day,
                                   event_type=q0_patient_MSK_allo$event_type,
                                   event_competing_day=q0_patient_MSK_allo$event_competing_day, 
                                   event_competing_type=q0_patient_MSK_allo$event_competing_type,
                                   institution="MSK_allo"),
                        data.frame(patient_id=q0_patient_Regensburg$patient_id,
                                   event_day=q0_patient_Regensburg$event_day,
                                   event_type=q0_patient_Regensburg$event_type,
                                   event_competing_day=q0_patient_Regensburg$event_competing_day,
                                   event_competing_type=q0_patient_Regensburg$event_competing_type,
                                   institution="Regensburg"),
                        data.frame(patient_id=q0_patient_Duke$pid,
                                   event_day=q0_patient_Duke$event_day,
                                   event_type=q0_patient_Duke$event_type,
                                   event_competing_type=q0_patient_Duke$event_competing_type,
                                   event_competing_day=q0_patient_Duke$event_competing_day,
                                   institution="Duke"));
    
    #Add `factor` label for competing events.
    competing_event_groups = c("censor",
                               "relapse/pod",
                               "no engraftment",
                               "GRM_death",
                               "other_death");
    competing_event_groups_as_f = factor(competing_event_groups,
                                         levels=unique(competing_event_groups));
    
    d_set_event$event_competing_type_as_f =competing_event_groups_as_f[d_set_event$event_competing_type+1];
    
    q = d_set_event;
  }
  
  if(query_id==13.11){
    #Colecting data for survival Analysis MSK only!
    
    #filter types: 
    #1. Use only `Allo` cases; Also 
    filter_index = args[1][[1]]; #Filter 0 (or null, nothing, length 0): Does not filter anything; `Allo` only is filtered by default;
    #bit1: Remove TCD from MSK data.
    #bit2: 
    
    if(is.null(filter_index)){
      filter_index=1;
      warning("using default filter_index=1. Allo only and does not remove TCD from MSK");
    }
    filter_index  = as.integer(filter_index);
    
    q0_patient_MSK_allo = get_data_from_query_OTU(8); #get_data_from_query_OTU(0,"patient_allo_ag");
    q0_patient_MSK_allo$source_simple <- ifelse((grepl("Cord", q0_patient_MSK_allo$source, ignore.case = T)),"cord",
                                                ifelse((grepl("SBA|CD34", q0_patient_MSK_allo$source, ignore.case = T)),"TCD", "Adult.unmodified"));
    
    if(intToBits(filter_index)[2]==1){
      #Filter TCD cases out of MSK data;
      #q0_patient_MSK_allo$source_simple <- ifelse((grepl("Cord", q0_patient_MSK_allo$source, ignore.case = T)),"cord",
      #                                            ifelse((grepl("SBA|CD34", q0_patient_MSK_allo$source, ignore.case = T)),"TCD", "Adult.unmodified"));
      q0_patient_MSK_allo = q0_patient_MSK_allo[q0_patient_MSK_allo$source_simple!="TCD",];
      cat("Filter out TCD cases");
      
    }
    
    #q0_patient_MSK_allo = get_data_from_query_OTU(0,"patient_allo_ag");
    
    #Some patients have duplicated field. Select only the first HCT case.
    q0_patient_MSK_allo = q0_patient_MSK_allo[order(q0_patient_MSK_allo$mrn, q0_patient_MSK_allo$hct),];
    q0_patient_MSK_allo = q0_patient_MSK_allo[!duplicated(q0_patient_MSK_allo$mrn),];
    
    q0_patient_MSK_allo$event_day = q0_patient_MSK_allo$last_contact - q0_patient_MSK_allo$hct;
    q0_patient_MSK_allo$event_type = q0_patient_MSK_allo$vital_status=="Dead";
    
    #Event competing
    # 0. survival
    # 1. relapse
    # 2. no engraftment
    # 3. GVHD related mortality
    # 4. other
    
    event_GRM_competing_set = c("censor",
                                "relapse",
                                "no engraftment",
                                "GVHD related mortality",
                                "other");
    
    q0_patient_MSK_allo$event_GRM_competing_type = 4; #Initialize with all values equal to 4;
    ind_alive = q0_patient_MSK_allo$vital_status == "Alive";
    q0_patient_MSK_allo$event_GRM_competing_type[ind_alive]=0;
    if( ! (sum( is.na(q0_patient_MSK_allo$cod) & !ind_alive)<=12 | !any( is.na(q0_patient_MSK_allo$cod) & !ind_alive) )){
      stop("I was expecting one or two `NA`s for cod field or none (when filtering out TCD cases). If this was updated, this should be fixed.");
    }
    q0_patient_MSK_allo$event_GRM_competing_type[grepl("graft",q0_patient_MSK_allo$cod,ignore.case = T)] = 2;
    q0_patient_MSK_allo$event_GRM_competing_type[grepl("gvhd",q0_patient_MSK_allo$cod,ignore.case = T)] = 3;
    ind_relapse = q0_patient_MSK_allo$relapse=="Y" | q0_patient_MSK_allo$pod=="Y";
    q0_patient_MSK_allo$event_GRM_competing_type[ind_relapse ] = 1; #There are two cases without a labeled `cod` that is relapse;
    
    #Update `event_GRM_competing_day` for cases of relapse
    q0_patient_MSK_allo$event_GRM_competing_day = q0_patient_MSK_allo$event_day;
    ind_relapse=q0_patient_MSK_allo$event_GRM_competing_type==1; 
    pre_competing_day_relapse1 = q0_patient_MSK_allo$relapse_date[ind_relapse] - q0_patient_MSK_allo$hct[ind_relapse];
    pre_competing_day_relapse2 = q0_patient_MSK_allo$pod_date[ind_relapse] - q0_patient_MSK_allo$hct[ind_relapse];
    competing_day_relapse = apply(cbind(pre_competing_day_relapse1, pre_competing_day_relapse2),1, function(x) {min(x, na.rm = T)})
    q0_patient_MSK_allo$event_GRM_competing_day[ind_relapse] = competing_day_relapse;
    
    q0_patient_MSK_allo$event_GRM_competing_type_as_f = factor(event_GRM_competing_set[q0_patient_MSK_allo$event_GRM_competing_type+1],
                                                               levels =event_GRM_competing_set);
    
   #Compute competing risk for GVHD.
    # 0 = alive without GVHD nor relapse
    # 1 = clinical diagnosis of GVHD
    # 2 = relapse without GVHD
    # 3 = death without GVHD
    
    d_set_competing_raw = data.frame(censor = as.numeric(q0_patient_MSK_allo$last_contact - q0_patient_MSK_allo$hct),
                                     GVHD_onset = as.numeric(q0_patient_MSK_allo$onset - q0_patient_MSK_allo$hct),
                                     relapse = as.numeric(q0_patient_MSK_allo$relapse_date - q0_patient_MSK_allo$hct),
                                     pod = as.numeric(q0_patient_MSK_allo$pod_date - q0_patient_MSK_allo$hct),
                                     day_of_death=as.numeric(NA));
    
    d_set_competing_raw$day_of_death[q0_patient_MSK_allo$vital_status=="Dead"] = d_set_competing_raw$censor[q0_patient_MSK_allo$vital_status=="Dead"];
    d_set_competing_raw$censor[q0_patient_MSK_allo$vital_status=="Dead"] = Inf;
    d_set_competing_raw$relapse_no_GVHD = d_set_competing_raw$relapse;
    d_set_competing_raw$relapse_no_GVHD[q0_patient_MSK_allo$gvhd=="Y" & !is.na(q0_patient_MSK_allo$gvhd)] = NA;
    d_set_competing_raw$pod_no_GVHD = d_set_competing_raw$pod;
    d_set_competing_raw$pod_no_GVHD[q0_patient_MSK_allo$gvhd=="Y" & !is.na(q0_patient_MSK_allo$gvhd)] = NA;
    d_set_competing_raw[is.na(d_set_competing_raw)] = Inf;
    
    #Doing checks
    if( !all(d_set_competing_raw$relapse[q0_patient_MSK_allo$relapse!="Y"]==Inf) |
        !all(d_set_competing_raw$pod[q0_patient_MSK_allo$pod!="Y"]==Inf) |
        !all(d_set_competing_raw$GVHD_onset[q0_patient_MSK_allo$gvhd!="Y" | is.na(q0_patient_MSK_allo$gvhd)]==Inf) ){
      stop("error! Missing data or annotation for pod/relapse/gvhd is wrong.");
    }
    
    f_event_GVHD_competing_day = function(d_set_competing_raw,hierarchical_risks, hierarchical_risks_value,hierarchical_risks_value_str) {
      #It must contain a field called censor!
      d_set_competing = data.table(d_set_competing_raw[,hierarchical_risks]);
      event_type = apply(as.matrix(d_set_competing) , 1, which.min);
      index_min = apply(d_set_competing, 1, which.min);
      
      event_competing_day = as.data.frame(d_set_competing)[cbind(1:nrow(d_set_competing),index_min)];
      event_competing_type = hierarchical_risks_value[index_min];
      event_competing_type_as_f = factor(hierarchical_risks_value_str[event_competing_type+1], levels = hierarchical_risks_value_str);
      
      return(data.frame(event_competing_day=event_competing_day,
                        event_competing_type=event_competing_type,
                        event_competing_type_as_f=event_competing_type_as_f));
    }
   
    hierarchical_risks_fields = c("censor", "GVHD_onset", "relapse_no_GVHD", "pod_no_GVHD", "day_of_death");
    hierarchical_risks_value=c(0,1,2,2,3);
    hierarchical_risks_value_str = c("censor","GVHD_onset","relapse/pod(noGVHD)","death");
    
    d_set_GVHD_competing = f_event_GVHD_competing_day(d_set_competing_raw,
                                                      hierarchical_risks_fields, 
                                                      hierarchical_risks_value,
                                                      hierarchical_risks_value_str);
    
    q0_patient_MSK_allo$event_GVHD_competing_day = d_set_GVHD_competing$event_competing_day;
    q0_patient_MSK_allo$event_GVHD_competing_type = d_set_GVHD_competing$event_competing_type;
    q0_patient_MSK_allo$event_GVHD_competing_type_as_f = d_set_GVHD_competing$event_competing_type_as_f;
    
    #Compute competing risk for GVHD grade 2-4.
    # 0 = alive without GVHD nor relapse, GVHD grade < 2
    # 1 = clinical diagnosis of GVHD grade 2-4
    # 2 = relapse without GVHD
    # 3 = death without GVHD
    d_set_competing_raw$GVHD_onset_gt1=d_set_competing_raw$GVHD_onset;
    d_set_competing_raw$GVHD_onset_gt1[q0_patient_MSK_allo$grade==1 | is.na(q0_patient_MSK_allo$grade)] = Inf;
    
    hierarchical_risks_fields = c("censor", "GVHD_onset_gt1", "relapse_no_GVHD", "pod_no_GVHD", "day_of_death");
    hierarchical_risks_value=c(0,1,2,2,3);
    hierarchical_risks_value_str = c("censor","GVHD_onset_gt1","relapse/pod(noGVHD)","death");
    
    d_set_GVHD_competing = f_event_GVHD_competing_day(d_set_competing_raw,
                                                      hierarchical_risks_fields, 
                                                      hierarchical_risks_value, 
                                                      hierarchical_risks_value_str);
    q0_patient_MSK_allo$event_GVHDgt1_competing_day = d_set_GVHD_competing$event_competing_day;
    q0_patient_MSK_allo$event_GVHDgt1_competing_type = d_set_GVHD_competing$event_competing_type;
    q0_patient_MSK_allo$event_GVHDgt1_competing_type_as_f = d_set_GVHD_competing$event_competing_type_as_f;
    
    #Compute competing risk for GVHD.
    # 0 = alive without GVHD nor relapse, GVHD grade < 3
    # 1 = clinical diagnosis of GVHD grade 3-4
    # 2 = relapse without GVHD
    # 3 = death without GVHD
    d_set_competing_raw$GVHD_onset_gt2=d_set_competing_raw$GVHD_onset;
    d_set_competing_raw$GVHD_onset_gt2[q0_patient_MSK_allo$grade<=2 | is.na(q0_patient_MSK_allo$grade)] = Inf;

    hierarchical_risks_fields = c("censor", "GVHD_onset_gt2", "relapse_no_GVHD", "pod_no_GVHD", "day_of_death");
    hierarchical_risks_value=c(0,1,2,2,3);
    hierarchical_risks_value_str = c("censor","GVHD_onset_gt2","relapse/pod(noGVHD)","death");
    
    d_set_GVHD_competing = f_event_GVHD_competing_day(d_set_competing_raw,
                                                      hierarchical_risks_fields, 
                                                      hierarchical_risks_value,
                                                      hierarchical_risks_value_str);
    q0_patient_MSK_allo$event_GVHDgt2_competing_day = d_set_GVHD_competing$event_competing_day;
    q0_patient_MSK_allo$event_GVHDgt2_competing_type = d_set_GVHD_competing$event_competing_type;
    q0_patient_MSK_allo$event_GVHDgt2_competing_type_as_f = d_set_GVHD_competing$event_competing_type_as_f;
    
    #Compute competing risk for TRM.
    # 0 = alive
    # 1 = relapse/pod
    # 2 = death
    hierarchical_risks_fields = c("censor", "relapse", "pod", "day_of_death");
    hierarchical_risks_value=c(0,1,1,2);
    hierarchical_risks_value_str= c("censor","relapse/pod","TRM_death");
    
    d_set_TRM_competing = f_event_GVHD_competing_day(d_set_competing_raw,
                                                     hierarchical_risks_fields,
                                                     hierarchical_risks_value, 
                                                     hierarchical_risks_value_str);
    q0_patient_MSK_allo$event_TRM_competing_day = d_set_TRM_competing$event_competing_day;
    q0_patient_MSK_allo$event_TRM_competing_type = d_set_TRM_competing$event_competing_type;
    q0_patient_MSK_allo$event_TRM_competing_type_as_f = d_set_TRM_competing$event_competing_type_as_f;
    
    #Adding engrafment day information:
    q8 = get_data_from_query_OTU(8);
    if(!all(q0_patient_MSK_allo$mrn %in% q8$mrn)){
      print(sprintf("%2.2f of patients are not listed in Doris/Djamila HCT file", 1- mean(q0_patient_MSK_allo$mrn %in% q8$mrn)));
    }
    m = merge(q0_patient_MSK_allo,
              data.frame(mrn=as.numeric(as.character(q8$mrn)),
                         engrafment_day_relative_to_HCT=q8$tt_anc_500),
              all.x=T);
    
    q = m;
    
    
  }
  
  if(query_id==13.12){
    #Colecting data for survival Analysis Regensburg only!
    
    #filter types: 
    #1. Use only `Allo` cases; Also 
    filter_index = args[1][[1]]; #Filter 0 (or null, nothing, length 0): Does not filter anything; `Allo` only is filtered by default;
    #bit1: Remove TCD from MSK data.
    #bit2: 
    
    if(is.null(filter_index)){
      filter_index=3;
    #  warning("using default filter_index=3. Allo only and remove TCD from MSK");
    }
    filter_index  = as.integer(filter_index);
    
    q0_patient_Regensburg = get_data_from_query_OTU(0,"patient_regensburg_ag");
    
    tp_type_set <- c("Allo", "Auto");
    
    if(intToBits(filter_index)[1]==1){
      #Allo only is filtered by default. I have this as a check.
      # - Filter for "Allo" only; 
      # - Filter msk for non-TCD cases.
      tp_type_set <- c("Allo", "Auto");
      if( !all(q0_patient_Regensburg$kind_of_transplantation %in% tp_type_set) ){
        stop("Unexpected Tptype. I was expecting Allo or Auto only");
      }
      q0_patient_Regensburg = q0_patient_Regensburg[q0_patient_Regensburg$kind_of_transplantation=="Allo",];
      cat("select only Allo");
    }
    
    if(intToBits(filter_index)[2]==1){
      #Filter TCD cases out of MSK data;
    }
    
    
    #Regensburg
    #q0_patient_Regensburg = get_data_from_query_OTU(0,"patient_regensburg_ag");
    warning("clean patient_Regensburg_ag table. Remove rows that are NAs in patient_id or hct`")
    q0_patient_Regensburg = q0_patient_Regensburg[!is.na(q0_patient_Regensburg$patient_id),];
    #Some patients have duplicated field. Select only the first HCT case.
    q0_patient_Regensburg = q0_patient_Regensburg[order(q0_patient_Regensburg$patient_id, q0_patient_Regensburg$hct),];
    q0_patient_Regensburg = q0_patient_Regensburg[!duplicated(q0_patient_Regensburg$patient_id),];
    
    q0_patient_Regensburg$last_contact_clean = q0_patient_Regensburg$last_follow_up; #I update the database. last_contact_clean field no longer exists, but it is unified in last_follow_up;
    
    if(0){
      #I updated a clean patients table in which I created a unified date_of_death field.
      ind_is_dead =!is.na(q0_patient_Regensburg$date_of_death__updated_july_2017);
      if(any (!is.na(q0_patient_Regensburg$last_contact_clean[ind_is_dead]) ) ){
        warning("last_contact_clean incompatible with date_of_death_updated_july_2017");
      }
      if(any(is.na(q0_patient_Regensburg$cause_of_death__updated_july_2017[ind_is_dead])) |
         any(!is.na(q0_patient_Regensburg$cause_of_death__updated_july_2017[!ind_is_dead])) ){
        stop("Incompatible alive/death cases have/don't have a value in cod");
      }
    }
    ind_is_dead =(q0_patient_Regensburg$vital_status_alive=="No");
    if(any(q0_patient_Regensburg$last_contact_clean[ind_is_dead] != q0_patient_Regensburg$date_of_death[ind_is_dead])) {
      warning("last_contact_clean incompatible with date_of_death");
    }
    if(any(is.na(q0_patient_Regensburg$cause_of_death[ind_is_dead])) ){
      stop("Incompatible alive/death cases have/don't have a value in cod");
    }
    q0_patient_Regensburg$event_day =  q0_patient_Regensburg$date_of_death -
      q0_patient_Regensburg$hct;
    q0_patient_Regensburg$event_day[!ind_is_dead] = q0_patient_Regensburg$last_contact_clean[!ind_is_dead] -
      q0_patient_Regensburg$hct[!ind_is_dead];
    q0_patient_Regensburg$event_type = ind_is_dead; #1 means death, 0 means alive.
    
    q0_patient_Regensburg$event_competing_type = 4; #Initialize with all values equal to 4;
    q0_patient_Regensburg$event_competing_type[!ind_is_dead]=0;
    q0_patient_Regensburg$event_competing_type[grepl("gvhd",q0_patient_Regensburg$cause_of_death,ignore.case = T)] = 3;
    q0_patient_Regensburg$event_competing_type[grepl("graft",q0_patient_Regensburg$cause_of_death,ignore.case = T)] = 2;
    q0_patient_Regensburg$event_competing_type[q0_patient_Regensburg$relapse_pod=="Yes"] = 1;
    
    #Update `event_competing_day` for cases of relapse
    q0_patient_Regensburg$event_competing_day = q0_patient_Regensburg$event_day;
    ind_relapse=q0_patient_Regensburg$event_competing_type==1; 
    warning("Relapse day has to be fixed! Some is not there, others are MONTH/YEAR");
    #q0_patient_Regensburg$relapse__pod_date[q0_patient_Regensburg$relapse__pod_date=="May 2015"] = "2015-05-15";
    #q0_patient_Regensburg$relapse__pod_date[q0_patient_Regensburg$relapse__pod_date=="June 2015"] = "2015-06-15";
    #q0_patient_Regensburg$relapse__pod_date[q0_patient_Regensburg$relapse__pod_date=="March 2015"] = "2015-03-15";
    
    q0_patient_Regensburg$event_competing_day = q0_patient_Regensburg$event_day;
    competing_day_relapse = as.Date(q0_patient_Regensburg$relapse__pod_date[ind_relapse], "%Y-%m-%d") - q0_patient_Regensburg$hct[ind_relapse];
    q0_patient_Regensburg$event_competing_day[ind_relapse] = competing_day_relapse;
    
    d_set_event = data.frame(patient_id=q0_patient_Regensburg$patient_id,
                             event_day=q0_patient_Regensburg$event_day,
                             event_type=q0_patient_Regensburg$event_type,
                             event_competing_day=q0_patient_Regensburg$event_competing_day,
                             event_competing_type=q0_patient_Regensburg$event_competing_type,
                             institution="Regensburg")
    
    
    #Add `factor` label for competing events.
    competing_event_groups = c("censor",
                               "relapse",
                               "no engraftment",
                               "GVHD related mortality",
                               "other");
    competing_event_groups_as_f = factor(competing_event_groups,
                                         levels=unique(competing_event_groups));
    
    d_set_event$event_competing_type_as_f =competing_event_groups_as_f[d_set_event$event_competing_type+1];
    
    q = d_set_event;
    
    
  }
  
  if(query_id==13.13){
    #Colecting data for survival Analysis Duke only!
    
    q0_patient_Duke = get_data_from_query_OTU(0,"patient_duke_ag");
    q0_patient_Duke = q0_patient_Duke[q0_patient_Duke$tp_type=="Allo",];
    cat("select only Allo");
    tp_type_set <- c("Allo", "Auto", "N/A", "not transplanted", "Not transplanted", NA);
    if( !all(q0_patient_Duke$tp_type %in% tp_type_set) ){
      #!all(q0_patient_Regensburg$kind_of_transplantation %in% tp_type_set) ){
      stop("Unexpected Tptype. I was expecting Allo or Auto only");
    }
    
    #Duke
    #q0_patient_Duke = get_data_from_query_OTU(0,"patient_duke_ag");
    #q0_patient_Duke = q0_patient_Duke[q0_patient_Duke$tp_type=="Allo",];
    
    #Some patients have duplicated field. Select only the first HCT case.
    q0_patient_Duke = q0_patient_Duke[order(q0_patient_Duke$pid, q0_patient_Duke$transplant_date),];
    q0_patient_Duke = q0_patient_Duke[!duplicated(q0_patient_Duke$pid),];
    warning("removing 2nd or more Allo transplants! You may want to censor competing risk survival on those points");
    
    q0_patient_Duke$last_contact_inferred = q0_patient_Duke$last_contact;
    q0_patient_Duke$last_contact_inferred[is.na(q0_patient_Duke$last_contact)] = as.Date("08/29/2017","%m/%d/%Y");
    q0_patient_Duke$last_contact_inferred[q0_patient_Duke$last_contact_inferred < q0_patient_Duke$transplant_date] = as.Date("08/29/2017","%m/%d/%Y");
    q0_patient_Duke$last_contact[is.na(q0_patient_Duke$last_contact)] = as.Date("08/29/2017","%m/%d/%Y");
    q0_patient_Duke$event_day = q0_patient_Duke$last_contact_inferred - q0_patient_Duke$transplant_date;
    ind_is_dead = q0_patient_Duke$vital_status=="Dead";
    q0_patient_Duke$cause_of_death[is.na(q0_patient_Duke$cause_of_death)] ="";
    q0_patient_Duke$cause_of_death[q0_patient_Duke$cause_of_death=="NA"] ="";
    q0_patient_Duke$event_type = ind_is_dead+0;
    #    q0_patient_Duke = q0_patient_Duke[q0_patient_Duke$tp_type=="Allo",]
    if(any(q0_patient_Duke$cause_of_death=="" & ind_is_dead)  | 
       any(q0_patient_Duke$cause_of_death!="" & !ind_is_dead) ){
       #any(is.na(q0_patient_Regensburg$cause_of_death__updated_july_2017) ) ){
      warning("Fix cod in DUKE");
      warning("Incompatible alive/death cases have/don't have a value in cod");
      warning("cod will be grouped as `others`");
      warning("Fix cod in DUKE");
    }
    #q0_patient_Duke = q0_patient_Duke[!(q0_patient_Duke$cause_of_death=="" & ind_is_dead),];
    q0_patient_Duke$event_competing_type = 4; #Initialize with all values equal to 4;
    ind_is_dead = q0_patient_Duke$vital_status=="Dead";
    q0_patient_Duke$event_competing_type[!ind_is_dead]=0; #Alive cases
    q0_patient_Duke$event_competing_type[grepl("gvhd",q0_patient_Duke$cause_of_death,ignore.case = T)] = 3;
    q0_patient_Duke$event_competing_type[grepl("graft",q0_patient_Duke$cause_of_death,ignore.case = T)] = 2;
    ind_relapse = !is.na(q0_patient_Duke$transplant_relapse) |
      !is.na(q0_patient_Duke$progressive_disease_date) |
      !is.na(q0_patient_Duke$dbtc_relapse_date) |
      grepl("Recurrent",q0_patient_Duke$cause_of_death,ignore.case = T);
    q0_patient_Duke$event_competing_type[ind_relapse] = 1;
    
    q0_patient_Duke$event_competing_day = q0_patient_Duke$event_day;
    
    
    pre_competing_day_relapse1 = q0_patient_Duke$transplant_relapse[ind_relapse] - q0_patient_Duke$transplant_date[ind_relapse];
    pre_competing_day_relapse2 = q0_patient_Duke$dbtc_relapse_date[ind_relapse] - q0_patient_Duke$transplant_date[ind_relapse];
    pre_competing_day_relapse3 = q0_patient_Duke$progressive_disease_date[ind_relapse] -q0_patient_Duke$transplant_date[ind_relapse];
    competing_day_relapse = apply(cbind(pre_competing_day_relapse1,
                                        pre_competing_day_relapse2,
                                        pre_competing_day_relapse3),
                                  1, 
                                  function(x) {min(x, na.rm = T)})
    
    q0_patient_Duke$event_competing_day[ind_relapse] = competing_day_relapse;
    
    d_set_event = data.frame(patient_id=q0_patient_Duke$pid,
                             event_day=q0_patient_Duke$event_day,
                             event_type=q0_patient_Duke$event_type,
                             event_competing_type=q0_patient_Duke$event_competing_type,
                             event_competing_day=q0_patient_Duke$event_competing_day,
                             institution="Duke");
    
    #Add `factor` label for competing events.
    competing_event_groups = c("censor",
                               "relapse",
                               "no engraftment",
                               "GVHD related mortality",
                               "other");
    competing_event_groups_as_f = factor(competing_event_groups,
                                         levels=unique(competing_event_groups));
    
    d_set_event$event_competing_type_as_f =competing_event_groups_as_f[d_set_event$event_competing_type+1];
    
    q = d_set_event;
    
  }
  
  if(query_id==13.131){
    #Colecting data for survival Analysis Duke only! (Using `patient_duke_simplified_ag`)
    
    q0_patient_Duke = get_data_from_query_OTU(0,"patient_duke_simplified_ag");
    q0_patient_Duke = q0_patient_Duke[q0_patient_Duke$tp_type=="Allo",];
    cat("select only Allo");
    tp_type_set <- c("Allo", "Auto", "N/A", "not transplanted", "Not transplanted", NA);
    if( !all(q0_patient_Duke$tp_type %in% tp_type_set) ){
      #!all(q0_patient_Regensburg$kind_of_transplantation %in% tp_type_set) ){
      stop("Unexpected Tptype. I was expecting Allo or Auto only");
    }
    
    #Duke
    #q0_patient_Duke = get_data_from_query_OTU(0,"patient_duke_ag");
    #q0_patient_Duke = q0_patient_Duke[q0_patient_Duke$tp_type=="Allo",];
    
    #Some patients have duplicated field. Select only the first HCT case.
    q0_patient_Duke = q0_patient_Duke[order(q0_patient_Duke$pid, q0_patient_Duke$transplant_date),];
    q0_patient_Duke = q0_patient_Duke[!duplicated(q0_patient_Duke$pid),];
    warning("removing 2nd or more Allo transplants! You may want to censor competing risk survival on those points");
    
    q0_patient_Duke$last_contact_inferred = q0_patient_Duke$last_contact;
    q0_patient_Duke$last_contact_inferred[is.na(q0_patient_Duke$last_contact)] = as.Date("08/29/2017","%m/%d/%Y");
    q0_patient_Duke$last_contact_inferred[q0_patient_Duke$last_contact_inferred < q0_patient_Duke$transplant_date] = as.Date("08/29/2017","%m/%d/%Y");
    q0_patient_Duke$last_contact[is.na(q0_patient_Duke$last_contact)] = as.Date("08/29/2017","%m/%d/%Y");
    q0_patient_Duke$event_day = q0_patient_Duke$last_contact_inferred - q0_patient_Duke$transplant_date;
    ind_is_dead = q0_patient_Duke$vital_status=="Dead";
    q0_patient_Duke$cause_of_death[is.na(q0_patient_Duke$cause_of_death)] ="";
    q0_patient_Duke$cause_of_death[q0_patient_Duke$cause_of_death=="NA"] ="";
    q0_patient_Duke$event_type = ind_is_dead+0;
    d_set_event = data.frame(patient_id=q0_patient_Duke$pid,
                             event_day=q0_patient_Duke$event_day,
                             event_type=q0_patient_Duke$event_type,
                             institution="Duke");
    q = d_set_event;
    
  }
  
  if(query_id==13.14){
    #Colecting data for survival Analysis Duke only!
    
    q0_patient_hokkaido = get_data_from_query_OTU(0,"patient_hokkaido_ag");
    
    #Hokkaido
    #q0_patient_Duke = get_data_from_query_OTU(0,"patient_duke_ag");
    #q0_patient_Duke = q0_patient_Duke[q0_patient_Duke$tp_type=="Allo",];
    
    #Some patients have duplicated field. Select only the first HCT case.
    q0_patient_hokkaido$transplant_date = q0_patient_hokkaido$hct;
    q0_patient_hokkaido = q0_patient_hokkaido[order(q0_patient_hokkaido$patientid, q0_patient_hokkaido$transplant_date),]
    q0_patient_hokkaido = q0_patient_hokkaido[!duplicated(q0_patient_hokkaido$patientid),];
    
    q0_patient_hokkaido$last_contact_inferred = q0_patient_hokkaido$last_contact_date;
    #q0_patient_hokkaido$last_contact_inferred[is.na(q0_patient_hokkaido$last_contact)] = as.Date("08/29/2017","%m/%d/%Y");
    #q0_patient_hokkaido$last_contact_inferred[q0_patient_hokkaido$last_contact_inferred < q0_patient_hokkaido$transplant_date] = as.Date("08/29/2017","%m/%d/%Y");
    #q0_patient_hokkaido$last_contact[is.na(q0_patient_hokkaido$last_contact)] = as.Date("08/29/2017","%m/%d/%Y");
    q0_patient_hokkaido$event_day = q0_patient_hokkaido$last_contact_inferred - q0_patient_hokkaido$transplant_date;
    ind_is_dead = q0_patient_hokkaido$vital_status=="dead";
    q0_patient_hokkaido$cause_of_death = q0_patient_hokkaido$cod_copeland;
    q0_patient_hokkaido$cause_of_death[is.na(q0_patient_hokkaido$cause_of_death)] ="";
    q0_patient_hokkaido$cause_of_death[q0_patient_hokkaido$cause_of_death=="NA"] ="";
    q0_patient_hokkaido$event_type = ind_is_dead+0;
    #    q0_patient_hokkaido = q0_patient_hokkaido[q0_patient_hokkaido$tp_type=="Allo",]
    if(any(q0_patient_hokkaido$cause_of_death=="" & ind_is_dead)  | 
       any(q0_patient_hokkaido$cause_of_death!="" & !ind_is_dead) ){
      #any(is.na(q0_patient_Regensburg$cause_of_death__updated_july_2017) ) ){
      warning("Fix cod in Hokkaido");
      #warning("Incompatible alive/death cases have/don't have a value in cod");
      #warning("cod will be grouped as `others`");
      #warning("Fix cod in DUKE");
    }
    #q0_patient_hokkaido = q0_patient_hokkaido[!(q0_patient_hokkaido$cause_of_death=="" & ind_is_dead),];
    q0_patient_hokkaido$event_competing_type = 4; #Initialize with all values equal to 4;
    ind_is_dead = q0_patient_hokkaido$vital_status=="dead";
    q0_patient_hokkaido$event_competing_type[!ind_is_dead]=0; #Alive cases
    q0_patient_hokkaido$event_competing_type[grepl("gvhd",q0_patient_hokkaido$cause_of_death,ignore.case = T)] = 3;
    q0_patient_hokkaido$event_competing_type[grepl("graft",q0_patient_hokkaido$cause_of_death,ignore.case = T)] = 2;
    #q0_patient_hokkaido$relapse_date[!grepl("relapse",q0_patient_hokkaido$cause_of_death,ignore.case = T)];
    q0_patient_hokkaido$transplant_relapse = q0_patient_hokkaido$relapse_date;
    q0_patient_hokkaido$progressive_disease_date = q0_patient_hokkaido$pod_date;
    ind_relapse = !is.na(q0_patient_hokkaido$transplant_relapse) |
      !is.na(q0_patient_hokkaido$progressive_disease_date);
      #!is.na(q0_patient_hokkaido$relapse_pod) | #This is removed! Hokkaido group uses value "n" when relapse_pod does not occur instead of "NA".
      #grepl("Recurrent",q0_patient_hokkaido$cause_of_death,ignore.case = T);
    q0_patient_hokkaido$event_competing_type[ind_relapse] = 1;
    
    q0_patient_hokkaido$event_competing_day = q0_patient_hokkaido$event_day;
    
    
    pre_competing_day_relapse1 = q0_patient_hokkaido$transplant_relapse[ind_relapse] - q0_patient_hokkaido$transplant_date[ind_relapse];
    #pre_competing_day_relapse2 = q0_patient_hokkaido$dbtc_relapse_date[ind_relapse] - q0_patient_hokkaido$transplant_date[ind_relapse];
    pre_competing_day_relapse3 = q0_patient_hokkaido$progressive_disease_date[ind_relapse] -q0_patient_hokkaido$transplant_date[ind_relapse];
    competing_day_relapse = apply(cbind(pre_competing_day_relapse1,
    #                                    pre_competing_day_relapse2,
                                        pre_competing_day_relapse3),
                                  1, 
                                  function(x) {min(x, na.rm = T)})
    
    q0_patient_hokkaido$event_competing_day[ind_relapse] = competing_day_relapse;
    
    d_set_event = data.frame(patient_id=q0_patient_hokkaido$patientid,
                             event_day=q0_patient_hokkaido$event_day,
                             event_type=q0_patient_hokkaido$event_type,
                             event_competing_type=q0_patient_hokkaido$event_competing_type,
                             event_competing_day=q0_patient_hokkaido$event_competing_day,
                             institution="Hokkaido");
    
    #Add `factor` label for competing events.
    competing_event_groups = c("censor",
                               "relapse",
                               "no engraftment",
                               "GVHD related mortality",
                               "other");
    competing_event_groups_as_f = factor(competing_event_groups,
                                         levels=unique(competing_event_groups));
    
    d_set_event$event_competing_type_as_f =competing_event_groups_as_f[d_set_event$event_competing_type+1];
    
    q = d_set_event;
    
  }
  
  if(query_id==13.2){
    q13p21 = get_data_from_query_OTU(13.21);
    q13p22 = get_data_from_query_OTU(13.22);
    q13p23 = get_data_from_query_OTU(13.23);
    q13p24 = get_data_from_query_OTU(13.24);
    
    q = rbind(q13p21,
              q13p22,
              q13p23,
              q13p24);
    
  }
  if(query_id== -13.2){
    stop("obsolete (as indicated by the negative value)! I am using subqueries 13.21... 13.24 to make it.")
    #Collecting data for survival analysis
    #13.2: Linking patient_id, sampleid and day_relative_to_hsct.
    
    q11p1 = get_data_from_query_OTU(11.1); #Get castori (MSK) data.
    q0_samples_regensburg = get_data_from_query_OTU(0, "samples_regensburg_ag");
    #q0_samples_duke  = get_data_from_query_OTU(0, "samples_duke_ag");
    qm1_dukecleaning = get_data_from_query_OTU(-1,"select sd.sampleid, sd.patient_uid, collection_date, sd.timepoint,  transplant_date,  tp_type, collection_date-transplant_date as day_relative_to_hct from samples_duke_clean_ag as sd join patient_duke_ag pd on sd.patient_uid=pd.pid  order by sd.sampleid, pd.pid, transplant_date");
    qm1_dukecleaning = qm1_dukecleaning[qm1_dukecleaning$tp_type=="Allo",];
    warning("Considering only Allo cases to map sampleid and transplant date for Duke patients");
    qm1_dukecleaning= qm1_dukecleaning[order(qm1_dukecleaning$patient_uid,qm1_dukecleaning$transplant_date),];
    qm1_dukecleaning = qm1_dukecleaning[!duplicated(qm1_dukecleaning$sampleid),];
    qm1_dukecleaning = qm1_dukecleaning[!qm1_dukecleaning$transplant_date=="5555-05-05",]; #Remove the cases without transplant information.
    #qm1_dukecleaning[order(-qm1_dukecleaning$day_relative_to_hct),][1:5,]
    
    qm1_dukecleaning$day_relative_to_hsct = qm1_dukecleaning$collection_date-qm1_dukecleaning$transplant_date;
    
    d_set_sample = rbind(data.frame(patient_id = q11p1$mrn,
                                    sampleid = q11p1$sampleid,
                                    day_relative_to_hsct = q11p1$day_relative_to_hsct,
                                    institution="MSK_allo"),
                         data.frame(patient_id=q0_samples_regensburg$patient_uid,
                                    sampleid=q0_samples_regensburg$sampleid,
                                    day_relative_to_hsct=q0_samples_regensburg$timepoint_as_integer,
                                    institution="Regensburg"),
                         #data.frame(patient_id=q0_samples_duke$patient_uid,
                         #            sampleid=q0_samples_duke$sampleid,
                         #           day_relative_to_hsct=q0_samples_duke$timepoint_as_integer,
                         #          institution="Duke"));
                         data.frame(patient_id=qm1_dukecleaning$patient_uid,
                                    sampleid=qm1_dukecleaning$sampleid,
                                    day_relative_to_hsct=qm1_dukecleaning$day_relative_to_hsct,
                                    institution="Duke"));
    
    q = d_set_sample;
  }
  
  if(query_id==13.21){
    #Collecting data for survival analysis for MSK only.
    #13.21: Linking patient_id, sampleid and day_relative_to_hsct
    
    q11p1 = get_data_from_query_OTU(11.1); #Get castori (MSK) data.

    d_set_sample = data.frame(patient_id = q11p1$mrn,
                              sampleid = q11p1$sampleid,
                              day_relative_to_hsct = q11p1$day_relative_to_hsct,
                              institution="MSK_allo");
    
    q = d_set_sample;
  }
  
  if(query_id==13.22){
    #Collecting data for survival analysis for Regensburg only.
    #13.21: Linking patient_id, sampleid and day_relative_to_hsct
    
    q0_samples_regensburg = get_data_from_query_OTU(0, "samples_regensburg_ag");
    
    d_set_sample = data.frame(patient_id=q0_samples_regensburg$patient_uid,
                              sampleid=q0_samples_regensburg$sampleid,
                              day_relative_to_hsct=q0_samples_regensburg$timepoint_as_integer,
                              institution="Regensburg");
    
    q = d_set_sample;
    
  }
  
  if(query_id==13.23){
    #Collecting data for survival analysis for Duke only.
    #13.23: 
    
    qm1_dukecleaning = get_data_from_query_OTU(-1,"select sd.sampleid, sd.patient_uid, collection_date, sd.timepoint,  transplant_date,  tp_type, collection_date-transplant_date as day_relative_to_hct from samples_duke_clean_ag as sd join patient_duke_ag pd on sd.patient_uid=pd.pid  order by sd.sampleid, pd.pid, transplant_date");
    qm1_dukecleaning = qm1_dukecleaning[qm1_dukecleaning$tp_type=="Allo",];
    warning("Considering only Allo cases to map sampleid and transplant date for Duke patients");
    qm1_dukecleaning= qm1_dukecleaning[order(qm1_dukecleaning$patient_uid,qm1_dukecleaning$transplant_date),];
    qm1_dukecleaning = qm1_dukecleaning[!duplicated(qm1_dukecleaning$sampleid),];
    qm1_dukecleaning = qm1_dukecleaning[!qm1_dukecleaning$transplant_date=="5555-05-05",]; #Remove the cases without transplant information.
    #qm1_dukecleaning[order(-qm1_dukecleaning$day_relative_to_hct),][1:5,]
    
    d_set_sample = data.frame(patient_id=qm1_dukecleaning$patient_uid,
                              sampleid=qm1_dukecleaning$sampleid,
                              day_relative_to_hsct=qm1_dukecleaning$day_relative_to_hct,
                              institution="Duke");
    
    q = d_set_sample;

  }
  
  if(query_id==13.231){
    #Collecting data for survival analysis for Duke only. Consider only `patient_duke_simplified_ag` 
    #13.23: 
    
    qm1_dukecleaning = get_data_from_query_OTU(-1, "select sd.sampleid, sd.patient_uid, collection_date, sd.timepoint,  transplant_date,  tp_type, collection_date-transplant_date as day_relative_to_hct from samples_duke_clean_ag as sd join patient_duke_simplified_ag pd on sd.patient_uid=pd.pid  order by sd.sampleid, pd.pid, transplant_date");
    qm1_dukecleaning = qm1_dukecleaning[qm1_dukecleaning$tp_type=="Allo",];
    warning("Considering only Allo cases to map sampleid and transplant date for Duke patients");
    qm1_dukecleaning= qm1_dukecleaning[order(qm1_dukecleaning$patient_uid,qm1_dukecleaning$transplant_date),];
    qm1_dukecleaning = qm1_dukecleaning[!duplicated(qm1_dukecleaning$sampleid),];
    qm1_dukecleaning = qm1_dukecleaning[!qm1_dukecleaning$transplant_date=="5555-05-05",]; #Remove the cases without transplant information.
    #qm1_dukecleaning[order(-qm1_dukecleaning$day_relative_to_hct),][1:5,]
    
    d_set_sample = data.frame(patient_id=qm1_dukecleaning$patient_uid,
                              sampleid=qm1_dukecleaning$sampleid,
                              day_relative_to_hsct=qm1_dukecleaning$day_relative_to_hct,
                              institution="Duke");
    
    q = d_set_sample;
    
  }
  
  if(query_id==13.24){
    #Collecting data for survival analysis for Hokkaido only.
    #13.24: 
    
    q0_samples_hokkaido = get_data_from_query_OTU(0, "samples_hokkaido_ag");
    
    d_set_sample = data.frame(patient_id=q0_samples_hokkaido$patientid,
                              sampleid=q0_samples_hokkaido$sampleid,
                              day_relative_to_hsct=q0_samples_hokkaido$day_sample_collection,
                              institution="Hokkaido");
    
    q = d_set_sample;
    
  }
  if(query!=''){
    
    q = dbGetQuery(con, query);
  }
  
  
  #####
  #This part makes special changes to output!
  #####
  if(query_id==2){
    q = q[order(q$key),];
  }
  if(query_id=="2_str"){
    q = q[order(q$key),];
  }
  
  if(query_id==3){
    q$count_relative = q$count/q$count_total;
  }
  if(query_id==3.1){
    warning("I made this query to get total_count and total_enterococcus per sample. I may do this is in a better way. As it is, key and uploaded_date gets duplicated.");
    q$genus=as.factor(q$genus);
    q$sampleid=as.factor(q$sampleid);
    a_total_per_sampleid= aggregate( data.frame(total_count=q$count), by=data.frame(sampleid=q$sampleid),sum);
    a_total_per_sampleid_enterococcus= aggregate( data.frame(total_enterococcus=q$count[q$genus=="Enterococcus"]), by=data.frame(sampleid=q$sampleid[q$genus=="Enterococcus"]),sum);
    m1 = merge(q, a_total_per_sampleid);
    m2 = merge(m1, a_total_per_sampleid_enterococcus, all.x=T);
    q=m2;
  }
  
  if(query_id==12){
    #Some sampleids have "..pool###", this part creates a `clean` version of it.
    q$sampleid_clean = q$sampleid; #Some sampleids have "..pool###" 
    q$sampleid_clean = sapply(strsplit(q$sampleid,"\\.\\."),
                                function(x) {return(x[[1]][1])});
    
  }
  
  return(data.table(q));
  
}

#Compare alpha diversity from count with what is output by pipeline.
#q3 = get_data_from_query_OTU(3);
#q3 = data.table(q3);
#a=q3[, .(simpsr = 1/(sum ( (count/count_total)^2))), by="sampleid"];
#m = merge(q12, a, by.x="sampleid", by.y="sampleid")
