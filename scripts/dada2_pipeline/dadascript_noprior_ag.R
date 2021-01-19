#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(dada2); packageVersion("dada2")
library(ShortRead); #To subsample fastq sequences.
#forward
path <- args[1];
filter_trunclen = c(180,180); #c(180,180) is default
asv_len=seq(300,360); #c(300,360) is default

# filter_trunclen = c(145,151); #c(180,180) is default
# asv_len=seq(240,260); #c(300,360) is default


#fnFs <- Sys.glob(file.path(path, "*R1*.fastq.gz"))
fnFs <- Sys.glob(file.path(path, "*_R1*.fastq.gz"))
fnFs <- fnFs[lapply(fnFs,function(x) length(grep("_unsigned",x,value=FALSE))) == 0]

#backward
#fnRs <- Sys.glob(file.path(path, "*R2*.fastq.gz"))
fnRs <- Sys.glob(file.path(path, "*_R2*.fastq.gz"))
#remove unsigned
fnRs <- fnRs[lapply(fnRs,function(x) length(grep("_unsigned",x,value=FALSE))) == 0]

print(fnRs);

forwardfileinfo <- file.info(fnFs)
forwardsizes <- forwardfileinfo[,"size"]
reversefileinfo <- file.info(fnRs)
reversesizes <- reversefileinfo[,"size"]#forwardfileinfo[,"size"]

fnFs <- fnFs[which(forwardsizes>10000 | reversesizes>10000)]
fnRs <- fnRs[which(forwardsizes>10000 | reversesizes>10000)]
#get sample names (update depending on filename conventions)
#sample.names <- sapply(strsplit(sapply(strsplit(basename(fnFs), ".fastq.gz"), `[`, 1), "fastq_"), `[`,2)
#sample.names <- sapply(strsplit(sapply(strsplit(basename(fnFs), ".fastq"), `[`, 1), "fastq_"), `[`,2)
sample.names = sapply(strsplit(basename(fnFs), "_R1.fastq.gz"), `[`, 1);


#temporary files for storing filtered data
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

multithread_var=8;#F;#T; #F; #F;#T; #FALSE;
randomize_var=F; #T;

temp_log_file = sprintf("%s/temp_log.txt",path);
write.table(sprintf("%s - before filterAndTrim Fs",Sys.time()), temp_log_file);

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=filter_trunclen,
                     maxN=0, maxEE=2, truncQ=2,
                     compress=TRUE, multithread=multithread_var) #TRUE)

#Capping sampleids up to 50k reads.
write.table(sprintf("%s - before capping",Sys.time()), temp_log_file);

Nread_cap=1e5;
#Nread_cap=1e10;

file.remove(dir( dirname(filtFs[1]),pattern="temp",full.names=T ));
do_capping=T;#F;

if(do_capping){
  
  for(i in 1:length(filtFs) ){
    R1_fq_file=filtFs[i];
    R2_fq_file=filtRs[i];
    
    sampler_R1 <- FastqSampler(R1_fq_file,Nread_cap);
    sampler_R2 <- FastqSampler(R2_fq_file,Nread_cap);
    
    set.seed(123); sampler_R1_y <- yield(sampler_R1);
    set.seed(123); sampler_R2_y <- yield(sampler_R2);	
    
    temp_R1_fq_file=gsub("fastq.gz","fastq.gz.temp",R1_fq_file);
    temp_R2_fq_file=gsub("fastq.gz","fastq.gz.temp",R2_fq_file);
    
    writeFastq(sampler_R1_y, temp_R1_fq_file);
    writeFastq(sampler_R2_y, temp_R2_fq_file);
    
    file.rename(temp_R1_fq_file, R1_fq_file);
    file.rename(temp_R2_fq_file, R2_fq_file);
    
    close(sampler_R1);
    close(sampler_R2);
    
  }
  
}

write.table(sprintf("%s - before learnErrors Fs",Sys.time()), temp_log_file, append=T);
errF <- learnErrors(filtFs, multithread=multithread_var, randomize=randomize_var);#imultithread=multithread_var)

write.table(sprintf("%s - before learnErrors Rs",Sys.time()), temp_log_file, append=T);
errR <- learnErrors(filtRs, multithread=multithread_var, randomize=randomize_var);#TRUE)

write.table(sprintf("%s - before derep F",Sys.time()), temp_log_file, append=T);
derepFs <- derepFastq(filtFs);#, verbose=TRUE)

write.table(sprintf("%s - before derep R",Sys.time()), temp_log_file, append=T);
derepRs <- derepFastq(filtRs);#, verbose=TRUE)

#names(derepFs) <- sample.names
#names(derepRs) <- sample.names

write.table(sprintf("%s - before dada F",Sys.time()), temp_log_file, append=T);
dadaFs <- dada(derepFs, err=errF,multithread=multithread_var);#TRUE)

write.table(sprintf("%s - before dada R",Sys.time()), temp_log_file, append=T);
dadaRs <- dada(derepRs, err=errR,multithread=multithread_var);#TRUE)

write.table(sprintf("%s - before merge",Sys.time()), temp_log_file, append=T);
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs);#, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)

write.table(sprintf("%s - before removebimera",Sys.time()), temp_log_file, append=T);

seqtab_prebimera_file = sprintf("%s/seqtab_prebimera.rds",path);
saveRDS(seqtab, seqtab_prebimera_file);
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=multithread_var)#, verbose=TRUE)

seqtabfinal <- seqtab.nochim[,nchar(colnames(seqtab.nochim)) %in% asv_len, drop=F]
#rownames(seqtabfinal) = names(derepFs);
rownames(seqtabfinal) = sample.names; # I changed names to sample.names to avoid error when dealing with 1 dimension data.

outfile <- file.path(path,'seqtab_noprior.csv');
write.csv(seqtabfinal,outfile);

outfile <- file.path(path, 'seqtab_noprior.rds')
saveRDS(seqtabfinal, file=outfile)

ASV_sequences <- colnames(seqtabfinal);

fasta_file= file.path(path,'ASV_sequences.fasta');
library(seqinr);
write.fasta(sequences=as.list(ASV_sequences), names=sprintf("ASV_temp_%d",1:length(ASV_sequences)), fasta_file);

seqtabfinal_counts_raw = as.data.frame(seqtabfinal); #Name table column as `ASV_temp_#` instead of full ASV sequence.
colnames(seqtabfinal_counts_raw) = sprintf("ASV_temp_%d",1:ncol(seqtabfinal_counts_raw));
library(reshape2);
seqtabfinal_counts_raw$oligos_id = rownames(seqtabfinal);
seqtabfinal_counts = melt(seqtabfinal_counts_raw, id.vars="oligos_id", variable.name="asv_temp_id",value.name="count" )
seqtabfinal_counts = seqtabfinal_counts[seqtabfinal_counts$count!=0,];

count_table_outfile <- file.path(path, "asv_counts.csv");
write.csv(seqtabfinal_counts, count_table_outfile, row.names = F);

outfile <- file.path(path, 'errR.rds')
saveRDS(errR, file=outfile)

outfile <- file.path(path, 'errF.rds')
saveRDS(errF, file=outfile)

outfile <- file.path(path, 'derepFs.rds')
saveRDS(derepFs, file=outfile)

outfile <- file.path(path, 'derepRs.rds')
saveRDS(derepRs, file=outfile)

write.table(sprintf("%s - ending script",Sys.time()), temp_log_file, append=T);
