rm(list=ls())
library(survival)
setwd("/Volumes/vandenBrinkLab/Sarah/Project_auto/figures/Scripts for figures/For survival check/")

#### NOTES
#### make sure to check date format

#ptdata2 <- read.csv("msk duke pttable 09 05 19.csv")
#ptdata2 <- read.csv("msk duke pttable 01 18 20.csv")
ptdata2 <- read.csv("../msk duke pttable 03 02 20.csv")
any(is.na(ptdata2$mrn))
any(duplicated(ptdata2$mrn))


ptdata2$ISS.Stage <- NULL

###################### adding in the clean myeloma cohort
#cleanmyeloma <- read.csv("msk duke cohort 1 pts 09 05.19.csv")
cleanmyeloma <- read.csv("../clean myeloma pttable 02 29 20.csv")

ptdata2$CleanMyeloma <- rep(0, dim(ptdata2)[1])
ptdata2$CleanMyeloma[ptdata2$mrn %in% cleanmyeloma$mrn] <- 1
table(ptdata2$CleanMyeloma)

cleanmyeloma$ISS.Stage
ISS.StageUpdated <- cleanmyeloma$ISS.Stage
ISS.StageUpdated[cleanmyeloma$ISS.Stage ==4] <-  NA  ################### ASK NILU


MelDoseNewFile <- cleanmyeloma$Dose.Level
table(cleanmyeloma$Dose.Level, exclude=NULL)


Morethan1Line <- 1*(cleanmyeloma$Lines.Pre.HCT >1)


myelomaISSmerge <- data.frame(mrn=cleanmyeloma$mrn, ISSUpdated=ISS.StageUpdated, MelDoseNewFile=MelDoseNewFile,Morethan1Line=Morethan1Line )

dim(ptdata2)
dim(myelomaISSmerge)


ptdata <- merge(ptdata2, myelomaISSmerge,  all.x=TRUE)
dim(ptdata)

table(ptdata$MelDoseNewFile[ptdata2$CleanMyeloma==1   ], ptdata$Dose.Level[ptdata2$CleanMyeloma==1   ], exclude=NULL)

#####################
HCT <- as.Date(ptdata$HCT, "%m/%d/%Y")
Relapse.POD.Date <- as.Date(ptdata$Relapse.POD.Date, "%m/%d/%Y")
LastFU <-  as.Date(ptdata$Last.Contact.DOD, "%m/%d/%Y")
table(Relapse.POD.Date <=  LastFU)

OSEvent <- 1*(ptdata$Vital.Status == "Dead")
PFSEvent <- OSEvent
PFSEvent[!is.na(Relapse.POD.Date)] <- 1
table(PFSEvent, OSEvent)

OSmonths <- as.numeric(LastFU - HCT)/(365/12)
PFSmonths <- as.numeric(pmin(LastFU, Relapse.POD.Date, na.rm=TRUE) - HCT)/(365/12)
table(PFSmonths < OSmonths, !is.na(Relapse.POD.Date), exclude=NULL)

ptdata$OSmonths <- OSmonths
ptdata$OSEvent <- OSEvent

ptdata$PFSmonths <- PFSmonths
ptdata$PFSEvent <- PFSEvent


PFSlandmark1 <- ptdata$PFSmonths- 17/(365/12)
PFSlandmark1[PFSlandmark1 <=0] <- NA
 
PFSEventlandmark1 <- PFSEvent

tt <- 24
PFSEventlandmark1[PFSlandmark1 > tt & !is.na(PFSlandmark1)] <- 0 
PFSlandmark1[PFSlandmark1 > tt & !is.na(PFSlandmark1)] <- tt

ptdata$PFSlandmark1 <- PFSlandmark1
ptdata$PFSEventlandmark1 <- PFSEventlandmark1

OSlandmark1 <- ptdata$OSmonths -   17/(365/12)
OSlandmark1[OSlandmark1 <=0] <- NA

OSEventlandmark1 <- OSEvent

tt <- 24
OSEventlandmark1[OSlandmark1 > tt] <- 0 
OSlandmark1[OSlandmark1 > tt] <- tt

ptdata$OSlandmark1 <- OSlandmark1
ptdata$OSEventlandmark1 <- OSEventlandmark1

Melphalan <- 1*(ptdata$Regimen_cat %in% c(3,4))
table(ptdata$Regimen_cat, Melphalan, exclude=NULL)
ptdata$Melphalan <- Melphalan

CR.nCR <- 1*(ptdata$Pt.Status_cat %in% c("CR/near CR", "CR1"))
table(ptdata$Pt.Status_cat, CR.nCR)
ptdata$CR.nCR <- CR.nCR


MyelomaAmyloid <- 1*(ptdata$disease_cat  %in% c(2,3))
ptdata$MyelomaAmyloid <- MyelomaAmyloid

table(ptdata$Disease, ptdata$disease_cat,exclude=NULL )


ptdata$OSmonths <- OSmonths
ptdata$OSEvent <- OSEvent


###########################################
###########################################
###########################################
###########################################   get diversity 
###########################################


ptdata2$TT.ANC[is.na(ptdata2$TT.ANC)] <- 8


#sampledata <- read.csv("msk duke sdtable 09 05 19.csv")
sampledata <- read.csv("msk duke sdtable 03 02 20.csv")

sampledata <- sampledata[sampledata$day.post.txp >= 9 & sampledata$day.post.txp <= 16, ]

divseritybymrn <- NULL
checkday <- NULL
for(ii in 1:dim(ptdata)[1]){
  sampledataH <- sampledata[sampledata$mrn %in%  ptdata$mrn[ii] , ]
#  sampledataH <- sampledataH[sampledataH$day.post.txp >  ptdata2$TT.ANC[ii]  , ]
  
 if(dim(sampledataH)[1] ==0 ){
   divseritybymrn <- rbind(divseritybymrn, c( ptdata$mrn[ii], NA,NA, 0 ))
   checkday <- c(checkday, NA)
 }else{
  
   divseritybymrn <- rbind(divseritybymrn, c( ptdata$mrn[ii], median(sampledataH$simpson_reciprocal),9999,dim(sampledataH)[1] ))
   checkday <- c(checkday, min(sampledataH$day.post.txp))
 }
}

divseritybymrn <- data.frame(divseritybymrn)
colnames(divseritybymrn) <- c("mrn", "simpson", "simpson.hctday", "SampleNumber")



table(!is.na(divseritybymrn$simpson))
#### 240 stool sample available


#divseritybymrnfulllist <- divseritybymrn

dim(divseritybymrn)
dim(ptdata)

combineddata <- merge(ptdata,divseritybymrn)
dim(combineddata)

combineddata <- combineddata[!is.na(combineddata$simpson),]
dim(combineddata)

combineddata$SimpsonMED <- cut(combineddata$simpson, c(-Inf, quantile(combineddata$simpson, c(1/2), na.rm=TRUE) ,  Inf))

###########################################################################
###########################################################################  progression-free survival
###########################################################################
###########################################################################

 
write.csv(data.frame(MRN=combineddata$mrn,MyelomaAmyloid=combineddata$MyelomaAmyloid, CleanMyeloma=combineddata$CleanMyeloma, 
                     Simpson= combineddata$simpson   ), paste("MyelomaAmyloid.Lymphoma.MRNWithSample.",Sys.Date(),".csv", sep=""))


mod1 <- coxph(Surv(PFSlandmark1,PFSEventlandmark1 )~log(simpson) +strata(institution), data=combineddata)
mod2 <- coxph(Surv(PFSlandmark1,PFSEventlandmark1 )~SimpsonMED +strata(institution), data=combineddata)
mod3 <- coxph(Surv(PFSlandmark1,PFSEventlandmark1 )~ MyelomaAmyloid +strata(institution), data=combineddata)
mod4 <- coxph(Surv(PFSlandmark1,PFSEventlandmark1 )~Melphalan +strata(institution), data=combineddata)
mod5 <- coxph(Surv(PFSlandmark1,PFSEventlandmark1 )~CR.nCR +strata(institution), data=combineddata)
mod6 <- coxph(Surv(PFSlandmark1,PFSEventlandmark1 )~Age +strata(institution), data=combineddata)



modH <- mod1
estimates1 <- c(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

modH <- mod2
estimates2 <- c(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

modH <- mod3
estimates3 <- c(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

modH <- mod4
estimates4 <- c(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

modH <- mod5
estimates5 <- c(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

modH <- mod6
estimates6 <- c(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))


write.csv(rbind(estimates1, estimates2, estimates3, estimates4, estimates5,estimates6), paste("uniPFSresults.",Sys.Date(),".csv",sep=""))


mod6 <- coxph(Surv(PFSlandmark1,PFSEventlandmark1 )~log(simpson) +MyelomaAmyloid+CR.nCR+ strata(institution) , data=combineddata)
mod7 <- coxph(Surv(PFSlandmark1,PFSEventlandmark1)~SimpsonMED  +MyelomaAmyloid+CR.nCR+ strata(institution) , data=combineddata)


modH <- mod6
estimates6 <- cbind(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

modH <- mod7
estimates7 <- cbind(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

write.csv(rbind(estimates6, estimates7), paste("MVPFSresults.",Sys.Date(),".csv",sep=""))





#################### just myeloma cohort 

combineddata$MelDoseAbove140 <- 1*(combineddata$MelDoseNewFile > 140)


table(combineddata$Dose.Level[combineddata$CleanMyeloma==1], exclude=NULL)

ISSCAT  <- rep(NA, length(combineddata$ISSUpdated))
ISSCAT[combineddata$ISSUpdated %in%c(1) & !is.na(combineddata$ISSUpdated)] <- "Stage 1"
ISSCAT[combineddata$ISSUpdated %in% c(2,3) & !is.na(combineddata$ISSUpdated)] <- "Stage 2.3"
combineddata$ISSCAT <- ISSCAT

mod1 <- coxph(Surv(PFSlandmark1,PFSEventlandmark1 )~log(simpson) +strata(institution), data=combineddata, subset=(CleanMyeloma ==1 ))
mod2 <- coxph(Surv(PFSlandmark1,PFSEventlandmark1 )~SimpsonMED +strata(institution), data=combineddata, subset=(CleanMyeloma ==1 ))
mod3 <- coxph(Surv(PFSlandmark1,PFSEventlandmark1 )~ISSCAT +strata(institution), data=combineddata, subset=(CleanMyeloma ==1 ))
mod4 <- coxph(Surv(PFSlandmark1,PFSEventlandmark1 )~MelDoseAbove140 +strata(institution), data=combineddata, subset=(CleanMyeloma ==1 ))
mod5 <- coxph(Surv(PFSlandmark1,PFSEventlandmark1 )~CR.nCR +strata(institution), data=combineddata, subset=(CleanMyeloma ==1 ))
mod6 <- coxph(Surv(PFSlandmark1,PFSEventlandmark1 )~Age +strata(institution), data=combineddata, subset=(CleanMyeloma ==1 ))
mod7 <- coxph(Surv(PFSlandmark1,PFSEventlandmark1 )~Morethan1Line +strata(institution), data=combineddata, subset=(CleanMyeloma ==1 ))




combineddata$simpson[combineddata$CleanMyeloma ==1]
combineddata$MelDoseAbove140[combineddata$CleanMyeloma ==1]


modH <- mod1
estimates1 <- c(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

modH <- mod2
estimates2 <- c(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

modH <- mod3
estimates3 <- c(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

modH <- mod4
estimates4 <- c(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

modH <- mod5
estimates5 <- c(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

modH <- mod6
estimates6 <- c(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))


modH <- mod7
estimates7 <- c(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))


write.csv(rbind(estimates1, estimates2, estimates3,estimates4, estimates5,estimates6,estimates7), paste("cleanmyelomauniPFSresults.",Sys.Date(),".csv",sep=""))


write.csv(combineddata[combineddata$CleanMyeloma ==1, c("mrn", "simpson", "SimpsonMED" )], "cleanmyelomawithsamplecheck.012120.csv")

#mod6 <- coxph(Surv(PFSlandmark1,PFSEventlandmark1 )~log(simpson) +MelDoseAbove140+ strata(institution) , data=combineddata, subset=(CleanMyeloma ==1 ))
#mod7 <- coxph(Surv(PFSlandmark1,PFSEventlandmark1)~SimpsonMED  +MelDoseAbove140+ strata(institution) , data=combineddata, subset=(CleanMyeloma ==1 ))

 

mod6 <- coxph(Surv(PFSlandmark1,PFSEventlandmark1 )~log(simpson) +CR.nCR+strata(institution) , data=combineddata, subset=(CleanMyeloma ==1 ))
mod7 <- coxph(Surv(PFSlandmark1,PFSEventlandmark1)~SimpsonMED  +CR.nCR+ strata(institution) , data=combineddata, subset=(CleanMyeloma ==1 ))


modH <- mod6
estimates6 <- cbind(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

modH <- mod7
estimates7 <- cbind(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

write.csv(rbind(estimates6, estimates7), paste("cleanmyelomaMVPFSresultsDiseaseStage.",Sys.Date(),".csv",sep=""))



###########################################################################
###########################################################################  overall survival
###########################################################################
###########################################################################

mod1 <- coxph(Surv(OSlandmark1,OSEventlandmark1 )~log(simpson) +strata(institution), data=combineddata)
mod2 <- coxph(Surv(OSlandmark1,OSEventlandmark1 )~SimpsonMED +strata(institution), data=combineddata)
mod3 <- coxph(Surv(OSlandmark1,OSEventlandmark1 )~ MyelomaAmyloid +strata(institution), data=combineddata)
mod4 <- coxph(Surv(OSlandmark1,OSEventlandmark1 )~Melphalan +strata(institution), data=combineddata)
mod5 <- coxph(Surv(OSlandmark1,OSEventlandmark1 )~CR.nCR +strata(institution), data=combineddata)
mod6 <- coxph(Surv(OSlandmark1,OSEventlandmark1 )~Age +strata(institution), data=combineddata)


modH <- mod1
estimates1 <- c(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

modH <- mod2
estimates2 <- cbind(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

modH <- mod3
estimates3 <- c(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

modH <- mod4
estimates4 <- c(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

modH <- mod5
estimates5 <- c(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

modH <- mod6
estimates6 <- c(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))


write.csv(rbind(estimates1, estimates2, estimates3, estimates4, estimates5,estimates6), paste("uniOSresults",Sys.Date(),".csv",sep=""))


mod6 <- coxph(Surv(OSlandmark1,OSEventlandmark1)~log(simpson) +MyelomaAmyloid+CR.nCR+ strata(institution) , data=combineddata)
mod7 <- coxph(Surv(OSlandmark1,OSEventlandmark1)~SimpsonMED  +MyelomaAmyloid+CR.nCR+ strata(institution) , data=combineddata)


modH <- mod6
estimates6 <- cbind(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

modH <- mod7
estimates7 <- cbind(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

write.csv(rbind(estimates6, estimates7), paste("MVOSresults.",Sys.Date(),".csv",sep=""))



#################### just myeloma cohort 
mod1 <- coxph(Surv(OSlandmark1,OSEventlandmark1 )~log(simpson) +strata(institution), data=combineddata, subset=(CleanMyeloma ==1 ))
mod2 <- coxph(Surv(OSlandmark1,OSEventlandmark1 )~SimpsonMED +strata(institution), data=combineddata, subset=(CleanMyeloma ==1 ))
mod3 <- coxph(Surv(OSlandmark1,OSEventlandmark1 )~ISSCAT +strata(institution), data=combineddata, subset=(CleanMyeloma ==1 ))
mod4 <- coxph(Surv(OSlandmark1,OSEventlandmark1 )~MelDoseAbove140 +strata(institution), data=combineddata, subset=(CleanMyeloma ==1 ))
mod5 <- coxph(Surv(OSlandmark1,OSEventlandmark1 )~CR.nCR +strata(institution), data=combineddata, subset=(CleanMyeloma ==1 ))
mod6 <- coxph(Surv(OSlandmark1,OSEventlandmark1 )~Age +strata(institution), data=combineddata, subset=(CleanMyeloma ==1 ))

modH <- mod1
estimates1 <- c(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

modH <- mod2
estimates2 <- c(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

modH <- mod3
estimates3 <- c(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

modH <- mod4
estimates4 <- c(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

modH <- mod5
estimates5 <- c(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

modH <- mod6
estimates6 <- c(paste(round(summary(modH)$conf.int[,1],2), " (",round(summary(modH)$conf.int[,3],2), "-", round(summary(modH)$conf.int[,4],2), ")", sep=""),round(summary(modH)$coefficients[,5],3))

write.csv(rbind(estimates1, estimates2, estimates3,estimates4, estimates5,estimates6), paste("cleanmyelomauniOSresults.",Sys.Date(),".csv",sep=""))




 #############################################################  PFS
 #############################################################  figures 
 #############################################################
 #############################################################
 #############################################################
 

 pdf(paste("Figure1.PFS.",Sys.Date(),".pdf",sep=""), height=5, width=3*6) 
 
 xmax <- 24
 stepatrisk <- 6
 legendxpos <- -4

par(mfrow=c(1,3))
 
 s2 = survfit(Surv(PFSlandmark1 ,PFSEventlandmark1 )~SimpsonMED, data=combineddata) 
 p1 <-  summary(coxph(Surv(PFSlandmark1,PFSEventlandmark1 )~SimpsonMED +strata(institution), data=combineddata))$coefficients[,5]
 

 lgd=c("Below Median","Above Median")
 

 par(mar=c(8,12,4.1,2.1))
 plot(s2,col=c(2:6), conf.int=F, xlim=c(0,xmax),ylim=c(0,1.05),yaxt="n", xaxs="i",yaxs="i",   bty="n",
      xaxt="n",frame.plot=F,xlab="",ylab="",cex.lab=1.2,cex.main=1.2, mark.time=TRUE,
      main="",cex=.6 , lwd=2)
 axis(side=1,seq(0,xmax,by=2))
 axis(side=2,seq(0,1,by=0.2))
 text(xmax*0.75, 0.05, paste("P-value = ", round(p1, 3), sep=""))
 
 mtext("Combined Cohorts", 3, line=1, cex=1)
 mtext("Months from Day +20 Landmark", 1, line=3.5, cex=1)
 mtext("Progression-free Survival", 2, line=3, cex=1) 
 
 nlgd = length(lgd)
 mtext("No. at risk",side=1,line=5,at=legendxpos,font=2,col=1,cex=1)
 for (i in 1:nlgd){
    mtext(lgd[i],side=1,line=(5+i),at=legendxpos,
          font=4,col=(i+1),cex=0.8)
 } 
 
 
 sout=s2
 maxtlist=xmax/stepatrisk   
 slist=rep(NA,maxtlist)
 
 for (k in 1:nlgd) {
    mtext(sout[k]$n.risk[1],side=1,line=(5+k),at=0,cex=0.8)
    for (i in 1:maxtlist) {
       j=i*stepatrisk;
       slist[i] = sout[k]$n.risk[which(sout[k]$time>=j)[1]]
       slist[i] = ifelse(is.na(slist[i]),0,slist[i])
       mtext(slist[i],side=1,line=(5+k),at=j,cex=0.8)
    }
 }
 
 
 
 
 s2 = survfit(Surv(PFSlandmark1 ,PFSEventlandmark1 )~SimpsonMED, data=combineddata, subset=(CleanMyeloma ==1 )) 
 p1 <-  summary(coxph(Surv(PFSlandmark1,PFSEventlandmark1 )~SimpsonMED +strata(institution), data=combineddata, subset=(CleanMyeloma ==1 )))$coefficients[,5]
 
 lgd=c("Below Median","Above Median")
 
 
 par(mar=c(8,12,4.1,2.1))
 plot(s2,col=c(2:6), conf.int=F, xlim=c(0,xmax),ylim=c(0,1.05), yaxt="n", xaxs="i",yaxs="i",   bty="n",
      xaxt="n",frame.plot=F,xlab="",ylab="",cex.lab=1.2, mark.time=TRUE,
      main="",cex=.6 , lwd=2)
 axis(side=1,seq(0,xmax,by=2))
 axis(side=2,seq(0,1,by=0.2))
 text(xmax*0.75, 0.05, paste("P-value = ", round(p1, 3), sep=""))
 
 
 mtext("Myeloma", 3, line=1, cex=1)
 mtext("Months from Day +20 Landmark", 1, line=3.5, cex=1)
 mtext("Progression-free Survival", 2, line=3, cex=1) 
 
 nlgd = length(lgd)
 mtext("No. at risk",side=1,line=5,at=legendxpos,font=2,col=1,cex=1)
 for (i in 1:nlgd){
    mtext(lgd[i],side=1,line=(5+i),at=legendxpos,
          font=4,col=(i+1),cex=0.8)
 } 
 
 
 sout=s2
 maxtlist=xmax/stepatrisk   
 slist=rep(NA,maxtlist)
 
 for (k in 1:nlgd) {
    mtext(sout[k]$n.risk[1],side=1,line=(5+k),at=0,cex=0.8)
    for (i in 1:maxtlist) {
       j=i*stepatrisk;
       slist[i] = sout[k]$n.risk[which(sout[k]$time>=j)[1]]
       slist[i] = ifelse(is.na(slist[i]),0,slist[i])
       mtext(slist[i],side=1,line=(5+k),at=j,cex=0.8)
    }
 }
 
 
 
 s2 = survfit(Surv(PFSlandmark1 ,PFSEventlandmark1 )~SimpsonMED, data=combineddata, subset=(MyelomaAmyloid ==0 )) 
 p1 <-  summary(coxph(Surv(PFSlandmark1,PFSEventlandmark1 )~SimpsonMED +strata(institution), data=combineddata, subset=(MyelomaAmyloid ==0 )))$coefficients[,5]
 
 lgd=c("Below Median","Above Median")
 
 par(mar=c(8,12,4.1,2.1))
 plot(s2,col=c(2:6), conf.int=F, xlim=c(0,xmax),ylim=c(0,1.05),yaxt="n", xaxs="i",yaxs="i",   bty="n",
      xaxt="n",frame.plot=F,xlab="",ylab="",cex.lab=1.2, mark.time=TRUE,
      main="",cex=.6 , lwd=2)
 axis(side=1,seq(0,xmax,by=2))
 axis(side=2,seq(0,1,by=0.2))
 text(xmax*0.75, 0.05, paste("P-value = ", round(p1, 3), sep=""))
 
 mtext("Lymphoma", 3, line=1, cex=1)
 mtext("Months from Day +20 Landmark", 1, line=3.5, cex=1)
 mtext("Progression-free Survival", 2, line=3, cex=1) 
 
 nlgd = length(lgd)
 mtext("No. at risk",side=1,line=5,at=legendxpos,font=2,col=1,cex=1)
 for (i in 1:nlgd){
    mtext(lgd[i],side=1,line=(5+i),at=legendxpos,
          font=4,col=(i+1),cex=0.8)
 } 
 
 
 sout=s2
 maxtlist=xmax/stepatrisk   
 slist=rep(NA,maxtlist)
 
 for (k in 1:nlgd) {
    mtext(sout[k]$n.risk[1],side=1,line=(5+k),at=0,cex=0.8)
    for (i in 1:maxtlist) {
       j=i*stepatrisk;
       slist[i] = sout[k]$n.risk[which(sout[k]$time>=j)[1]]
       slist[i] = ifelse(is.na(slist[i]),0,slist[i])
       mtext(slist[i],side=1,line=(5+k),at=j,cex=0.8)
    }
 }
 
 
 
 
 dev.off()
 
 

 #############################################################  OS
 #############################################################  figures 
 #############################################################
 #############################################################
 #############################################################
 

 
 pdf( paste("Figure1.OS.",Sys.Date(),".pdf",sep=""), height=6, width=3*7) 
 
 xmax <- 24
 stepatrisk <- 6
 legendxpos <- -4
 
 par(mfrow=c(1,3))
 
 s2 = survfit(Surv(OSlandmark1 ,OSEventlandmark1 )~SimpsonMED, data=combineddata) 
 p1 <-  summary(coxph(Surv(OSlandmark1,OSEventlandmark1 )~SimpsonMED +strata(institution), data=combineddata))$coefficients[,5]
 
 
 lgd=c("Below Median","Above Median")
 
 
 par(mar=c(8,12,4.1,2.1))
 plot(s2,col=c(2:6), conf.int=F, xlim=c(0,xmax),ylim=c(0,1.05),yaxt="n", xaxs="i",yaxs="i",   bty="n",
      xaxt="n",frame.plot=F,xlab=" ",ylab="",cex.lab=1.2, mark.time=TRUE,
      main="",cex=.6 , lwd=2)
 axis(side=1,seq(0,xmax,by=2))
 axis(side=2,seq(0,1,by=0.2))
 text(xmax*0.75, 0.05, paste("P-value = ", round(p1, 3), sep=""))
 
 mtext("Combined Cohorts", 3, line=1, cex=1)
 mtext("Months from Day +20 Landmark", 1, line=3.5, cex=1)
 mtext("Overall Survival", 2, line=3, cex=1) 
 
 nlgd = length(lgd)
 mtext("No. at risk",side=1,line=5,at=legendxpos,font=2,col=1,cex=1)
 for (i in 1:nlgd){
   mtext(lgd[i],side=1,line=(5+i),at=legendxpos,
         font=4,col=(i+1),cex=0.8)
 } 
 
 
 sout=s2
 maxtlist=xmax/stepatrisk   
 slist=rep(NA,maxtlist)
 
 for (k in 1:nlgd) {
   mtext(sout[k]$n.risk[1],side=1,line=(5+k),at=0,cex=0.8)
   for (i in 1:maxtlist) {
     j=i*stepatrisk;
     slist[i] = sout[k]$n.risk[which(sout[k]$time>=j)[1]]
     slist[i] = ifelse(is.na(slist[i]),0,slist[i])
     mtext(slist[i],side=1,line=(5+k),at=j,cex=0.8)
   }
 }
 
 
 
 
 s2 = survfit(Surv(OSlandmark1 ,OSEventlandmark1 )~SimpsonMED, data=combineddata, subset=(CleanMyeloma ==1 )) 
 p1 <-  summary(coxph(Surv(OSlandmark1,OSEventlandmark1 )~SimpsonMED +strata(institution), data=combineddata, subset=(CleanMyeloma ==1 )))$coefficients[,5]
 
 lgd=c("Below Median","Above Median")
 
 
 par(mar=c(8,12,4.1,2.1))
 plot(s2,col=c(2:6), conf.int=F, xlim=c(0,xmax),ylim=c(0,1.05),yaxt="n", xaxs="i",yaxs="i",   bty="n",
      xaxt="n",frame.plot=F,xlab="",ylab="",cex.lab=1.2, mark.time=TRUE,
      main="",cex=.6 , lwd=2)
 axis(side=1,seq(0,xmax,by=2))
 axis(side=2,seq(0,1,by=0.2))
 text(xmax*0.75, 0.05, paste("P-value = ", round(p1, 3), sep=""))
 
 
 mtext("Myeloma", 3, line=1, cex=1)
 mtext("Months from Day +20 Landmark", 1, line=3.5, cex=1)
 mtext("Overall Survival", 2, line=3, cex=1) 
 
 nlgd = length(lgd)
 mtext("No. at risk",side=1,line=5,at=legendxpos,font=2,col=1,cex=1)
 for (i in 1:nlgd){
    mtext(lgd[i],side=1,line=(5+i),at=legendxpos,
          font=4,col=(i+1),cex=0.8)
 } 
 
 
 sout=s2
 maxtlist=xmax/stepatrisk   
 slist=rep(NA,maxtlist)
 
 for (k in 1:nlgd) {
    mtext(sout[k]$n.risk[1],side=1,line=(5+k),at=0,cex=0.8)
    for (i in 1:maxtlist) {
       j=i*stepatrisk;
       slist[i] = sout[k]$n.risk[which(sout[k]$time>=j)[1]]
       slist[i] = ifelse(is.na(slist[i]),0,slist[i])
       mtext(slist[i],side=1,line=(5+k),at=j,cex=0.8)
    }
 }
 
 
 
 
 s2 = survfit(Surv(OSlandmark1 ,OSEventlandmark1 )~SimpsonMED, data=combineddata, subset=(MyelomaAmyloid ==0 )) 
 p1 <-  summary(coxph(Surv(OSlandmark1,OSEventlandmark1 )~SimpsonMED +strata(institution), data=combineddata, subset=(MyelomaAmyloid ==0 )))$coefficients[,5]
 
 lgd=c("Below Median","Above Median")
 
 par(mar=c(8,12,4.1,2.1))
 plot(s2,col=c(2:6), conf.int=F, xlim=c(0,xmax),ylim=c(0,1.05),yaxt="n", xaxs="i",yaxs="i",   bty="n",
      xaxt="n",frame.plot=F,xlab="",ylab="",cex.lab=1.2, mark.time=TRUE,
      main="",cex=.6 , lwd=2)
 axis(side=1,seq(0,xmax,by=2))
 axis(side=2,seq(0,1,by=0.2))
 text(xmax*0.75, 0.05, paste("P-value = ", round(p1, 3), sep=""))
 
 
 mtext("Lymphoma", 3, line=1, cex=1)
 mtext("Months from Day +20 Landmark", 1, line=3.5, cex=1)
 mtext("Overall Survival", 2, line=3, cex=1) 
 
 nlgd = length(lgd)
 mtext("No. at risk",side=1,line=5,at=legendxpos,font=2,col=1,cex=1)
 for (i in 1:nlgd){
    mtext(lgd[i],side=1,line=(5+i),at=legendxpos,
          font=4,col=(i+1),cex=0.8)
 } 
 
 
 sout=s2
 maxtlist=xmax/stepatrisk   
 slist=rep(NA,maxtlist)
 
 for (k in 1:nlgd) {
    mtext(sout[k]$n.risk[1],side=1,line=(5+k),at=0,cex=0.8)
    for (i in 1:maxtlist) {
       j=i*stepatrisk;
       slist[i] = sout[k]$n.risk[which(sout[k]$time>=j)[1]]
       slist[i] = ifelse(is.na(slist[i]),0,slist[i])
       mtext(slist[i],side=1,line=(5+k),at=j,cex=0.8)
    }
 }
 

 
 dev.off()
 
 
 
 #############################################################  new GEE model 
 #############################################################   
 #############################################################
 #############################################################
 #############################################################
 
 
 #sampledata <- read.csv("msk duke sdtable 09 05 19.csv")
# sampledata <- read.csv("msk duke sdtable 01 18 20.csv")
 sampledata <- read.csv("../msk duke sdtable 03 02 20.csv")
 
 dim(sampledata)
 dim(ptdata)
 
 sampledata <- sampledata[sampledata$day.post.txp >= -10 & sampledata$day.post.txp <= 30,]
 dim(sampledata)
 
 #sampledata$Log_simpson_reciprocal <- log( sampledata$simpson_reciprocal)
 
library(splines)

library(geepack)


mod1 <- geese(simpson_reciprocal~  day.post.txp*institution,id=mrn ,data=sampledata, family=gaussian(),cor.link = "identity")
summary(mod1) 

paste(round(mod1$beta,3), " (", round( mod1$beta -1.96*sqrt(diag(mod1$vbeta)),3) , " to ", round( mod1$beta+ 1.96*sqrt(diag(mod1$vbeta)),3) , ")", sep="")


sampledata$institutionB <- 1*(sampledata$institution=="Duke Auto")
mod1 <- geese(simpson_reciprocal~  day.post.txp*institutionB,id=mrn ,data=sampledata, family=gaussian(),cor.link = "identity")
summary(mod1) 

paste(round(mod1$beta,3), " (", round( mod1$beta -1.96*sqrt(diag(mod1$vbeta)),3) , " to ", round( mod1$beta+ 1.96*sqrt(diag(mod1$vbeta)),3) , ")", sep="")



mod1 <- geese(simpson_reciprocal~  day.post.txp,id=mrn ,data=sampledata, family=gaussian(),cor.link = "identity", 
              subset=( institution =="Duke Auto"))
summary(mod1) 

mod1 <- geese(simpson_reciprocal~  day.post.txp,id=mrn ,data=sampledata, family=gaussian(),cor.link = "identity", 
              subset=( institution =="MSK Auto"))
summary(mod1)


table(sampledata$institution)



 #############################################################  consistency by simpson_reciprocal
 #############################################################   
 #############################################################
 #############################################################
 #############################################################
 
  
#sampledata <- read.csv("msk duke sdtable 01 18 20.csv")
sampledata <- read.csv("../msk duke sdtable 03 02 20.csv")


sampledata <- sampledata[sampledata$day.post.txp >= 5 & sampledata$day.post.txp <= 15, ]
 dim(sampledata)
 dim(ptdata)
 sampledata$consistency <- as.character(sampledata$consistency)
 sampledata <- sampledata[!is.na(sampledata$consistency),]
 sampledata <- sampledata[ sampledata$consistency != "",]
 
 table(sampledata$consistency, exclude=NULL)
 
 
 
 ptdataREDUCED <- ptdata[, c("mrn", "Melphalan", "op.vs.ip" , "MyelomaAmyloid")]
 dim(sampledata)
 dim(ptdataREDUCED)
 
 consistencyanalysis <- merge( sampledata,ptdataREDUCED )
 dim(consistencyanalysis)
 
 ORDoutcome <- rep(NA, length(consistencyanalysis$consistency ))
 ORDoutcome[consistencyanalysis$consistency %in% c("Formed stool", "formed stool")] <- 1
 ORDoutcome[consistencyanalysis$consistency == "semi-formed"] <- 2
 ORDoutcome[consistencyanalysis$consistency == "liquid"] <- 3
 
 table(consistencyanalysis$consistency,ORDoutcome, exclude=NULL)
 
 consistencyanalysis$ORDoutcome <-  ORDoutcome 
 consistencyanalysis$Log_simpson_reciprocal <- log( consistencyanalysis$simpson_reciprocal)
 
 install.packages("geepack")
library(geepack)
mod1 <-summary(ordgee(ordered(ORDoutcome) ~ Log_simpson_reciprocal+day.post.txp ,id=mrn ,  data=consistencyanalysis,corstr = "independence"))  
mod2 <- summary(ordgee(ordered(ORDoutcome) ~ Log_simpson_reciprocal +MyelomaAmyloid+day.post.txp,id=mrn ,  data=consistencyanalysis,corstr = "independence"))

consistencyanalysis$ORDoutcomeCAT <- ordered(consistencyanalysis$ORDoutcome,
                                 levels = c(1,2, 3),
                                 labels = c("Formed, 96 Samples", "Semi-formed, 169 Samples", "Liquid, 180 Samples"))




pdf(paste("StoolConsistency.",Sys.Date(),".pdf",sep=""), height=7, width=7)

boxplot(Log_simpson_reciprocal~ ORDoutcomeCAT,  data=consistencyanalysis,  main="Stool Consistency \n Days +5 to +15", ylim=c(0,4),
        xlab="", ylab="Simpson's reciprocal (log-transformed)",frame=FALSE, xaxt="n")

axis(side=1, line=0, c(1,2,3), c("Formed \n  96 Samples", "Semi-formed \n  169 Samples", "Liquid  \n 180 Samples"), tick=FALSE)

posH=4
text(0.5,4, "GEE for Ordinal Consistency Outcome", pos=posH, cex=0.8)
text(0.5,3.8, "O.R. per log diversity: 0.57 (95% CI: 0.44-0.73); P-value: <0.001", pos=posH, cex=0.8)
text(0.7,3.7, "Adjusted for HCT day", pos=posH, cex=0.8)
text(0.5,3.5, "O.R. per log diversity: 0.73 (95% CI: 0.55-0.97); P-value: 0.028", pos=posH, cex=0.8)
text(0.7,3.4, "Adjusted for HCT day and disease", pos=posH, cex=0.8)

dev.off()
 
 

#############################################################   below is not updated.    
#############################################################   below is not updated.        
#############################################################   below is not updated.  
#############################################################   below is not updated.  
#############################################################   below is not updated.  

#############################################################   below is not updated.    
#############################################################   below is not updated.        
#############################################################   below is not updated.  
#############################################################   below is not updated.  
#############################################################   below is not updated.  

#############################################################   below is not updated.    
#############################################################   below is not updated.        
#############################################################   below is not updated.  
#############################################################   below is not updated.  
#############################################################   below is not updated.  

#############################################################   below is not updated.    
#############################################################   below is not updated.        
#############################################################   below is not updated.  
#############################################################   below is not updated.  
#############################################################   below is not updated.  

#############################################################  
#############################################################      
#############################################################   same analysis, but including 
#############################################################   non-dominate
#############################################################

monodomination <- read.csv("monodomination_by_consistency 07 01 19.csv")
table(monodomination$butyrate)

monodomination <- read.csv("monodomination by consistency 09 08 19.csv")
table(monodomination$butyrate[monodomination$dominant != "not-dominant"])


monodomination <- read.csv("monodomination by consistency 09 08 19 for SD.csv")



monodomination <- monodomination[!is.na(monodomination$butyrate),]
monodomination$butyrate <- as.character(monodomination$butyrate)
#monodomination <-  monodomination[monodomination$butyrate != "non-dominant",]




ORDoutcome <- rep(NA, length(monodomination$consistency ))
ORDoutcome[monodomination$consistency %in% c("Formed stool", "formed stool")] <- 1
ORDoutcome[monodomination$consistency == "semi-formed"] <- 2
ORDoutcome[monodomination$consistency == "liquid"] <- 3
monodomination$ORDoutcome <- ORDoutcome


monodomination$NonDom <- 1*(monodomination$butyrate == "not-dominant")
monodomination$butyrateDom <- 1*(monodomination$butyrate == "butyrate")

summary(ordgee(ordered(ORDoutcome) ~ NonDom +butyrateDom+day_relative_to_hct_bin,id=patient_id ,  data=monodomination,corstr = "independence")) 

 
prop.table(table(monodomination$butyrate, monodomination$consistency), 2)

r1B <- 100*mean(monodomination$butyrate[monodomination$ORDoutcome==1] == "butyrate", na.rm = TRUE)
r2B <- 100*mean(monodomination$butyrate[monodomination$ORDoutcome==2] == "butyrate", na.rm = TRUE)
r3B <- 100*mean(monodomination$butyrate[monodomination$ORDoutcome==3] == "butyrate", na.rm = TRUE)

r1NB <- 100*mean(monodomination$butyrate[monodomination$ORDoutcome==1] == "non butyrate", na.rm = TRUE)
r2NB <- 100*mean(monodomination$butyrate[monodomination$ORDoutcome==2] == "non butyrate", na.rm = TRUE)
r3NB <- 100*mean(monodomination$butyrate[monodomination$ORDoutcome==3] == "non butyrate", na.rm = TRUE)

pdf(paste("ConsistencyByButyrateDominationALLSamples.",Sys.Date(),".pdf",sep=""), height=8, width=8)
factH <- 0.4
plot(c(0,0),c(0,0), typ="n", ylab="Percentage of Samples (%)", xlab="", ylim=c(0,110),xlim=c(factH,5.25), yaxt="n",xaxt="n",bty="n")

axis(side=1, pos=-1, c(1,2,3), c("Formed \n  240 Samples", "Semi-formed \n  186 Samples", "Liquid  \n 137 Samples"), tick=FALSE)
axis(side=2, seq(0,100, by=20))

rect(1-factH,0,1+factH,100, col="gainsboro")
rect(1-factH,0,1+factH,r1B+r1NB, col="darkorange")
rect(1-factH,0,1+factH,r1B, col="blue")
text(1, r1B+r1NB+2, paste("+", round(r1NB),"%" ,sep="") ,col="darkorange")
text(1, r1B+2, paste("+",round(r1B), "%",sep=""), col="blue")

rect(2-factH,0,2+factH,100, col="gainsboro")
rect(2-factH,0,2+factH,r2B+r2NB, col="darkorange")
rect(2-factH,0,2+factH,r2B, col="blue")
text(2, r2B+r2NB+2, paste("+",round(r2NB), "%",sep=""), col="darkorange")
text(2, r2B+2, paste("+",round(r2B), "%",sep=""), col="blue")

 
rect(3-factH,0,3+factH,100, col="gainsboro")
rect(3-factH,0,3+factH,r3B+r3NB, col="darkorange")
rect(3-factH,0,3+factH,r3B, col="blue")
text(3, r3B+r3NB+2, paste("+",round(r3NB), "%",sep=""), col="darkorange")
text(3, r3B+2, paste("+",round(r3B), "%",sep=""), col="blue")

legend(3.4, 20, rev(c("Butyrate Dominated","Non-butyrate Dominated", "Non-dominated")), col=rev(c("blue", "darkorange","grey")),
       fill=rev(c("blue", "darkorange","grey")), bty="n", cex=1.15)
 
 
dev.off()
 






