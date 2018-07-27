#####################################################
# Comparison between RNASeq and Arrays for AS detection
#####################################################
library(EventPointer)
library(aroma.affymetrix)
library(limma)



split3 <- function(gr) {
  
  exons <- (gr$type == "exon")
  
  # Refs
  Refs <- grepl("^Ref", gr$transcript_id)
  onlyexons <-  exons & Refs
  Refgrl <-  split(gr[onlyexons,], gr[onlyexons,]$transcript_id)
  
  # Path A
  Path_A <- grepl("^A", gr$transcript_id)
  onlyexons <-  exons & Path_A
  PathAgrl <-  split(gr[onlyexons,], gr[onlyexons,]$transcript_id)
  
  # Path B
  Path_B <- grepl("^B", gr$transcript_id)
  onlyexons <-  exons & Path_B
  PathBgrl <-  split(gr[onlyexons,], gr[onlyexons,]$transcript_id)
  
  # Return a list with the three of them
  return(list(Refgrl = Refgrl, PathAgrl = PathAgrl, PathBgrl=PathBgrl))
}


ProspliceMatching<-function(HTAfile,RNAfile,Class,Events_HTA_Prosplice,Events_RNA_Prosplice)
{
  
  RNAgr <- import.gff(RNAfile) 
  HTAgr <- import.gff(HTAfile) 
  
  RNAgrl3 <- split3(RNAgr)
  HTAgrl3 <- split3(HTAgr)
  
  ## A-A match
  AA_RH <- findSpliceOverlaps(RNAgrl3$PathAgrl, HTAgrl3$PathAgrl, ignore.strand = T)
  
  
  hits <- as.data.frame(AA_RH)
  
  M_AA_RH <- sparseMatrix(hits$queryHits[hits$compatible], hits$subjectHits[hits$compatible], x = 1, 
                          dims = c(length(RNAgrl3$PathAgrl), length(HTAgrl3$PathAgrl)))
  
  M_AA_Light <- sparseMatrix(hits$queryHits, hits$subjectHits, x = 1, 
                             dims = c(length(RNAgrl3$PathAgrl), length(HTAgrl3$PathAgrl)))
  
  AA_HR <- findSpliceOverlaps(HTAgrl3$PathAgrl, RNAgrl3$PathAgrl, ignore.strand = T)
  
  
  hits <- as.data.frame(AA_HR)
  M_AA_HR <- t(sparseMatrix(hits$queryHits[hits$compatible], hits$subjectHits[hits$compatible], x = 1, 
                            dims = c(length(HTAgrl3$PathAgrl),length(RNAgrl3$PathAgrl))))
  
  M_AA <- M_AA_HR + M_AA_RH
  M_AA@x[M_AA@x>1] <- 1
  
  rownames(M_AA)<-names(RNAgrl3$PathAgrl)
  colnames(M_AA)<-names(HTAgrl3$PathAgrl)
  
  rownames(M_AA_Light)<-names(RNAgrl3$PathAgrl)
  colnames(M_AA_Light)<-names(HTAgrl3$PathAgrl)
  
  ## B-B match
  BB_RH <- findSpliceOverlaps(RNAgrl3$PathBgrl, HTAgrl3$PathBgrl, ignore.strand = T)
  hits <- as.data.frame(BB_RH)
  M_BB_RH <- sparseMatrix(hits$queryHits[hits$compatible], hits$subjectHits[hits$compatible], x = 1, 
                          dims = c(length(RNAgrl3$PathBgrl), length(HTAgrl3$PathBgrl)))
  
  M_BB_Light <- sparseMatrix(hits$queryHits, hits$subjectHits, x = 1, 
                             dims = c(length(RNAgrl3$PathBgrl), length(HTAgrl3$PathBgrl)))
  
  BB_HR <- findSpliceOverlaps(HTAgrl3$PathBgrl, RNAgrl3$PathBgrl, ignore.strand = T)
  hits <- as.data.frame(BB_HR)
  M_BB_HR <- t(sparseMatrix(hits$queryHits[hits$compatible], hits$subjectHits[hits$compatible], x = 1, 
                            dims = c(length(HTAgrl3$PathBgrl),length(RNAgrl3$PathBgrl))))
  
  M_BB <- M_BB_HR + M_BB_RH
  M_BB@x[M_BB@x>1] <- 1
  
  rownames(M_BB)<-names(RNAgrl3$PathAgrl)
  colnames(M_BB)<-names(HTAgrl3$PathAgrl)
  
  rownames(M_BB_Light)<-names(RNAgrl3$PathAgrl)
  colnames(M_BB_Light)<-names(HTAgrl3$PathAgrl)
  
  ## R-R match
  
  RR_RH <- findOverlaps(RNAgrl3$Refgrl, HTAgrl3$Refgrl, ignore.strand = T, type = "any")
  hits <- as.data.frame(RR_RH)
  M_RR <- sparseMatrix(hits[,1], hits[,2], x = 1,
                       dims = c(length(RNAgrl3$Refgrl), length(HTAgrl3$Refgrl)))
  
  M_RR@x[M_RR@x>1] <- 1
  
  
  EventosRNA <- unlist(lapply(strsplit(names(RNAgrl3$Refgrl),"_"),function(x) {paste(tail(x,2),collapse="_")}))
  EventosRNAunique <- unique(EventosRNA)
  ix  <- match(EventosRNA,EventosRNAunique)
  iy <- 1:length(EventosRNA)
  EventsRNAXRef <- sparseMatrix(ix, iy, x = 1,
                                dims = c(length(EventosRNAunique), length(EventosRNA)))
  colnames(EventsRNAXRef) <- names(RNAgrl3$Refgrl)
  rownames(EventsRNAXRef) <- EventosRNAunique
  
  
  
  EventosHTA <- unlist(lapply(strsplit(names(HTAgrl3$Refgrl),"_"),function(x) {paste(tail(x,2),collapse="_")}))
  EventosHTAunique <- unique(EventosHTA)
  ix  <- match(EventosHTA,EventosHTAunique)
  iy <- 1:length(EventosHTA)
  EventsHTAXRef <- sparseMatrix(ix, iy, x = 1,
                                dims = c(length(EventosHTAunique), length(EventosHTA)))
  colnames(EventsHTAXRef) <- names(HTAgrl3$Refgrl)
  rownames(EventsHTAXRef) <- EventosHTAunique
  
  
  M_RR<-(EventsRNAXRef%*%M_RR)%*%t(EventsHTAXRef)
  M_RR@x[M_RR@x>1] <- 1
  
  NamesRefRNA<-unlist(lapply(strsplit(names(RNAgrl3$Refgrl),"_"),function(x) {paste(tail(x,2),collapse="_")}))
  NamesRNA<-unlist(lapply(strsplit(names(RNAgrl3$PathAgrl),"_"),function(x) {paste(tail(x,2),collapse="_")}))
  
  NamesRefHTA<-unlist(lapply(strsplit(names(HTAgrl3$Refgrl),"_"),function(x) {paste(tail(x,2),collapse="_")}))
  NamesHTA<-unlist(lapply(strsplit(names(HTAgrl3$PathAgrl),"_"),function(x) {paste(tail(x,2),collapse="_")}))
  
  
  RA_RH <- findOverlaps(RNAgrl3$Refgrl, HTAgrl3$PathAgrl, ignore.strand = T, type = "any",minoverlap = 3)
  hits <- as.data.frame(RA_RH)
  M_RA_RH <- sparseMatrix(hits[,1], hits[,2], x = 1,
                          dims = c(length(RNAgrl3$Refgrl), length(HTAgrl3$PathAgrl)))
  rownames(M_RA_RH)<-NamesRefRNA
  colnames(M_RA_RH)<-NamesHTA
  
  RB_RH <- findOverlaps(RNAgrl3$Refgrl, HTAgrl3$PathBgrl, ignore.strand = T, type = "any",minoverlap = 3)
  hits <- as.data.frame(RB_RH)
  M_RB_RH <- sparseMatrix(hits[,1], hits[,2], x = 1,
                          dims = c(length(RNAgrl3$Refgrl), length(HTAgrl3$PathAgrl)))
  rownames(M_RB_RH)<-NamesRefRNA
  colnames(M_RB_RH)<-NamesHTA
  
  RA_HR<- findOverlaps(HTAgrl3$Refgrl, RNAgrl3$PathAgrl, ignore.strand = T, type = "any",minoverlap = 3)
  hits <- as.data.frame(RA_HR)
  M_RA_HR <- sparseMatrix(hits[,1], hits[,2], x = 1,
                          dims = c(length(HTAgrl3$Refgrl), length(RNAgrl3$PathAgrl)))
  rownames(M_RA_HR)<-NamesRefHTA
  colnames(M_RA_HR)<-NamesRNA
  
  RB_HR <- findOverlaps(HTAgrl3$Refgrl, RNAgrl3$PathBgrl, ignore.strand = T, type = "any",minoverlap = 3)
  hits <- as.data.frame(RB_HR)
  M_RB_HR <- sparseMatrix(hits[,1], hits[,2], x = 1,
                          dims = c(length(HTAgrl3$Refgrl), length(RNAgrl3$PathAgrl)))
  rownames(M_RB_HR)<-NamesRefHTA
  colnames(M_RB_HR)<-NamesRNA
  
  M_RA_RH<-EventsRNAXRef%*%M_RA_RH
  M_RB_RH<-EventsRNAXRef%*%M_RB_RH
  
  M_RAB_RH<-M_RA_RH+M_RB_RH
  
  M_RA_HR<-EventsHTAXRef%*%M_RA_HR
  M_RB_HR<-EventsHTAXRef%*%M_RB_HR
  
  M_RAB_HR<-t(M_RA_HR+M_RB_HR)
  
  M_RAB<-M_RAB_RH+M_RAB_HR
  
  rownames(M_RAB)<-names(RNAgrl3$PathAgrl)
  colnames(M_RAB)<-names(HTAgrl3$PathAgrl)
  
  M_RAB@x[M_RAB@x>1] <- 1
  ##############################################################
  
  
  ## A-B match
  AB_RH <- findSpliceOverlaps(RNAgrl3$PathAgrl, HTAgrl3$PathBgrl, ignore.strand = T)
  hits <- as.data.frame(AB_RH)
  M_AB_RH <- sparseMatrix(hits$queryHits[hits$compatible], hits$subjectHits[hits$compatible], x = 1, 
                          dims = c(length(RNAgrl3$PathAgrl), length(HTAgrl3$PathBgrl)))
  
  M_AB_Light <- sparseMatrix(hits$queryHits, hits$subjectHits, x = 1, 
                             dims = c(length(RNAgrl3$PathAgrl), length(HTAgrl3$PathBgrl)))
  
  AB_HR <- findSpliceOverlaps(HTAgrl3$PathBgrl, RNAgrl3$PathAgrl, ignore.strand = T)
  hits <- as.data.frame(AB_HR)
  M_AB_HR <- t(sparseMatrix(hits$queryHits[hits$compatible], hits$subjectHits[hits$compatible], x = 1, 
                            dims = c(length(HTAgrl3$PathAgrl),length(RNAgrl3$PathBgrl))))
  
  M_AB <- M_AB_HR + M_AB_RH
  M_AB@x[M_AB@x>1] <- 1
  
  rownames(M_AB)<-names(RNAgrl3$PathAgrl)
  colnames(M_AB)<-names(HTAgrl3$PathBgrl)
  
  rownames(M_AB_Light)<-names(RNAgrl3$PathAgrl)
  colnames(M_AB_Light)<-names(HTAgrl3$PathBgrl)
  
  ## B-A match
  BA_RH <- findSpliceOverlaps(RNAgrl3$PathBgrl, HTAgrl3$PathAgrl, ignore.strand = T)
  hits <- as.data.frame(BA_RH)
  M_BA_RH <- sparseMatrix(hits$queryHits[hits$compatible], hits$subjectHits[hits$compatible], x = 1, 
                          dims = c(length(RNAgrl3$PathBgrl), length(HTAgrl3$PathAgrl)))
  
  M_BA_Light <- sparseMatrix(hits$queryHits, hits$subjectHits, x = 1, 
                             dims = c(length(RNAgrl3$PathBgrl), length(HTAgrl3$PathAgrl)))
  
  BA_HR <- findSpliceOverlaps(HTAgrl3$PathAgrl, RNAgrl3$PathBgrl, ignore.strand = T)
  hits <- as.data.frame(BA_HR)
  M_BA_HR <- t(sparseMatrix(hits$queryHits[hits$compatible], hits$subjectHits[hits$compatible], x = 1, 
                            dims = c(length(HTAgrl3$PathBgrl),length(RNAgrl3$PathAgrl))))
  
  M_BA <- M_BA_HR + M_BA_RH
  M_BA@x[M_BA@x>1] <- 1
  
  rownames(M_BA)<-names(RNAgrl3$PathBgrl)
  colnames(M_BA)<-names(HTAgrl3$PathAgrl)
  
  rownames(M_BA_Light)<-names(RNAgrl3$PathBgrl)
  colnames(M_BA_Light)<-names(HTAgrl3$PathAgrl)
  
  
  ### Choose Which HT to work with
  
  if(Class=="Normal")
  {
    Good <- (M_AA * M_BB * M_RR)-M_RAB
    Good2 <- (M_AB * M_BA * M_RR)-M_RAB
    
  }else if(Class=="LooseLight")
  {
    Good <- (M_AA_Light * M_BB_Light)-M_RAB
    Good2 <- (M_AB_Light * M_BA_Light)-M_RAB
    
  }else if(Class=="Loose")
  {
    Good <- (M_AA_Light * M_BB_Light* M_RR)-M_RAB
    Good2 <- (M_AB_Light * M_BA_Light* M_RR)-M_RAB
    
  }
  
  
  HT1  <- which(as.matrix(Good > 0), arr.ind = T)
  HT2<- which(as.matrix(Good2>0), arr.ind = T)
  HT<-unique(rbind(HT1,HT2))
  
  
  NamesRNA<-names(RNAgrl3$Refgrl)
  NamesRNA<-sapply(NamesRNA,function(X){A<-tail(unlist(strsplit(X,"_")),2);A<-paste(A[1],A[2],sep="_");return(A)})
  names(NamesRNA)<-NULL
  NamesHTA<-names(HTAgrl3$Refgrl)
  NamesHTA<-sapply(NamesHTA,function(X){A<-tail(unlist(strsplit(X,"_")),2);A<-paste(A[1],A[2],sep="_");return(A)})
  names(NamesHTA)<-NULL
  RNA<-NamesRNA[HT[,1]]
  HTA<-NamesHTA[HT[,2]]
  
  
  if(Class=="Normal")
  {
    XX<-c(rep("AA_BB",nrow(HT1)),rep("AB_BA",nrow(HT2)))
    
  }else if((Class=="Loose") | (Class=="LooseLight"))
  {
    XX<-rep("AA_BB",nrow(HT))
    
  }
  
  Prosplice_Matching<-cbind(HTA,RNA)
  HTA_Match<-match(Prosplice_Matching[,1],rownames(Events_HTA_Prosplice))
  RNA_Match<-match(Prosplice_Matching[,2],rownames(Events_RNASeq_Prosplice))
  
  
  RNASeqID<-rownames(Events_RNASeq_Prosplice)[RNA_Match]
  HTAID<-rownames(Events_HTA_Prosplice)[HTA_Match]
  
  Prosplice_Matching<-cbind(HTAID,Events_HTA_Prosplice[HTA_Match,],HTA_Match,HTA_Match/nrow(Events_HTA_Prosplice),RNASeqID,Events_RNASeq_Prosplice[RNA_Match,],RNA_Match,RNA_Match/nrow(Events_RNASeq_Prosplice),XX)
  
  colnames(Prosplice_Matching)<-c("ID_HTA","Gene_HTA","EventType_HTA","Region_HTA","Zvalue_HTA","Pvalue_HTA","PSI_HTA","Rank_HTA","RankPercent_HTA",
                                  "ID_RNA","Gene_RNA","EventType_RNA","Region_RNA","Pvalue_RNA","Zvalue_RNA","PSI_RNA","Rank_RNA","RankPercent_RNA","Match_Class")
  
  Prosplice_Matching[which(Prosplice_Matching[,"Match_Class"]=="AB_BA"),"Zvalue_RNA"]<-Prosplice_Matching[which(Prosplice_Matching[,"Match_Class"]=="AB_BA"),"Zvalue_RNA"]*-1
  Prosplice_Matching[which(Prosplice_Matching[,"Match_Class"]=="AB_BA"),"PSI_RNA"]<-Prosplice_Matching[which(Prosplice_Matching[,"Match_Class"]=="AB_BA"),"PSI_RNA"]*-1
  
  return(Prosplice_Matching)
  
  
  
}

GetFDR<-function(Pvalues,alpha=0.001,span=c(0.3,0.5))
{
  malos<-sum((Pvalues>span[1]) & (Pvalues<span[2]))
  buenos<-sum((Pvalues<alpha))
  FDR<-malos/(span[2]-span[1])*(alpha /buenos)
  Pi_zero<-(malos/(span[2]-span[1]))/length(Pvalues)
  return(c(FDR,Pi_zero))
}


######################
# Download Data
######################
## In this section, we specfiy the accession codes to download the data for both
## arrays and RNASeq. Also, the corresponding dropbox link to download specific files
## required to run EventPointer with microarrays.

# All data is available in a SuperSeries at GEO with accession code GSE104974

# For microarrays, this is the corresponding dropbox link:
# https://www.dropbox.com/sh/wpwz1jx0l112icw/AAD4yrEY4HG1fExUmtoBmrOWa/HTA%202.0?dl=0

######################
# Directory Structure
######################

## We assume that the files are organized in the following way:

# /SplicingComparison/
# -------------------/RNASeq
# -------------------/------/BAM/  All .BAM files of the experiment
# -------------------/Microarrays
# -------------------/-----------/EP_Files/   Files downloaded from dropbox
# -------------------/-----------/aroma/  Required directory structure for aroma affymetrix (http://www.aroma-project.org/docs/HowDataFilesAndDataSetsAreLocated/)

######################
# Event Pointer RNASeq
######################


  
Samples<-dir("/SplicingComparison/RNASeq/BAM/",pattern=".bam$")
  
SamplePath<-"/SplicingComparison/RNASeq/BAM/"
fileTransc<-"/SplicingComparison/Microarrays/EP_Files/HTA-2_0_GTF_Transcripts.gtf"
TxtPath<-"/SplicingComparison/RNASeq/"


### This matrices are for all cell lines together
SS<-matrix(unlist(strsplit(Samples,"_")),ncol=5,byrow=T)[,1:2]
colnames(SS)<-c("CellLine","Treatment")
SS<-as.data.frame(SS)
Design<-model.matrix(~.,data=SS)
Design[,4]<-(Design[,4]+1)%%2
Design[,3]<-(Design[,1]+Design[,2]+Design[,3])%%2
Design<-Design[,c(1,3,2,4)]
colnames(Design)<-c("(Intercept)","CellLineMDA231","CellLineMDA468","TreatmentCX4945")
Contrast<-t(t(c(0,0,0,1)))

SG_RNASeq<-PrepareBam_EP(Samples,SamplePath,Ref_Transc="GTF",fileTransc=fileTransc,cores=5,Alpha=2)
AllEvents_RNASeq_EP<-EventDetection(SG_RNASeq,15,TxtPath)
Events_RNASeq_Prosplice<-EventPointer_RNASeq(AllEvents_RNASeq_EP,Design,Contrast,Statistic="LogFC",PSI=TRUE)

## Generate GTF files

EventPointer_RNASeq_IGV(Events_RNASeq_Prosplice, SG_RNASeq, "/SplicingComparison/RNASeq/EventsFound_RNASeq.txt", "/SplicingComparison/RNASeq/")

######################
# Event Pointer Arrays
######################

setwd("/SplicingComparison/Microarrays/aroma/")

verbose <- Arguments$getVerbose(-8);
timestampOn(verbose);
projectName <- "Prosplice"
chipType <- "EP_HTA-2_0"
cdfGFile <- "EP_HTA-2_0,r"
cdfG <- AffymetrixCdfFile$byChipType(cdfGFile)
cs <- AffymetrixCelSet$byName(projectName, cdf=cdfG)
bc <- NormExpBackgroundCorrection(cs, method="mle", tag=c("*","r11"));
csBC <- process(bc,verbose=verbose,ram=20);
qn <- QuantileNormalization(csBC, typesToUpdate="pm");
csN <- process(qn,verbose=verbose,ram=20);
plmEx <- ExonRmaPlm(csN, mergeGroups=FALSE)
fit(plmEx, verbose=verbose)
cesEx <- getChipEffectSet(plmEx)

### This matrices are for all cell lines together

# Design Matrix 
ExFit <- extractDataFrame(cesEx, addNames = TRUE)
nombres  <- colnames(ExFit)[6:ncol(ExFit)] 
Dmatrix <- matrix(1,nrow = length(nombres), ncol=4) 
Dmatrix[,2] <- grepl("MDA231",nombres, fixed=TRUE) 
Dmatrix[,3] <- grepl("MDA468",nombres, fixed=TRUE) 
Dmatrix[,4] <- grepl("CX4945",nombres, fixed=TRUE) 

# Contrast Matrix
Cmatrix <- t(t(c(0,0,0,1)))


# Main EventPointer Function

Events_HTA_Prosplice<-EventPointer(Design=Dmatrix,Contrast=Cmatrix,ExFit=ExFit,Eventstxt=EventstxtFile,Filter=T,Qn=0.25,Statistic="LogFC",PSI=T)

# Generate GTF files required to perform the matching between events 
EventPointer_IGV(Events_HTA_Prosplice, "CustomGTF", inputFile = "/SplicingComparison/Microarrays/EP_Files/HTA-2_0_GTF_Transcripts.gtf", 
                 PSR="/SplicingComparison/Microarrays/EP_Files/HTA-2_0_PSR_Probes.txt", 
                 Junc="/SplicingComparison/Microarrays/EP_Files/HTA-2_0_Junc_Probes.txt", 
                 PathGTF="/SplicingComparison/Microarrays/", 
                 EventsFile="/SplicingComparison/Microarrays/EP_Files/EventsFound.txt", 
                 microarray = "HTA-2_0")

#############################
# Matching between events
#############################

library(rtracklayer)
library(GenomicAlignments)
library(SRAdb)
library(Matrix)


## Arrays > 25 & RNASeq alpha = 2 (Normal Match Criterion)

## The gtf files generated by EventPointer arent sorted based on genomic coordinates,
## we sorted them using samtools sort command

HTAfile<-"/SplicingComparison/Microarrays/paths_Arrays.sorted.gtf"
RNAfile<-"/SplicingComparison/RNASeq/paths_RNASeq.sorted.gtf"

Prosplice_Matching<-ProspliceMatching(HTAfile,RNAfile,"Normal",Events_HTA_Prosplice,Events_RNA_Prosplice)

###########################
# Separte into groups
############################


G1H<-rownames(Events_HTA_Prosplice)[which(Events_HTA_Prosplice[,"Splicing Pvalue"]<1e-3)]
G1H<-Events_HTA_Prosplice[G1H[which(is.na(match(G1H,Prosplice_Matching[,"ID_HTA"])))],]

G1R<-rownames(Events_RNASeq_Prosplice)[which(Events_RNASeq_Prosplice[,"Pvalue"]<1e-3)]
G1R<-Events_RNASeq_Prosplice[which(is.na(match(G1R,Prosplice_Matching[,"ID_RNA"]))),]

G4<-Prosplice_Matching[which(Prosplice_Matching[,"Pvalue_HTA"]<1e-3 & 
                               Prosplice_Matching[,"Pvalue_RNA"]<1e-3 &
                               sign(Prosplice_Matching[,"Zvalue_HTA"])==sign(Prosplice_Matching[,"Zvalue_RNA"]) &
                               abs(Prosplice_Matching[,"PSI_HTA"]) > 0.1 &
                               abs(Prosplice_Matching[,"PSI_RNA"]) > 0.1),]

G5<-Prosplice_Matching[which(Prosplice_Matching[,"Pvalue_HTA"]<1e-3 & 
                               Prosplice_Matching[,"Pvalue_RNA"]<1e-3 &
                               sign(Prosplice_Matching[,"Zvalue_HTA"])!=sign(Prosplice_Matching[,"Zvalue_RNA"])),]

G2H<-Prosplice_Matching[which(Prosplice_Matching[,"Pvalue_HTA"]<1e-3 & 
                                Prosplice_Matching[,"Pvalue_RNA"]>0.2 &
                                abs(Prosplice_Matching[,"PSI_HTA"]) > 0.1 &
                                sign(Prosplice_Matching[,"Zvalue_HTA"])==sign(Prosplice_Matching[,"Zvalue_RNA"])),]

G3H<-Prosplice_Matching[which(Prosplice_Matching[,"Pvalue_HTA"]<1e-3 & 
                                Prosplice_Matching[,"Pvalue_RNA"]>0.2 &
                                abs(Prosplice_Matching[,"PSI_HTA"]) > 0.1 &
                                sign(Prosplice_Matching[,"Zvalue_HTA"])!=sign(Prosplice_Matching[,"Zvalue_RNA"])),]

G2R<-Prosplice_Matching[which(Prosplice_Matching[,"Pvalue_RNA"]<1e-3 & 
                                Prosplice_Matching[,"Pvalue_HTA"]>0.2 &
                                abs(Prosplice_Matching[,"PSI_RNA"]) > 0.1 &
                                sign(Prosplice_Matching[,"Zvalue_RNA"])==sign(Prosplice_Matching[,"Zvalue_HTA"])),]


G3R<-Prosplice_Matching[which(Prosplice_Matching[,"Pvalue_RNA"]<1e-3 & 
                                Prosplice_Matching[,"Pvalue_HTA"]>0.2 &
                                abs(Prosplice_Matching[,"PSI_RNA"]) > 0.1 &
                                sign(Prosplice_Matching[,"Zvalue_RNA"])!=sign(Prosplice_Matching[,"Zvalue_HTA"])),]

G4<-G4[order(G4[,"RankPercent_RNA"]+G4[,"RankPercent_HTA"]),]
G5<-G5[order(G5[,"RankPercent_RNA"]+G5[,"RankPercent_HTA"]),]
G1H<-G1H[order(G1H[,"Splicing Pvalue"]),]
G2H<-G2H[order(G2H[,"Pvalue_HTA"]),]
G3H<-G3H[order(G3H[,"Pvalue_HTA"]),]
G1R<-G1R[order(G1R[,"Pvalue"]),]
G2R<-G2R[order(G2R[,"Pvalue_RNA"]),]
G3R<-G3R[order(G3R[,"Pvalue_RNA"]),]


