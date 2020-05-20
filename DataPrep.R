# This is the script to prepare for the training set and testing set for signature gene prediction model
# Please note that this script used pre-calculated gene statistics using functions in SINCERA
# Author: Shuyang Zhao (shuyang.zhao@cchmc.org)

source("https://raw.githubusercontent.com/xu-lab/SINCERA/umi/R/functions.R")
library(plyr)
library(data.table)
library(stringr)
library(tidyverse)
#Cell type to be predicted
db <- "AT1"
celltype <- db
#load pre-calculated gene peak file
load(file="../../peaks/peaks.rda")
#load pre-calculated gene statistics(gs) files
path <- getwd()
list.data<-list()
file.names <- dir(path, pattern =".gs")
for(i in 1:length(file.names)){
  list.data[[i]] <- read.table(file.names[i], header=T, sep="\t")
  rownames(list.data[[i]]) <- list.data[[i]]$gene
}
names(list.data) <- unlist(lapply(file.names, extract_field, field=c(1,2,6), delim="."))
print(names(list.data))

#define negative training gene filtering thresholds for Drop-seq samples
pv.d <- 0.05
fc.d <- 2.5
avg.d <- -Inf
cnt.d <- 100
pct.d <- -Inf
recall.d <- 0.95
#define negative training gene filtering thresholds for C1-fluidigm samples
pv.c <- 0.05
fc.c <- 4
avg.c <- -Inf
cnt.c <- -Inf
pct.c <- -Inf
recall.c <- 1

#gene filtering
list.filter <- list()
list.neggenes <- list()
gs.list <- grep("gs",names(list.data),value=T)
for(j in 1:length(gs.list))
{
  tmp <- substr(gs.list[j],1,2)
  if(tmp == "DP"){
    list.filter[[j]] <- GetSigs(list.data[[gs.list[j]]],groups = substr(colnames(list.data[[gs.list[j]]])[2:3],7,max(nchar(colnames(list.data[[gs.list[j]]])[2:3]))),
                                criteria=c("binom","efsize", "ident.avg","ident.cnt","ident.pct","ident.recall"),
                                thresh=c(pv.d,fc.d,avg.d,cnt.d,pct.d,recall.d), 
                                op=c(-2, 1, 1, 1, 1, 1), do.fdr = F)
  }
  if(tmp == "C1"){
    list.filter[[j]] <- GetSigs(list.data[[gs.list[j]]],groups = substr(colnames(list.data[[gs.list[j]]])[2:3],7,max(nchar(colnames(list.data[[gs.list[j]]])[2:3]))),
                                criteria = c("welch","fc","ident.cnt","ident.pct","ident.recall"),
                                thresh=c(pv.c,fc.c,cnt.c,pct.c,recall.c), 
                                op=c(-2, 1, 1, 1, 1), do.fdr = F)
  }else{}
  list.neggenes[[j]] <- as.character(list.filter[[j]]$gene[which(list.filter[[j]]$group == "Other")])
}
names(list.neggenes) <- unlist(lapply(gs.list, extract_field, field=c(1,2), delim="."))

#find common negative training genes
neg.genes <- Reduce(intersect, list.neggenes)
if(length(neg.genes) < 10){
  cat(paste(length(neg.genes),"genes are selecting as negative genes",sep=" "))
  print(neg.genes)
}else{
  print("10 or more genes are selecting as negative genes")
}

#known celltype markers
#pre-stored markers file
markers <- read.table("Z:/Bioinformatics/Shuyang/Tools/known.celltype.markers.txt",sep="\t",header = T)
pos.genes <- as.character(markers$Gene.Mouse[which(markers$Subtype == db)])
if(length(pos.genes) > 10){
  #check marker qualities in testing set
  tmp <- data.frame(rep(NA,length(pos.genes)))
  rownames(tmp) <- pos.genes
  sig.list <- grep("sig",names(list.data),value=T)
  for(z in 1:length(sig.list)){
    for(v in 1:length(pos.genes)){
      a <- grep(pattern = pos.genes[v],x = as.character(list.data[[sig.list[z]]]$gene),value = T)
      if(length(a) == 0){
        tmp[v,z] <- "NO"
      }else{
        tmp[v,z] <-"YES"
      }
    }
    colnames(tmp)[z] <- sig.list[[z]]
  }  
  tmp$score <- rep(0,nrow(tmp))
  tmp$score <- rowSums(tmp == "YES")
  pos.genes <- pos.genes[which(tmp$score > length(sig.list) * 0.6)]
  #write.table(tmp, file="tmp.txt",sep="\t")
}else{
  print("All markers will be used for training")
}
#traning genes
pos.genes <- as.data.frame(pos.genes)
pos.genes$group <- rep("positive",nrow(pos.genes))
colnames(pos.genes) <- c("genes","group")
neg.genes <- as.data.frame(neg.genes)
neg.genes$group <- rep("negative",nrow(neg.genes))
colnames(neg.genes) <- c("genes","group")
genes.train <- rbind(pos.genes,neg.genes)
rownames(genes.train) <- genes.train$genes
View(genes.train)
#store gene information
write.table(genes.train,file="genes.train.txt",sep="\t",quote=F)


#######################################################################################
################################# TRAINING SET ########################################
#######################################################################################
list.train <- list()
for(j in 1:length(gs.list)){
  tmp <- substr(gs.list[j],1,2)
  if(tmp == "DP"){
    list.train[[j]] <- list.data[[gs.list[j]]][rownames(genes.train),c(1:5,10:13)]
    rownames(list.train[[j]]) <- paste(unlist(lapply(gs.list[j], extract_field, field=c(1,2), delim=".")),rownames(list.train[[j]]),sep=".")
    colnames(list.train[[j]]) <- c("gene","pv.tg","pv.other","fc.tg","fc.other","pct.tg","pct.other","recall.tg","recall.other")
  }
  if(tmp == "C1"){
    list.train[[j]] <- list.data[[gs.list[j]]][rownames(genes.train),c(1:5,8:11)]
    rownames(list.train[[j]]) <- paste(unlist(lapply(gs.list[j], extract_field, field=c(1,2), delim=".")),rownames(list.train[[j]]),sep=".")
    colnames(list.train[[j]]) <- c("gene","pv.tg","pv.other","fc.tg","fc.other","pct.tg","pct.other","recall.tg","recall.other")
  }else{}
}
#combine all data.frames into one data frame
trainset <- do.call(rbind,unname(list.train))
#PEAKS
trainset$peak.1 <- peaks[rownames(trainset),]$Peak1
trainset$peak.2 <- peaks[rownames(trainset),]$Peak2
trainset$peak.1 <- as.character(trainset$peak.1)
trainset$peak.2 <- as.character(trainset$peak.2)
trainset$peak <- rep(0,nrow(trainset))
#assign value
for(i in 1:nrow(trainset)){
  if(trainset[i,]$peak.1 == db & trainset[i,]$peak.2 == db ){
    trainset[i,]$peak <- 1
  }else if(trainset[i,]$peak.1 == db & trainset[i,]$peak.2 != db){
    trainset[i,]$peak <- 0.75
  }else if(trainset[i,]$peak.1 != db & trainset[i,]$peak.2 == db){
    trainset[i,]$peak <- 0.25
  }else{
    trainset[i,]$peak <- 0
  }
}
trainset <-trainset[,-which(colnames(trainset) %in% c("peak.1","peak.2"))]
trainset$peak <- as.factor(trainset$peak)
#create response variables
trainset$y <- rep(0,nrow(trainset))
trainset[which(as.character(trainset$gene) %in% as.character(genes.train$genes[which(genes.train$group == "positive")])),]$y <- 1
trainset$y <- as.factor(trainset$y)
#output training set
write.table(trainset,file="trainset.txt",sep="\t")


#######################################################################################
##################################### TESTING SET #####################################
#######################################################################################
#construct testing set
#filter genes expressed at least 10% in the given celltype across all samples
#positive test genes
#DropSeq
C1.idx <- which(sapply(strsplit(gs.list,"[.]"), `[`, 1) == "C1")
DP.idx <- which(sapply(strsplit(gs.list,"[.]"), `[`, 1) == "DP")
#Drop-seq
gene.list.DP <- list()
for(i in DP.idx){
  tmp <- list.data[[gs.list[i]]]
  genes <- as.character(tmp$gene[which(tmp$ident.pct.AT1 >= 0.1)])
  gene.list.DP[[i]] <- genes
  names(gene.list.DP)[i] <- gs.list[i]
}
gene.list.DP <- gene.list.DP[lapply(gene.list.DP,length)>0]
#C1-fluidigm
gene.list.C1 <- list()
for(i in C1.idx){
  tmp <- list.data[[gs.list[i]]]
  genes <- as.character(tmp$gene[which(tmp$ident.pct.AT1 >= 0.1)])
  gene.list.C1[[i]] <- genes
  names(gene.list.C1)[i] <- gs.list[i]
}
#find intersect
pos.test.genes.DP <- Reduce(intersect, gene.list.DP)
pos.test.genes.C1 <- Reduce(intersect, gene.list.C1)
pos.test.genes <- intersect(pos.test.genes.DP,pos.test.genes.C1)

#negative test genes
#Drop-seq
gene.list.DP <- list()
for(i in DP.idx){
  tmp <- list.data[[gs.list[i]]]
  genes <- as.character(tmp$gene[which(tmp$ident.recall.Other >= 0.95 & tmp$ident.cnt.Other > 20)])
  gene.list.DP[[i]] <- genes
  names(gene.list.DP)[i] <- gs.list[i]
}
gene.list.DP <- gene.list.DP[lapply(gene.list.DP,length)>0]
#C1-fluidigm
gene.list.C1 <- list()
for(i in C1.idx){
  tmp <- list.data[[gs.list[i]]]
  genes <- as.character(tmp$gene[which(tmp$ident.recall.Other >= 0.95)])
  gene.list.C1[[i]] <- genes
  names(gene.list.C1)[i] <- gs.list[i]
}
#find intersect
neg.test.genes.DP <- Reduce(intersect, gene.list.DP)
neg.test.genes.C1 <- Reduce(intersect, gene.list.C1)
neg.test.genes <- intersect(neg.test.genes.DP,neg.test.genes.C1)

#check traning genes qualities
#extract p values
#positive
tmp.p <- matrix(0, nrow = length(pos.test.genes))
tmp.p <- as.data.frame(tmp.p)
row.names(tmp.p) <- pos.test.genes
#negative
tmp.n <- matrix(0, nrow = length(neg.test.genes))
tmp.n <- as.data.frame(tmp.n)
row.names(tmp.n) <- neg.test.genes

for(i in 1:length(gs.list)){
  tmp <- list.data[[gs.list[i]]][pos.test.genes,2] #Given celltype P-values
  tmp.1 <- list.data[[gs.list[i]]][neg.test.genes,3]  #Other celltype P-values
  
  tmp.p <- cbind(tmp.p,tmp)
  tmp.n <- cbind(tmp.n,tmp.1)
  colnames(tmp.p)[i+1] <- gs.list[i]
  colnames(tmp.n)[i+1] <- gs.list[i]
  
}
tmp.p <- tmp.p[,-1]
tmp.n <- tmp.n[,-1]

#check significance
#positive
tmp.p$x <- rep(0,nrow(tmp.p))
tmp.p$y <- rep(0,nrow(tmp.p))
for(i in 1:nrow(tmp.p)){
  tmp.p[i,]$x <- as.numeric(length(which(tmp.p[i,c(1:length(gs.list))] < 0.1)))
  tmp.p[i,]$y <- as.numeric(tmp.p[i,]$x/length(gs.list))
}
#positive testing genes have to show sigfificance(p < 0.1) for targeted cell type in at least 80% of the samples
pos.test.genes <- rownames(tmp.p)[which(tmp.p$y > 0.8)]

#negative testing genes have to show significance(p < 0.05) for all other cell types in all samples
tmp.n$x <- rep(0,nrow(tmp.n))
tmp.n$y <- rep(0,nrow(tmp.n))
for(i in 1:nrow(tmp.n)){
  tmp.n[i,]$x <- as.numeric(length(which(tmp.n[i,c(1:length(gs.list))] < 0.05)))
  tmp.n[i,]$y <- as.numeric(tmp.n[i,]$x/length(gs.list))
}
neg.test.genes <- rownames(tmp.n)[which(tmp.n$y == 1)]
#output
genes.test <- c(pos.test.genes,neg.test.genes)
write.table(genes.test,file="genes.test.txt",sep="\t",quote = F)

#extract testset for each time points using this set of filtered genes
list.test <- list()
for(i in 1:length(gs.list)){
  if(substr(gs.list[i],1,2) =="DP"){
    list.test[[i]] <- list.data[[gs.list[i]]][genes.test,c(1:5,10:13)]
    rownames(list.test[[i]]) <- paste(unlist(lapply(gs.list[i], extract_field, field=c(1,2), delim=".")),genes.test,sep=".")
    colnames(list.test[[i]]) <- c("gene","pv.tg","pv.other","fc.tg","fc.other","pct.tg","pct.other","recall.tg","recall.other")
  }
  if(substr(gs.list[i],1,2) =="C1"){
    list.test[[i]] <- list.data[[gs.list[i]]][genes.test,c(1:5,8:11)]
    rownames(list.test[[i]]) <- paste(unlist(lapply(gs.list[i], extract_field, field=c(1,2), delim=".")),genes.test,sep=".")
    colnames(list.test[[i]]) <- c("gene","pv.tg","pv.other","fc.tg","fc.other","pct.tg","pct.other","recall.tg","recall.other")
  }else{}
  names(list.test)[i] <- gs.list[i]
}

#fill in peak info
for(i in 1:length(gs.list)){
  list.test[[gs.list[i]]]$peak.1 <- peaks[rownames(list.test[[gs.list[i]]]),]$Peak1
  list.test[[gs.list[i]]]$peak.2 <- peaks[rownames(list.test[[gs.list[i]]]),]$Peak2
  list.test[[gs.list[i]]]$peak.1 <- as.character(list.test[[gs.list[i]]]$peak.1)
  list.test[[gs.list[i]]]$peak.2 <- as.character(list.test[[gs.list[i]]]$peak.2)
  list.test[[gs.list[i]]]$peak <- rep(0,nrow(list.test[[gs.list[i]]])) 
  for(j in 1:nrow(list.test[[gs.list[i]]])){
    if(list.test[[gs.list[i]]][j,]$peak.1 == db & list.test[[gs.list[i]]][j,]$peak.2 == db ){
      list.test[[gs.list[i]]][j,]$peak <- 1
    }else if(list.test[[gs.list[i]]][j,]$peak.1 == db & list.test[[gs.list[i]]][j,]$peak.2 != db){
      list.test[[gs.list[i]]][j,]$peak <- 0.75
    }else if(list.test[[gs.list[i]]][j,]$peak.1 != db & list.test[[gs.list[i]]][j,]$peak.2 == db){
      list.test[[gs.list[i]]][j,]$peak <- 0.25
    }else{
      list.test[[gs.list[i]]][j,]$peak <- 0
    }
  }
  list.test[[gs.list[i]]] <- list.test[[gs.list[i]]][,-which(colnames(list.test[[gs.list[i]]])%in%c("peak.1","peak.2"))]
  list.test[[gs.list[i]]]$peak <- as.factor(list.test[[gs.list[i]]]$peak)
  #write.table(list.test[[gs.list[i]]],file=paste(gs.list[i],".test.txt",sep=""),sep="\t")
  
}

######### save all prepared datasets ###########
save(trainset,list.test,genes.train,file=paste(db,".model.data.package.Rda",sep=""))

