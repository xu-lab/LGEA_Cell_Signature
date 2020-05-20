# This is the script for predicting and ranking 
# Author: Shuyang Zhao (shuyang.zhao@cchmc.org)

source("https://raw.githubusercontent.com/xu-lab/SINCERA/umi/R/functions.R")
library(plyr)
library(dplyr)
library(Matrix)
library(reshape2)
library(openxlsx)
library(glmnet)
library(caret)
library(ROCR)
############################ Loading prepared data package ######################################
celltype <- "AT1"
#load pre-stored traiing and testing set
load(file=paste("GS.UMI/",celltype,"/",celltype,".model.data.package.Rda",sep=""))
#variable transformation
trainset$peak <- as.numeric(as.character(trainset$peak))
for(i in 1:length(list.test)){
  list.test[[i]]$peak <- as.numeric(as.character(list.test[[i]]$peak))
 }
dim(trainset)
class(trainset$peak)

##############################################################################################################
################################################  model fiting ###############################################
##############################################################################################################
trainset <- trainset[,-1]
#predictor variables
x <-  model.matrix(y~.,trainset)[,-1]
#response variable
y <- trainset$y

#############ELASTIC NET################
#set.seed(12345)
#"caret" train
elastic <- train(
  y ~., data = trainset, method = "glmnet",
  trControl = trainControl("cv", number = 10),
  tuneLength = 10
)
#model coefficients
elastic$bestTune$alpha
elastic$bestTune$lambda
#fitting the model
coef(elastic$finalModel, elastic$bestTune$lambda)
#using alpha, lambda value from "caret" and fit glmnet() model
elnet.model <- glmnet(x,y,family="binomial",
                 alpha=elastic$bestTune$alpha,
                 lambda =elastic$bestTune$lambda)

############################################################################################
#Looping throught testing sets
list.predict <- list()
for(z in 1:length(list.test)){
  testset <- list.test[[z]]

#get test data
testset <- testset[,-1]
testset$y <- rep(0,nrow(testset))
x_test <- model.matrix(y~.,testset)[,-1]
#dim(x_test)

#predict probability
elnet_prob <- predict(elnet.model,newx = x_test,type="response")
elnet_prob <- as.data.frame(elnet_prob)
genes <- unlist(lapply(rownames(elnet_prob), extract_field, field=3, delim="."))
elnet_prob$gene <- genes

#loop through all alphas and lambdas parameters
grid<-expand.grid(.alpha=seq(0,1,by=.1),.lambda=seq(0,1,by=.1))

for (i in unique(grid$.alpha)){
  for (j in unique(grid$.lambda)){
    
    model <- glmnet(x,y,family="binomial",
                          alpha=i,
                          lambda =j)
    prob <- predict(model,newx = x_test,type="response")
    prob <- as.data.frame(prob)
    col.name <- paste("alpha",i, "lambda",j, sep=".")
    elnet_prob[,col.name] <- prob$s0
    
  }
}
colnames(elnet_prob)[1] <- "alpha.best.lambda.best"
#save(elnet_prob,file=paste(celltype,"elnet_prob.rda",sep="."))

#ROC to test parameters performance
ROC.genes <- rownames(genes.train)[which(genes.train$group == "positive")]
index <- which(elnet_prob$gene %in% ROC.genes)
length(index)
elnet_prob.1 <- elnet_prob[,-2]
#rownames(elnet_prob.1)[index]

auc.perf <- data.frame(colnames(elnet_prob.1),auc = 0)
colnames(auc.perf) <- c("parameters","auc")
for (i in 1:ncol(elnet_prob.1)){
i.cells <- data.frame(p = elnet_prob.1[,i], s=0)
#get marker genes
i.cells[index,]$s <- 1
#AUC
i.pred <- prediction(i.cells$p, i.cells$s)
i.perf <- performance(i.pred,"auc")
i.auc <- i.perf@y.values[[1]]
auc.perf[i,]$auc <- i.auc
}
#locate the parameters with the best result
param.select <- as.character(auc.perf[which(auc.perf$auc == max(auc.perf$auc)),]$parameters)

###################################################################################################################
##################################################### RANKING #####################################################
###################################################################################################################
#print(param.select)
if(length(param.select) >1){
  print(paste("Converting multiple parameters for",names(list.test)[z],sep=" "))
  predict.best <- elnet_prob[,c("gene",param.select)]
  tmp.rank <- matrix(0,nrow=nrow(predict.best),ncol=length(param.select))
  rownames(tmp.rank) <- genes
  tmp.rank <- as.data.frame(tmp.rank)
  colnames(tmp.rank) <- param.select
  for(j in 2:ncol(predict.best)){
    predict.best <- predict.best[order(predict.best[,j],decreasing = T),]
    rank <- as.data.frame(as.matrix(cbind(predict.best$gene,c(1:nrow(predict.best)))))
    rownames(rank) <- rank$V1
    rank <- rank[rownames(tmp.rank),]
    tmp.rank[,j-1] <- as.numeric(as.character(rank$V2))
  }
  tmp.rank$avg.rank <- rowSums(tmp.rank)/length(param.select)
  tmp.rank <- tmp.rank[order(tmp.rank$avg.rank,decreasing = F),]
  tmp.rank$rank <- c(1:nrow(tmp.rank))
  tmp.rank$gene <- rownames(tmp.rank)
  tmp.rank$source <- rep(unlist(lapply(names(list.test)[z], extract_field, field=1, delim=".")),nrow(tmp.rank))
  tmp.rank$tps <- rep(unlist(lapply(names(list.test)[z], extract_field, field=2, delim=".")),nrow(tmp.rank))
  tmp.rank <- tmp.rank[,-which(colnames(tmp.rank)%in%param.select)]
  tmp.rank <- tmp.rank[,c(3,4,5,1,2)]
  rownames(tmp.rank) <- paste(tmp.rank$source,tmp.rank$tps,tmp.rank$gene)
  list.predict[[z]] <- tmp.rank
  names(list.predict)[z] <- names(list.test)[z]
  write.table(tmp.rank,paste(celltype,names(list.test)[z],"predict.best.txt",sep="."), sep="\t",quote = F)
  }
else{
  print(paste("The best parameters for",names(list.test)[z],"is",param.select,sep=" "))
  predict.best <- elnet_prob[,c("gene",param.select)]
  predict.best$source <- unlist(lapply(rownames(predict.best), extract_field, field=1, delim="."))
  predict.best$tps <- unlist(lapply(rownames(predict.best), extract_field, field=2, delim="."))
  predict.best <- predict.best[,c(1,3,4,2)]  
  colnames(predict.best) <- c("gene","source","tps","prob")
  predict.best <- predict.best[order(predict.best$prob,decreasing = T),]
  predict.best$rank <- c(1:nrow(predict.best))
  list.predict[[z]] <- predict.best
  names(list.predict)[z] <- names(list.test)[z]
  #output individual sample ranking
  #write.table(predict.best,paste(celltype,names(list.test)[z],"predict.best.txt",sep="."), sep="\t",quote = F)
}
  
}

#calculate the average rank
ranks <- matrix(0,nrow=length(genes),ncol=length(list.predict))
rownames(ranks) <- genes
ranks <- as.data.frame(ranks)
for(i in 1:length(list.predict)){
  tmp <- list.predict[[i]]
  rownames(tmp) <- as.character(tmp$gene)
  a <- tmp[rownames(ranks),which(colnames(tmp)=="rank")]
  ranks[,i] <- a
  colnames(ranks)[i] <- names(list.predict)[i]
}
ranks$sum.rank <- rowSums(ranks)
ranks$avg.rank <- ranks$sum.rank/length(list.test)
#order based on average ranks
ranks <- ranks[order(ranks$avg.rank,decreasing = F),]
#output final ranking 
write.table(ranks,paste(celltype,"ranks.txt",sep="."), sep="\t",quote = F)

