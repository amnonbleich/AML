#install.packages("kernlab")
#install.packages("seqinr")
setwd("./Uni/Master/aml/AML/task_05/")
library(seqinr)# for loading fasta
#library(kernlab)
#source("https://bioconductor.org/biocLite.R")
#biocLite("kebabs")
library(kebabs)
library(e1071)
#library(ROCR)
#library(caret)

stripe_upper<-function(element)
{

return(gsub('[[:lower:]]','',element[[1]]))
}

stripe_all<-function(element)
{
  
  return(tolower(element[[1]]))
}


PUM2pos <- read.fasta('./rna-binding/positive_PUM2.fasta',as.string=F,forceDNAtolower=F)
PUM2pos_striped <-sapply(PUM2pos,stripe_upper)
PUM2pos_all <- sapply(PUM2pos,stripe_all)

PUM2neg <- read.fasta('./rna-binding/negative_PUM2.fasta', as.string=T,forceDNAtolower=F)
PUM2neg_striped <- sapply(PUM2neg,stripe_upper)
PUM2neg_all <- sapply(PUM2neg,stripe_all)


#splitpercent<- 0.85

#positive_training_idx <- sample(length(PUM2pos_striped),floor(length(PUM2pos_striped)*splitpercent))
#negative_training_idx <- sample(length(PUM2neg_striped),floor(length(PUM2neg_striped)*splitpercent))

#positive_training_set <- PUM2pos_striped[positive_training_idx]
#positive_test_set     <- PUM2pos_striped[-positive_training_idx]

#negative_training_set <- PUM2neg_striped[negative_training_idx]
#negative_test_set <- PUM2neg_striped[-negative_training_idx]


# finaly create the training and test set as two columned data.frames
#training_set <- data.frame( lables=as.factor(c(rep(T,length(positive_training_set)),rep(F,length(negative_training_set)))))
#training_set$seq<- DNAStringSet(c(positive_training_set,negative_training_set))

#test_set <- data.frame(
#  seq=c(positive_test_set,negative_test_set),
#  lables=as.factor(c(rep(T,length(positive_test_set)),rep(F,length(negative_test_set)))))


dataset_wo_flanks<-  data.frame( lables=as.factor(c(rep(T,length(PUM2pos_striped)),rep(F,length(PUM2neg_striped)))))
dataset_wo_flanks$seq<- DNAStringSet(c(PUM2pos_striped,PUM2neg_striped))
#dataset_wo_flanks<- dataset_wo_flanks[sample(nrow(dataset_wo_flanks),nrow(dataset_wo_flanks)),]


dataset_flanks<-  data.frame( lables=as.factor(c(rep(T,length(PUM2pos_striped)),rep(F,length(PUM2neg_striped)))))
dataset_flanks$seq<- DNAStringSet(c(PUM2pos_all,PUM2neg_all))




####################
#     Task 2 /3    #
####################

#stringkernel <- stringdot(type="spectrum", length=2, normalized=TRUE)
#
#sk_model <- ksvm(lables~.,data=training_set,kernel=stringkernel)
#test <- predict()


acc<-list()
acc_min<-3
acc_max<-0
auc<-list()
auc_min<-3
auc_max<-0


ks<-c(1,2,3,4,5)
costs<-c(0.5,0.75,1,2,4,8,10)

for (k in ks)
{
  accC<-c()
  aucC<-c()
  
  for (cost in costs)
  {
  cat(paste(k,'-',cost,'\n'))
  specK <- spectrumKernel(k=k)
  model <- kbsvm(
    x=dataset_wo_flanks$seq,
    y=dataset_wo_flanks$lables,
    kernel=specK,
    pkg="e1071",
    svm="C-svc",
    perfParameters=c('ACC','AUC'),
    cross=10,
    showProgress=T,
    cost=cost)
    accC<-c(accC,model@cvResult@ACC)
    aucC<-c(aucC,model@cvResult@AUC)
  }
  names(accC)<-costs
  names(aucC)<-costs
  auc[[k]]<-aucC
  acc[[k]]<-accC
  auc_min<-min(auc_min,aucC)
  auc_max<-max(auc_max,aucC)
  acc_min<-min(acc_min,accC)
  acc_max<-max(acc_max,accC)
  
}

## Plot the AUC
plot(costs,auc[[1]],col=rainbow(5)[1],ylim=c(auc_min,auc_max),type='l',xlim=c(0,10),ylab='AUC- Area under the curve', main="AUC's in dependence of costs and kmer size - wo flanks")

for(i in 2:5)
{lines(costs,auc[[i]],col=rainbow(5)[i])}
legend("bottomright",legend=1:5,col=rainbow(5),pch='l')


## Plot the Accuracy
plot(costs,acc[[1]],col=rainbow(5)[1],ylim=c(acc_min,acc_max),type='l',xlim=c(0,10),ylab='Accuracy', main="Accuracy in dependence of costs and kmer size - wo flanks")
for(i in 2:5)
{lines(costs,acc[[i]],col=rainbow(5)[i])}
legend("bottomright",legend=1:5,col=rainbow(5),pch='l')


####################
#     Task 4       #
####################



accf<-list()
accf_min<-3
accf_max<-0
aucf<-list()
aucf_min<-3
aucf_max<-0


ks<-c(1,2,3,4,5)
costs<-c(0.5,0.75,1,2,4,8,10)

for (k in ks)
{
  accC<-c()
  aucC<-c()
  
  for (cost in costs)
  {
    cat(paste(k,'-',cost,'\n'))
    specK <- spectrumKernel(k=k)
    model <- kbsvm(
      x=dataset_flanks$seq,
      y=dataset_flanks$lables,
      kernel=specK,
      pkg="e1071",
      svm="C-svc",
      perfParameters=c('ACC','AUC'),
      cross=10,
      showProgress=T,
      cost=cost)
    accC<-c(accC,model@cvResult@ACC)
    aucC<-c(aucC,model@cvResult@AUC)
  }
  names(accC)<-costs
  names(aucC)<-costs
  aucf[[k]]<-aucC
  accf[[k]]<-accC
  aucf_min<-min(aucf_min,aucC)
  aucf_max<-max(aucf_max,aucC)
  accf_min<-min(accf_min,accC)
  accf_max<-max(accf_max,accC)
  
}

## Plot the AUC
plot(costs,aucf[[1]],col=rainbow(5)[1],ylim=c(aucf_min,aucf_max),type='l',xlim=c(0,10),ylab='AUC- Area under the curve', main="AUC's in dependence of costs and kmer size - with flanks")

for(i in 2:5)
{lines(costs,aucf[[i]],col=rainbow(5)[i])}
legend("bottomright",legend=1:5,col=rainbow(5),pch='l')


## Plot the Accuracy
plot(costs,accf[[1]],col=rainbow(5)[1],ylim=c(accf_min,accf_max),type='l',xlim=c(0,10),ylab='Accuracy', main="Accuracy in dependence of costs and kmer size - with flanks")
for(i in 2:5)
{lines(costs,accf[[i]],col=rainbow(5)[i])}
legend("bottomright",legend=1:5,col=rainbow(5),pch='l')













