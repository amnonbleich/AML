#install.packages("kernlab")
#install.packages("seqinr")
setwd("./Uni/Master/aml/AML/task_05/")
library(seqinr)# for loading fasta
#library(kernlab)
#source("https://bioconductor.org/biocLite.R")
#biocLite("kebabs")
library(kebabs)
library(e1071)
library(ROCR)
library(caret)

stripe_upper<-function(element)
{

return(gsub('[[:lower:]]','',element[[1]]))
}

strip_all<<-function(element)
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





####################
#     Task 2       #
####################

#stringkernel <- stringdot(type="spectrum", length=2, normalized=TRUE)
#
#sk_model <- ksvm(lables~.,data=training_set,kernel=stringkernel)
#test <- predict()
specK2 <- spectrumKernel(k=3)




# k-fold 
set.seed(2)
folds <- createFolds(dataset_wo_flanks$seq, k = 10)
all.perf<-c()

for (cost in c(0.5,1,2,3,4,5,6))
{
  cat(cost)
  lables<- c()
  predictions<- c()
  for (fold in folds)
  {
    model <- kbsvm(x=dataset_wo_flanks$seq[-fold],y=dataset_wo_flanks$lables[-fold],kernel=specK2,pkg="e1071",  svm="C-svc",  cost=cost)
    predictions<- c(predictions,predict(model, dataset_wo_flanks$seq[fold]))
    lables<- c(lables,dataset_wo_flanks$lables[fold])
  }
  all.perf<-c(all.perf,prediction(as.factor(as.logical(predictions)),lables))
  
}







# k-fold 

all.perfK<-c()

for (k in c(2,3,4,5))
{
  specK <- spectrumKernel(k=k)
  cat(k)
  lables<- c()
  predictions<- c()
  for (fold in folds)
  {
    model <- kbsvm(x=dataset_wo_flanks$seq[-fold],y=dataset_wo_flanks$lables[-fold],kernel=specK,pkg="e1071",  svm="C-svc")
    predictions<- c(predictions,as.logical(predict(model, dataset_wo_flanks$seq[fold])))
    lables<- c(lables,as.logical(dataset_wo_flanks$lables[fold]))
  }
  browser()
  all.perfK<-c(all.perfK,prediction(as.factor(as.logical(predictions)),lables))
  
}





