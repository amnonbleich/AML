#install.packages("kernlab")
#install.packages("seqinr")
setwd("./Uni/Master/aml/AML/task_05/")
library(seqinr)# for loading fasta
library(kernlab)


?stripe_upper<-function(element)
{

return(gsub('[[:lower:]]','',element[[1]]))
}


PUM2pos <- read.fasta('./rna-binding/positive_PUM2.fasta',as.string=T,forceDNAtolower=F)
PUM2pos_striped <-sapply(PUM2pos,stripe_upper)

PUM2neg <- read.fasta('./rna-binding/negative_PUM2.fasta', as.string=T,forceDNAtolower=F)
PUM2neg_striped <- sapply(PUM2neg,stripe_upper)


splitpercent<- 0.85

positive_training_idx <- sample(length(PUM2pos_striped),floor(length(PUM2pos_striped)*splitpercent))
negative_training_idx <- sample(length(PUM2neg_striped),floor(length(PUM2neg_striped)*splitpercent))

positive_training_set <- PUM2pos_striped[positive_training_idx]
positive_test_set     <- PUM2pos_striped[-positive_training_idx]

negative_training_set <- PUM2neg_striped[negative_training_idx]
negative_test_set <- PUM2neg_striped[-negative_training_idx]


# finaly create the training and test set as two columned data.frames
training_set <- data.frame(
  seq=c(positive_training_set,negative_training_set),
  lables=as.factor(c(rep(T,length(positive_training_set)),rep(F,length(negative_training_set)))))

test_set <- data.frame(
  seq=c(positive_test_set,negative_test_set),
  lables=as.factor(c(rep(T,length(positive_test_set)),rep(F,length(negative_test_set)))))


