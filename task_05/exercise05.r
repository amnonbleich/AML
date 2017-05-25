#####################################
# Exercise 5 - AML - 26.05.17       #
# Ben Wulf, Amnon Bleich, Lie Hong  #
#####################################


#install.packages("seqinr")
setwd("./Uni/Master/aml/AML/task_05/")
library(seqinr)# for loading fasta
#source("https://bioconductor.org/biocLite.R")
#biocLite("kebabs")
library(kebabs)
library(e1071)


####################
#     Task 1       #
####################


# We decided to use the  PUM2 dataset and the kebabs package. kernlab crashes caused by allocation errors


###
# This function removes all lower capitals from a given string in a list
###
stripe_upper<-function(element)
{
return(gsub('[[:lower:]]','',element[[1]]))
}


###
# This function give a string in a list back as lower letter string
###

stripe_all<-function(element)
{
return(tolower(element[[1]]))
}


###
# Load the positive sequences as vector of strings
###
PUM2pos <- read.fasta('./rna-binding/positive_PUM2.fasta',as.string=T,forceDNAtolower=F)
PUM2pos_striped <-sapply(PUM2pos,stripe_upper)  # only Upper letters
PUM2pos_all <- sapply(PUM2pos,stripe_all)       # full string aka binding site with flanks



###
# Load the negative sequences as vector of strings
###
PUM2neg <- read.fasta('./rna-binding/negative_PUM2.fasta', as.string=T,forceDNAtolower=F)
PUM2neg_striped <- sapply(PUM2neg,stripe_upper) # only Upper letters
PUM2neg_all <- sapply(PUM2neg,stripe_all)       # full string aka binding site with flanks

###
# create a dataframe containing the lable and the binding site sequence. Sequence is a DNAStringSet as required for the ksvm.
###
dataset_wo_flanks<-  data.frame( lables=as.factor(c(rep(T,length(PUM2pos_striped)),rep(F,length(PUM2neg_striped)))))
dataset_wo_flanks$seq<- DNAStringSet(c(PUM2pos_striped,PUM2neg_striped))


###
# create a dataframe containing the lable and the binding-site-sequence with the flanks. Sequence is a DNAStringSet as required for the ksvm.
###
dataset_flanks<-  data.frame( lables=as.factor(c(rep(T,length(PUM2pos_striped)),rep(F,length(PUM2neg_striped)))))
dataset_flanks$seq<- DNAStringSet(c(PUM2pos_all,PUM2neg_all))


####################
#     Task 2       #
####################

# We use the kebabs package

specK <- spectrumKernel(k=3) # create a spectrum kernel with kmer size 3

###
# Train the SVM with a subset of datapoints kmer size 3 and default costs
###

basic_model <- kbsvm( 
  x=dataset_wo_flanks$seq[-(10000:11692)],
  y=dataset_wo_flanks$lables[-(10000:11692)],
  kernel=specK,
  pkg="e1071",
  svm="C-svc")

# predict the rest of the datapoints with the trained svm
pred<- predict(model_init,dataset_wo_flanks$seq[(10000:11692)])



#############################################################
#     Task 3                                                #
# Tune your model: select the cost C and the k-mer length   #
# that produce the best classi  fication accuracy (or AUC). #
#############################################################


###
# with the ksvm it is possivle to perform a crossvalitation automaticly and to get the accuracy and the area under the curve
# we perform a 10-fold cv and reciving the auccuracy and the auc
###

# ACC
acc<-list() # stores the accuracy for the different k
acc_min<-3  # just for axis limits for the plot
acc_max<-0  # just for axis limits for the plot

# AUC
auc<-list() # stores the area under the curve for the different k
auc_min<-3  # just for axis limits for the plot
auc_max<-0  # just for axis limits for the plot


ks<-c(1,2,3,4,5)  # the different kmer sizes we tried out
costs<-c(0.5,0.75,1,2,4,8,10) # the different costs we use

for (k in ks) # for each kmer size
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
    accC<-c(accC,model@cvResult@ACC)  # save the accuracy value from 10-fold cv
    aucC<-c(aucC,model@cvResult@AUC)  # save the area under the curve value from 10-fold cv
  }
  names(accC)<-costs
  names(aucC)<-costs
  auc[[k]]<-aucC
  acc[[k]]<-accC
  
  # just plotting limits
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




####################
#     Task 5       #
####################


## not really what was asked... see lecture slide 26

model@featureWeights







