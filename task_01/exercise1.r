# aml -Exercise 01
#Lie Hong, Amnon Bleich, Ben Wulf
# Task 1
library(paprbag)
library(randomForest)
library(boot)
data('trainingData')

set.seed(42)  # set seed to get always the same random numbers

#Task 2
#Choose some random values.
postive_samples_positions<-sample.int(100,50)       #get 50 randomly chosen postiv labled datapoints
negative_samples_positions<-sample.int(100,50)+100  #get 50 randomly chosen negativ labled datapoints


trainingselection<- rep(F,200)  # Creats a vector with 200 FALSE entries
trainingselection[postive_samples_positions]<-TRUE
trainingselection[negative_samples_positions]<-TRUE

# Select the datapoints for training
data_training<-trainingData[trainingselection,]
# select the datapoints for testing
data_test<-trainingData[!trainingselection,]

# Train the randomForest model and test it automaticly 
rf_model <- randomForest(Labels ~ ., data=data_training)

# Get the confusion-matrix for the training dataset and calculate sensitivity an specifity.
cat('Training-Set confusion-matrix')
print(rf_model$confusion)
cat(paste("sensitivity train-set:",rf_model$confusion[2,2]/(rf_model$confusion[2,2]+rf_model$confusion[2,1]) ),'\n')
cat(paste("specifity train-set:",rf_model$confusion[1,1]/(rf_model$confusion[1,1]+rf_model$confusion[1,2]) ),'\n')


# Apply predict manualy 
model_test <- as.logical(predict(rf_model,data_test[,-1]))
TP<-sum(model_test==T &model_test==as.logical(data_test[,1]))
FP<-sum(model_test==T &model_test!=as.logical(data_test[,1]))
FN<-sum(model_test==F &model_test!=as.logical(data_test[,1]))
TN<-sum(model_test==F &model_test==as.logical(data_test[,1]))
conf<- matrix(c(TP,FN,FP,TN),ncol=2)
colnames(conf)<-c('pred=T','pred=F')
rownames(conf)<-c('ground=T','ground=F')
print (conf)
cat(paste("sensitivity test-set:",TP/(TP+FN),'\n'))
cat(paste("specifity test-set:",TN/(TN+FP) ),'\n')


# Task 3

LOOCV<- function(data,indices)
{
  
  subdata<- data[indices,]   # allow boot to select rows
  testresults<-matrix(c(0,0),ncol=2,nrow=2) # to store the results
  
  for (i in 1:(nrow(data))) # for each data point
    {
    rf_model <- randomForest(subdata[-i,-1],subdata[-i,1])  # remove it from the training Dataset
    model_test <- as.logical(predict(rf_model,subdata[i,-1])) # and predict for the single Datapoint

    
    matrix_position=2*as.logical(subdata[i,1])+model_test+1   # if reference TRUE select second row, otherwise first. if test TRUE select Second column otherwise first. + R offset for legal intervall 1:4

    
    #          P0  P1
    # Ground 0 TN  FP
    # Ground 0 FN  TP
    
    testresults[matrix_position]<-testresults[matrix_position]+1
  }

  sens<-testresults[4]/(testresults[4]+testresults[3]) # calculate sensitivtiy
  spec<-testresults[1]/(testresults[1]+testresults[2]) # calculate specifity
  
  cat('.')
  return(c(sens,spec))
  
}


bootstrapObj<- boot(data = trainingData[91:110], statistic=LOOCV, R=10)

# try to get confidence intervalls
sen<-boot.ci(bootstrapObj, index=1)
spe<-boot.ci(bootstrapObj, index=2)
