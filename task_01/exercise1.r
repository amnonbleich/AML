#loaddata
library(paprbag)
library(randomForest)
data('trainingData')




set.seed(42)  # set seed to get always the same random numbers


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
rf_model <- randomForest(Labels ~ ., data=data_training,xtest=data_test[,-1],ytest=data_test[,1])

cat('Training-Set Confusion matrix')
print(rf_model$confusion)
cat(paste("sensitivity train-set:",rf_model$confusion[2,2]/(rf_model$confusion[2,2]+rf_model$confusion[2,1]) ),'\n')
cat(paste("specifity train-set:",rf_model$confusion[1,1]/(rf_model$confusion[1,1]+rf_model$confusion[1,2]) ),'\n')



cat('Test-Set Confusion matrix')
print(rf_model$test$confusion)

cat(paste("sensitivity test-set:",rf_model$test$confusion[2,2]/(rf_model$test$confusion[2,2]+rf_model$test$confusion[2,1]) ),'\n')
cat(paste("specifity test-set:",rf_model$test$confusion[1,1]/(rf_model$test$confusion[1,1]+rf_model$test$confusion[1,2]) ),'\n')






# Apply predict manualy 
model_test <- predict(rf_model,data_test)


# Task 3

LOOCV<- function(formula,data,indices)
{
  print (indices)
  subdata<- data[indices,]   # allow boot to select rows
  testresults<-matrix(c(0,0),ncol=2,nrow=2)
  
  for (i in 1:(nrow(data)))
    {
    print(paste("iteration",i))
    rf_model <- randomForest(Labels ~ ., data=subdata[-i,])
    model_test <- as.logical(predict(rf_model,subdata[i,-1]))

    
     matrix_position=2*as.logical(subdata[i,1])+model_test+1   # if reference TRUE select second row, otherwise first. if test TRUE select Second column otherwise first. + R offset for legal intervall 1:4
    #print(as.logical(data[i,1]))
    #print(model_test)
    
    #          P0  P1
    # Ground 0 TN  FP
    # Ground 0 FN  TP
    
    #print(matrix_position)
    testresults[matrix_position]<-testresults[matrix_position]+1
    #print(testresults)
  }

  sens<-testresults[4]/(testresults[4]+testresults[3])
  spec<-testresults[1]/(testresults[1]+testresults[1])
  
  
  return(c(sens,spec))
  
}

                                        #subset shut be removed
bootstrapObj<- boot(data = trainingData[80:120,], statistic=LOOCV, R=1000, formula=Labels ~ .)

# Get confidence intervalls
sen<-boot.ci(results, index=1)
spe<-boot.ci(results, index=2)
