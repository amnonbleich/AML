#loaddata
library(paprbag)
library(randomForest)
data('trainingData')

confusion <-function (data, reference, mute=F)
{
  equal <- data==reference
  TP<- sum(data[equal])
  TN<- sum(!data[equal])
  FP<- sum(data[!equal])
  FN<- sum(!data[!equal])
  
  
}



set.seed(42)  # set seed to get always the same random numbers

postive_samples_positions<-sample.int(100,50)       #get 50 randomly chosen postiv labled datapoints
negative_samples_positions<-sample.int(100,50)+100  #get 50 randomly chosen negativ labled datapoints

trainingselection<- rep(F,200)  # Creats a vector with 200 FALSE entries
trainingselection[postive_samples_positions]<-TRUE
trainingselection[negative_samples_positions]<-TRUE

data_training<-trainingData[trainingselection,]


data_test<-trainingData[!trainingselection,]

rf_model <- randomForest(Labels ~ ., data=data_training,xtest=data_test[,-1],ytest=data_test[,1])

print_results <- function (x)
{
  cat("Training Set\n")
  cat(paste("Training set error:",mean(x$confusion[,3],'\n')))
  print(x$confusion)

  cat("Test Set")
  cat(paste("Test set error:",mean(x$test$confusion[,3])))
  print(x$confusion)
}

print_results(rf_model)

print(paste("Sensitivity Train set:",rf_model ))

model_test <- predict(rf_model,data=data_test)


