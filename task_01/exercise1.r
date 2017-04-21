#loaddata
library(paprbag)
library(randomForest)
data('trainigData')

set.seed(42)  # set seed to get always the same random numbers

postive_samples_positions<-sample.int(100,50)       #get 50 randomly chosen postiv labled datapoints
negative_samples_positions<-sample.int(100,50)+100  #get 50 randomly chosen negativ labled datapoints

trainingselection<- rep(F,200)  # Creats a vector with 200 FALSE entries
trainingselection[postive_samples_positions]<-TRUE
trainingselection[negative_samples_positions]<-TRUE

data_training<-trainingData[trainingselection,]
data_test<-trainingData[!trainingselection,]




rf <- randomForest(label ~ ., data=train)