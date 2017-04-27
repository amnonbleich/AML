#source("https://bioconductor.org/biocLite.R")
#biocLite() 

#install.packages("http://github.com/crarlus/paprbag")
#install.packages("devtools")
#library("devtools")
#source("https://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
#install.packages("iterators")
#install.packages("ade4")

#devtools::install_github("crarlus/paprbag")

#install.packages("randomForest")

# Step1
library("paprbag")
library("randomForest")

data("trainingData")
ind_pos <- sample(1:100,50,replace = F)
ind_neg <- sample(101:200,50,replace = F)
test <- trainingData[c(ind_pos, ind_neg), ]
training <- trainingData[-c(ind_pos,ind_neg), ]

#seperate label from data
label_test = test[, 1]
training_test = test[, -1]
label_training = training[, 1]
training_training = training[, -1]

#training
forest = randomForest(x = training_training, y = label_training)
#test / predicting 
test_prediction <- predict(forest, training_test)

# transform into boolean vector
test_prediction <- test_prediction == T
label_test <- label_test == T

#confusion matrix values
TP = sum(test_prediction[label_test] == T)
FP = sum(!test_prediction[label_test] == T)
TN = sum(!test_prediction[!label_test] == T)
FN = sum(test_prediction[!label_test] == T)
sensitivity = TP/(TP+FN)
specificity = TN/(TN+FP)
sensitivity
specificity


#Step 2
#Leave one out
cvTP <- 0
cvFP <- 0
cvTN <- 0
cvFN <- 0
labels <- trainingData[,1]
trainingData_no_label <- trainingData[,-1] 

#leave one out
for(i in 1:200){
  if(i %% 10 == 0){
    print(paste(i, "out of 200", sep = " "))
  }
  label_test = labels[i] == T
  training_test = trainingData_no_label[i, ]
  label_training = labels[-i]
  training_training = trainingData_no_label[-i, ]
  
  # tforest <- randomForest(x = trainingData[,-1], y = trainingData[,1])
  forest = randomForest(x = training_training, y = label_training)
  test_prediction <- predict(forest, training_test) == T
  cvTP <- cvTP + as.numeric(test_prediction && label_test)
  cvFP <- cvFP + as.numeric(test_prediction && !label_test)
  cvTN <- cvTN + as.numeric(!test_prediction && !label_test)
  cvFN <- cvFN + as.numeric(!test_prediction && label_test)
}

cvsensitivity = cvTP/(cvTP+cvFN)
cvspecificity = cvTN/(cvTN+cvFP)
cvsensitivity
cvspecificity



#Bootstrap 95% Confidence Interval
bTP <- NULL
bTN <- NULL
#n has to be divisible by 40
n <- 40

for(i in 1:n){
  if(i %% 10 == 0){
    print(paste(i, "out of", n, sep = " "))
  }
  
  ind <- sample(1:200,100,replace = F)
  test <- trainingData[ind, ]
  training <- trainingData[-ind, ]
  
  label_test = test[,1]
  training_test = test[,-1]
  label_training = training[,1]
  training_training = training[,-1]
  
  # tforest <- randomForest(x = trainingData[,-1], y = trainingData[,1])
  forest = randomForest(x = training_training, y = label_training)
  test_prediction <- predict(forest, training_test)
  
  # transform into boolean vector
  test_prediction <- test_prediction == T
  label_test <- label_test == T
  
  
  sum(test_prediction == label_test)
  
  bTP[i] = sum(test_prediction[label_test] == T)
  bTN[i] = sum(!test_prediction[!label_test] == T)
}
bFN <- 50 - bTP
bFP <- 50 - bTN
acc <- (bTP + bTN)/(bTP + bTN + bFN + bFP)
ind = order(acc)[((n*0.025)+1):(n*0.975)]
bsensitivity <- (bTP/(bTP + bFN))[ind]
bspecificity <- (bTN/(bTN + bFP))[ind]
baccurracy <- acc[ind]

min(bsensitivity)
min(bspecificity)
min(baccurracy)

max(bsensitivity)
max(bspecificity)
max(baccurracy)

