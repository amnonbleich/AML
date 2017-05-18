setwd('~/Workspace/Machine_Learning/AML/task_04/')

library(ROCR)
library(e1071)

# load the datafiles
clinical_data <- read.csv('./crc_clinical_sheet.txt', header=TRUE, sep='\t')

sum_na <- function (to_sum)
{
  return(sum(is.na(to_sum)))
}

clinical_data_filtered<-clinical_data
#Remove unneccessary columns:
clinical_data_filtered$patient<-NULL
clinical_data_filtered$lymphnodes_examined<-NULL
clinical_data_filtered$vital_status<-NULL


#Remove features with 30% or more N/A and patients with at least one N/A
num_na_feature<- apply(clinical_data_filtered,2,sum_na)
clinical_data_filtered<-clinical_data_filtered[,(num_na_feature<nrow(clinical_data)*0.3)]
num_na_patient<- apply(clinical_data_filtered,1,sum_na)
clinical_data_filtered<-clinical_data_filtered[num_na_patient==0,]

# varify all categorial are "factor"
d<-(lapply(clinical_data_filtered,class) == "factor" | lapply(clinical_data_filtered,class) == "numeric" | lapply(clinical_data_filtered,class) == "integer")
all(d)
# TRUE


# TODO tune cost and gamma
# success_rate=0
# for (i in (1:500))
# {
#   smp_size <- floor(0.75 * nrow(clinical_data_filtered))
#   training_idxs = sample(nrow(clinical_data_filtered), size = smp_size)
#   training_data = clinical_data_filtered[training_idxs,]
#   test_data = clinical_data_filtered[-training_idxs, ]
#   svm_model <- svm(vascular_invasion_present ~ ., data = training_data)
#   svm_pred <- predict(svm_model, test_data)
#   correct = test_data$vascular_invasion_present
#   success_rate=success_rate+sum(correct == svm_pred)/nrow(test_data)
# }
# success_rate=success_rate/i

#svm_model <- svm(vascular_invasion_present ~ ., data = training_data)
smp_size <- floor(0.75 * nrow(clinical_data_filtered))
training_idxs = sample(nrow(clinical_data_filtered), size = smp_size)
training_data = clinical_data_filtered[training_idxs,]
test_data = clinical_data_filtered[-training_idxs,]

# tune using 10-fild cv, find bets gamma/cost
set.seed(1)
tuned = tune.svm(vascular_invasion_present~., data = training_data, gamma = seq(0.005, 0.05, 0.005), cost = 10^(-3:3), tunecontrol=tune.control(cross=10),    kernel = "linear")
# best = best.svm(vascular_invasion_present~., data = training_data, gamma = seq(0.005, 0.05, 0.005), cost = 10^(0:3), tunecontrol=tune.control(cross=10),    kernel = "linear")

summary(tuned)
plot(tuned)

# - sampling method: 10-fold cross validation 
# 
# - best parameters:
#   gamma cost
# 0.005   1
# 
# - best performance: 0.1686813 error rate

# Parameter tuning of ‘svm’- radial kernel:
#   
#   - sampling method: 10-fold cross validation 
# 
# - best parameters:
#   gamma cost
# 0.005  100
# 
# - best performance: 0.182967 error rate
# 

svm_model <- svm(vascular_invasion_present ~ ., data = training_data, cost=tuned$best.parameters$cost ,gamma=tuned$best.parameters$gamma, kernel="linear")
# remove vascular invation columns for prediction
svm_pred <- predict(svm_model, test_data[,-ncol(test_data)])
correct = test_data$vascular_invasion_present
success_rate=sum(correct == svm_pred)/nrow(test_data)
