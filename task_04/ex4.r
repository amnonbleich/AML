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
# clinical_data_filtered$days_to_last_followup<-NULL
clinical_data_filtered$lymphnodes_examined<-NULL
# clinical_data_filtered$number_of_lymphnodes_examined<-NULL


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
success_rate=0
for (i in (1:100))
{
  smp_size <- floor(0.75 * nrow(clinical_data_filtered))
  training_idxs = sample(nrow(clinical_data_filtered), size = smp_size)
  training_data = clinical_data_filtered[training_idxs,]
  test_data = clinical_data_filtered[-training_idxs, ]
  svm_model <- svm(vascular_invasion_present ~ ., data = training_data)
  svm_pred <- predict(svm_model, test_data)
  correct = test_data$vascular_invasion_present
  success_rate=success_rate+sum(correct == svm_pred)/nrow(test_data)
}
success_rate=success_rate/100



# Useless

# clinical_data_filtered$cancer <- as.numeric(clinical_data_filtered$cancer=="Colon")
# 
# tmp_model <- model.matrix(~as.numeric(anatomic_organ_subdivision), clinical_data_filtered)
# tmp_model <- tmp_model[match(rownames(clinical_data_filtered), rownames(tmp_model)),] #fill in corresponding NA's
# clinical_data_filtered$anatomic_organ_subdivision <- tmp_model[,2]
# 
# tmp_model <- model.matrix(~as.numeric(tumor_site), clinical_data_filtered)
# tmp_model <- tmp_model[match(rownames(clinical_data_filtered), rownames(tmp_model)),] #fill in corresponding NA's
# clinical_data_filtered$tumor_site <- tmp_model[,2]
# 
# tmp_model <- model.matrix(~as.numeric(distant_metastasis_pathologic_spread), clinical_data_filtered)
# tmp_model <- tmp_model[match(rownames(clinical_data_filtered), rownames(tmp_model)),] #fill in corresponding NA's
# clinical_data_filtered$distant_metastasis_pathologic_spread <- tmp_model[,2]
# 
# clinical_data_filtered$gender <- as.numeric(clinical_data_filtered$gender=="MALE")
# 
# #clinical_data_filtered$histological_type <- as.numeric(clinical_data_filtered$histological_type)
# 
# tmp_model <- model.matrix(~as.numeric(histological_type), clinical_data_filtered)
# tmp_model <- tmp_model[match(rownames(clinical_data_filtered), rownames(tmp_model)),] #fill in corresponding NA's
# clinical_data_filtered$histological_type <- tmp_model[,2]
# 
# tmp_model <- model.matrix(~as.numeric(history_of_colon_polyps), clinical_data_filtered)
# tmp_model <- tmp_model[match(rownames(clinical_data_filtered), rownames(tmp_model)),] #fill in corresponding NA's
# clinical_data_filtered$history_of_colon_polyps <- tmp_model[,2]
# 
# clinical_data_filtered$icd_o_3_histology <- as.numeric(clinical_data_filtered$icd_o_3_histology=="8140/3")
# 
# tmp_model <- model.matrix(~as.numeric(icd_o_3_site), clinical_data_filtered)
# tmp_model <- tmp_model[match(rownames(clinical_data_filtered), rownames(tmp_model)),] #fill in corresponding NA's
# clinical_data_filtered$icd_o_3_site <- tmp_model[,2]
# 
# clinical_data_filtered$lymphatic_invasion_present <- as.numeric(clinical_data_filtered$lymphatic_invasion_present=="YES")
# 
# tmp_model <- model.matrix(~as.numeric(lymphnode_pathologic_spread), clinical_data_filtered)
# tmp_model <- tmp_model[match(rownames(clinical_data_filtered), rownames(tmp_model)),] #fill in corresponding NA's
# clinical_data_filtered$lymphnode_pathologic_spread <- tmp_model[,2]
# 
# clinical_data_filtered$lymphnodes_examined <- as.numeric(clinical_data_filtered$lymphnodes_examined=="YES")
# clinical_data_filtered$person_neoplasm_cancer_status <- as.numeric(clinical_data_filtered$person_neoplasm_cancer_status=="WITH TUMOR")
# 
# tmp_model <- model.matrix(~as.numeric(primary_tumor_pathologic_spread), clinical_data_filtered)
# tmp_model <- tmp_model[match(rownames(clinical_data_filtered), rownames(tmp_model)),] #fill in corresponding NA's
# clinical_data_filtered$primary_tumor_pathologic_spread <- tmp_model[,2]
# 
# clinical_data_filtered$prior_diagnosis <- as.numeric(clinical_data_filtered$prior_diagnosis=="YES")
# 
# tmp_model <- model.matrix(~as.numeric(residual_tumor), clinical_data_filtered)
# tmp_model <- tmp_model[match(rownames(clinical_data_filtered), rownames(tmp_model)),] #fill in corresponding NA's
# clinical_data_filtered$residual_tumor <- tmp_model[,2]
# 
# clinical_data_filtered$synchronous_colon_cancer_present <- as.numeric(clinical_data_filtered$synchronous_colon_cancer_present=="YES")
# 
# tmp_model <- model.matrix(~as.numeric(tumor_stage), clinical_data_filtered)
# tmp_model <- tmp_model[match(rownames(clinical_data_filtered), rownames(tmp_model)),] #fill in corresponding NA's
# clinical_data_filtered$tumor_stage <- tmp_model[,2]
# 
# 
# clinical_data_filtered$vascular_invasion_present <- as.numeric(clinical_data_filtered$vascular_invasion_present=="YES")

