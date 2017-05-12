setwd('~/Workspace/Machine_Learning/AML/task_04/')
# load the datafiles
clinical_data <- read.csv('./crc_clinical_sheet.txt', sep='\t')

sum_na <- function (to_sum)
{
  return(sum(is.na(to_sum)))
}

# create a vector with the number of NAs for each feature (row)
#num_na = sum(is.na(row)
num_na_feature<- apply(clinical_data,2,sum_na)
num_na_patient<- apply(clinical_data,1,sum_na)

#Remove features with 30% or more N/A and patients with at least one N/A
clinical_data_filtered<-clinical_data[(num_na_patient==0),(num_na_feature<nrow(clinical_data)*0.3)]
