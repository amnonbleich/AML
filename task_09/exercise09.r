library(nnet)
library(caret)
raw_mirna_data <- read.csv('./crc_mirna_datamatrix_rpm.tsv',sep='\t')
rownames(raw_mirna_data)<-raw_mirna_data$Gene
raw_mirna_data<-as.matrix(raw_mirna_data[,-1])


normalize<- function(x){
  {
    if(min(x)==max(x))
    {return(x*0)}
    return((x-min(x))/(max(x)-min(x)))}
}

rowwise_scaled<- apply(raw_mirna_data,1,normalize)






variances <- apply(rowwise_scaled,2,var,na.rm=T)
high_variances<-order(variances,decreasing=T)[1:70] # choose the top 10% of mirna's with the highest variance

high_Var_data <- rowwise_scaled[,high_variances]

# calculate correlation between mirna's
correlations<- cor(high_Var_data,use='pairwise')-diag(length(high_variances))
correlations[upper.tri(correlations)]<-0

low_cor_mirna<-apply(correlations,1,function(x){max(abs(x))<=0.5}) # choose only these mirnas with a correlation lower than 0.5
highVar_lowCor_data <- high_Var_data[,low_cor_mirna]

########
## Bind the tumortype to the dataset
########

clinical_sheet <- read.csv('./crc_clinical_sheet.txt',sep='\t')
clinical_sheet <- clinical_sheet[order(clinical_sheet$patient),]

analyzingdata <- as.data.frame(highVar_lowCor_data[order(row.names(highVar_lowCor_data)),])
rownames(analyzingdata)<-sapply(rownames(analyzingdata),function(x){gsub('\\.','-',substr(x,1,12))})
analyzingdata$tumor_site <- clinical_sheet$tumor_site[clinical_sheet$patient %in% rownames(analyzingdata)]

# Remove patients without classification
analyzingdata<-analyzingdata[!is.na(analyzingdata$tumor_site),]


#########
# Task 2
#########

print(ncol(analyzingdata)-1)
# We have 53 features. That means we will need 53 input neurons, for each feature (mirna) one

print(nrow(analyzingdata))
# We keep 250 patients


print(unique(analyzingdata$tumor_site))
# We have four Categories: 3 - left colon, 2 - transverse colon, 1 - right colon, 4 - rectum   this will mean that we need four output neurons.
#########
# Task 3
#########


folds <- createFolds(1:nrow(analyzingdata),k=5)

default_param_accuracy<-c()

for (i in 1:5)
{
  fold<-folds[[i]]
  names(fold)<-rownames(analyzingdata)[fold]
  
  res<-nnet(tumor_site ~ .,data=analyzingdata[-fold,],size=53,MaxNWts=84581)
  pred<-predict(res,analyzingdata[fold,-54])
  mostlikely<- colnames(pred)[max.col(pred)]
  tab<-table(truth=analyzingdata[fold,54],prediction=mostlikely)
  print(tab)
  print(paste('accuracy',round(sum(diag(tab))/sum(tab),3)))
  default_param_accuracy<-c(default_param_accuracy,sum(diag(tab))/sum(tab))
  
}

