setwd('C:/Users/Ben/Documents/Uni/Master/aml/AML/task_03/')

microarry_data <- read.csv('./data/crc_220_microarray.txt', sep='\t')
labels <- read.csv('./data/crc_patients_all.txt', sep='\t')
clinical_data <- read.csv('./data/crc_clinical_sheet.txt', sep='\t')


# preprocessing

microarray_matrix<- as.matrix(microarry_data[(-1:-2)])
rownames(microarray_matrix)<-microarry_data$CLID

sum_na <- function (row)
{
  return(sum(is.na(row)))
}


num_na<- apply(microarray_matrix,1,sum_na)
hist(num_na)

sum(num_na<4)

microarray_matrix_filtered<-microarray_matrix[(num_na<4),]

variance<-apply(microarray_matrix_filtered,1,var,na.rm=T)
plot(variance)
sum(variance>5)

top_var<-microarray_matrix_filtered[variance>5,]
plot(density(top_var[1,]),xlim=c(-14,10),ylim=c(0,0.4),lwd=2,main='Density of Genes with variance >5',xlab="Expression-Value",col=rainbow(16)[1])
for (i in 2:16)
  {lines(density(top_var[i,],na.rm = T),col=rainbow(16)[i],lwd=2)}

# Task 1.2

# PCA not able to handle NA

microarray_matrix_filtered<-microarray_matrix[(num_na<1),]

pca <- prcomp(t(microarray_matrix_filtered),center=T,scale=T)


# TO dOOOOO!

# Task 1.3

pc1<- pca$rotation[order(abs(pca$rotation[,1]),decreasing = T)[1:5],1]
pc2<- pca$rotation[order(abs(pca$rotation[,1]),decreasing = T)[1:5],2]

# Task 1.4
color_key<-function(patient)
{

  return(ifelse(labels[which(labels[,1]==patient)[1],2],"blue","red"))}

lables2<-list(labels$cancer)
names(labels2)<-labels$patient

paste0(strsplit("TCGA.AG.4015.01A.01R.1119.07",'\\.')[[1]][1:3], collapse = '-')

pca_plot<- function (comp1, comp2)
  {
  plot(pca[comp2],pca[comp2],col=lapply(colnames(microarray_matrix_filtered),color_key))
  
  }