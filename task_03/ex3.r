setwd('C:/Users/Ben/Documents/Uni/Master/aml/AML/task_03/')

microarry_data <- read.csv('./data/crc_220_microarray.txt', sep='\t')
labels <- read.csv('./data/crc_patients_all.txt', sep='\t')
clinical_data <- read.csv('./data/crc_clinical_sheet.txt', sep='\t')


# preprocessing

microarray_matrix<- as.matrix(microarry_data[(-1:-2)])


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
plot(density(top_var[1,]),xlim(?pc))
for (i in 2:16){lines(density(top_var[i,]))}

