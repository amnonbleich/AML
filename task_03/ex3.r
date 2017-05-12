# Ben Wulf, Lee Hong, Amnon Bleich

#setwd('AML/task_03/')

require(pls)
require(glmnet)
require(ROCR)

# load the give scripts
source('./pca-utils.r')
source('./pls-varimp.r')

# load the datafiles
microarry_data <- read.csv('./data/crc_220_microarray.txt', sep='\t')
labels <- read.csv('./data/crc_patients_all.txt', sep='\t')
clinical_data <- read.csv('./data/crc_clinical_sheet.txt', sep='\t')


# data preprocessing

# we create a matrix with the expression values
microarray_matrix<- as.matrix(microarry_data[(-1:-2)])

# add the gene names to the matrix
rownames(microarray_matrix)<-microarry_data$CLID


# function that gets a vector and returns the number of NA in it
sum_na <- function (row)
{
  return(sum(is.na(row)))
}

# create a vector with the number of NAs for each gene (row)
num_na<- apply(microarray_matrix,1,sum_na)

# how many genes have no NA values
sum(num_na==0) # 11910

# produce a submatrix which contains only genes without NA, we decide to use only genes without any na because the pca don't support na. So it's the easiest.
microarray_matrix_filtered<-microarray_matrix[(num_na==0),]


# The problem is that the patient data in the expression data are not equal to the clinical data. Thus, we have to remap them

# The function converts name in the expression for the patients to the format in the clinical dataset
rename<-function(elem)
{
  return(paste0(strsplit(elem,'\\.')[[1]][1:3], collapse = '-'))
}

# get alternitive names for the patients and create a dataframe with originale name and the new name
newnames<- data.frame(matrixLabels=colnames(microarray_matrix_filtered),annotationname=sapply(colnames(microarray_matrix_filtered),rename))

# merge the new names of patients with the names in the clinical dataset
merged <-merge(newnames,clinical_data,by.x = 'annotationname',by.y = 'patient')

# sort the dataset in lexicographic order for the old names -> makes coloring easier
merged<- merged[order(merged$matrixLabels),]



##############
## Task 1.1 ##
##############

# the genes with the biggest Variance are the most critical ones. if they don't have the variance of 1 and mean of 0 we have to correct them.

# get for each gene the variance
variance<-apply(microarray_matrix_filtered,1,var,na.rm=T)

# select the genes with a variance bigger then 5 (magic experimental number) and create a subset for them
top_var<-microarray_matrix_filtered[variance>5,]

# Plot for these genes the density
png('densitys.png')
plot(density(top_var[1,]),xlim=c(-10,10),ylim=c(0,0.4),lwd=2,main='Density of Genes (Top Variance) without NA and variance >5',xlab="Expression-Value",col=rainbow(7)[1])
for (i in 2:7)
{lines(density(top_var[i,],na.rm = T),col=rainbow(7)[i],lwd=2)}
dev.off()
# The variance is not 1 and the mean not 0 so we have to correct that for the pca

##############
## Task 1.2 ##
##############

# apply a pca for the data
# we have to transpose the matrix, otherwise we try to seperate the genes instead of the patients
# center set the mean to 0
# scale set the variance in the range -1,+1

pca <- prcomp(t(microarray_matrix_filtered),center=T,scale=T)


# produce the screeplot
png('screeplot.png')
screeplot_percent(pca)
dev.off()

# The screeplot shows that the first 10 varables declares only 40% of the variance.
# That means that the pca is not able to seperate the values with a good quality
#

##############
## Task 1.3 ##
##############

# get loadings from pc 1 und pc 2
# order the loadings by max positive and negative influence on the components
pc1<- pca$rotation[order(abs(pca$rotation[,1]),decreasing = T)[1:5],1]

# These are the genes with the biggest influence on PC 1
# A_23_P211600 A_23_P201529  A_23_P78976 A_32_P138556 A_23_P250044 
# 0.02452318   0.02408349   0.02336386   0.02319608  -0.02319463
# the loadings are small because they don't have a huge influence on pc!


pc2<- pca$rotation[order(abs(pca$rotation[,1]),decreasing = T)[1:5],2]
# These are the genes with the biggest influence on PC 2
# A_23_P211600  A_23_P201529   A_23_P78976  A_32_P138556  A_23_P250044 
# 4.219083e-05  5.781348e-04  4.194760e-04 -2.799018e-03 -7.689996e-03 

##############
## Task 1.4 ##
##############


# get different colors for the datapoints if the cancer type is colon we use a red dot otherwise a blue
# the order command is important to get the right color for each patient
colors<- ifelse(merged$cancer[order(rownames(pca$x))]=="Colon",'red','blue')

# perform the multiple pca plot
png('pca_cancer.png')
splom_pca(pca, col=colors)
dev.off()

# What can we see?
# The classes are not realy separeated, but there is a trend between PC1 and the others. The rectal cancer is tendential on the right side.


##############
## Task 1.5 ##
##############

# color the gender
colors2<- ifelse(merged$gender[order(rownames(pca$x))]=="MALE",'lightblue','pink')
png('pca_gender.png')
splom_pca(pca, npcs = 5, col = colors2)
dev.off()
# color age
ramp <- colorRampPalette(c("orange", "blue"))

colors3<- ramp(58)[merged$age_at_initial_pathologic_diagnosis[order(rownames(pca$x))]-34]

png('pca_age.png')
splom_pca(pca, npcs = 5, col = colors3)
dev.off()




##############
## Task 2.2 ##
##############

cancerTypedummy<-ifelse(merged$cancer[order(rownames(pca$x))]=='Colon',1,0)

cor_mat<-cor(t(microarray_matrix_filtered), y = cancerTypedummy)

max(cor_mat) #  0.3999388
min(cor_mat) # -0.5037711

# There are some genes which sceems to be correlated to the cancertype
# Thus, it shut be possible do find a linear correlation between some genes and the cancer type with plsr which asummes a linear correlation


##############
## Task 2.3 ##
############## 

# perform pls regression with 10 fold Crossvalidation
pls_regression<- plsr(cancerTypedummy ~ t(microarray_matrix_filtered), ncomp=10,scale=T, validation='CV',segments=10, segment.type="consecutive",center=T)


##############
## Task 2.4 ##
##############

# performance plots for lasso regression
lasso_regression <- glmnet(t(microarray_matrix_filtered), cancerTypedummy, family="gaussian", alpha=1)

# Perform the Principle component regression
pc_regression <- pcr(cancerTypedummy ~ t(microarray_matrix_filtered), ncomp=6, validation="CV",scale=T,center=T)

#########
# Plots #
#########

# prediction for Lasso
lasso_prob <- predict(lasso_regression, type="response", newx=t(microarray_matrix_filtered))
lasso_pred <- prediction(lasso_prob[,1], cancerTypedummy)

# prediction for PC-Regression
pcr_prob <- predict(pc_regression, type="response", newx=t(microarray_matrix_filtered))
pcr_pred <- prediction(pcr_prob[1:220], cancerTypedummy)

# prediction for PLS
pls_prob <- predict(pls_regression, type="response", newx=t(microarray_matrix_filtered))
pls_pred <- prediction(pls_prob[1:220], cancerTypedummy)


# ROC Curves
lasso_perf1 <- performance(lasso_pred, "tpr", "fpr")
pcr_perf1 <- performance(pcr_pred, "tpr", "fpr")
pls_perf1 <- performance(pls_pred, "tpr", "fpr")


png("ROC.png")
par(mfrow=c(1,3))
plot(lasso_perf1, main="ROC curve for lasso")
plot(pcr_perf1, main="ROC curve for PCR")
plot(pls_perf1, main="ROC curve for pls")
dev.off()



# sensitivity/specificity curve for lasso curves

lasso_perf2 <- performance(lasso_pred, "sens", "spec")
pcr_perf2 <- performance(pcr_pred, "sens", "spec")
pls_perf2 <- performance(pls_pred, "sens", "spec")

png("sens_spec.png")
par(mfrow=c(1,3))
plot(lasso_perf2, main="sensitivity/specificity for lasso")
plot(pcr_perf2, main="sensitivity/specificity for PCR")
plot(pls_perf2, main="sensitivity/specificity for pls")
dev.off()

