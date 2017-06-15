library(seqinr)
library(boot)
reads<-read.fasta(file = 'fna-files/fna-files/3.TCA.454Reads.fna',as.string=TRUE)

read_Cat<- sapply(reads,function(x){if(substr(x[1],1,10)=='acgagtgcgt'){return (1)} else if(substr(x[1],1,10)=='aggctcgaca'){return (2)}else {return (3)}})
clone_mix_reads_raw <- reads[read_Cat==1]
clone_mix_reads_seq <- sapply(clone_mix_reads_raw,function(x){substring(x[1],11)})

## Count duplicated reads for U
read_counts<-list()
read_counts[unique(clone_mix_reads_seq)]=1
for( i in clone_mix_reads_seq[duplicated(clone_mix_reads_seq)])
{
  read_counts[i]<-as.numeric(read_counts[i])+1
  
}
U<- sapply(read_counts,function(x){x})


references<- read.fasta('./clonal_sequences.fasta',as.string=T)
referencestrings <-sapply(references,function(x){x[1]})


#######
# 2
#######

occ<-matrix(nrow=length(U),ncol=length(referencestrings),dimnames=(list(read=names(U),reference=names(referencestrings))))

for (i in 1:10)
{
  print(i)
  occ[,i]<-sapply(names(U),function(read){grepl(read,referencestrings[i])})
}



# Subset the occurance matrix and the count vector of the sequences
occ_bool<-as.logical(apply(occ,1,sum))
occ_sub<-occ[occ_bool,]


u_sub <- U[names(U) %in% rownames(occ_sub)]
Ur <-u_sub[order(names(u_sub))]
occ_sub<-occ_sub[order(rownames(occ_sub)),]





############
# 3
############


ph_init <- apply(occ_sub,2,function(col){sum(col)/sum(occ_sub)})



our.EM<-function(read_map_matrix, ph,Ur,epsilon=1e-16)
  {
  Pr.R_r.cond.H_h<-apply(read_map_matrix,2,function(col)
    {
      if(sum(col))
      {
        return(1/sum(col))
      }
      else
      {
        return (0)
      }
    }
    )
  read_prob_matrix <- Pr.R_r.cond.H_h*t(read_map_matrix)
  likelihood <-c()
  
  
  while(T)
    {
    Pr.R_r<-apply(read_map_matrix,1,
                  function(row)
                  {
                    sum(
                      Pr.R_r.cond.H_h*row*ph
                    )})
   
    
    numerator<- (ph*read_prob_matrix)
    
    U_rh<- t(numerator)*(Ur/Pr.R_r)
  
    ph_head<-apply(U_rh,2,sum)/sum(Ur)
    if(sum(is.nan(ph_head)))
    {browser()}

    likelihood<-c(likelihood,sum(Ur*log(Pr.R_r)))
    if(all(abs(ph-ph_head)<epsilon))
    {
      return(list(ph=ph_head,likelihoods=likelihood,epsilon=epsilon))
    }
    else
    {
      ph=ph_head
      }
  
    
    }
  
  }
  

a=our.EM(occ_sub,ph_init,Ur)
png('./barlot.EM.png')
barplot(a$ph,las=2,ylab='Properbility',main='Properbilitys of haplotypes for the dataset (337 reads)')
dev.off()
png('./likelihoods.EM.png')
plot(a$likelihoods,ylab='likelihood',xlab='iteration')
dev.off()



#####
# Task 4
#####

bootstrap_resample_dataset <-function(data,indeces)
{ 
  bootstrap_sub<-data[sort(unique(indeces)),]
  
  count<- table(indeces)
  newUr<- as.vector(count)
  names(newUr)<-names(data)[as.numeric(names(count))]
  
  ph_i<- rep(1/10,10)
  
  r=our.EM(bootstrap_sub,ph_i,newUr,1e-5)

  return(r$ph)
}

bootobj<- boot(occ_sub,bootstrap_resample_dataset,1000,parallel ="multicore",ncpus=4)

res= matrix(ncol=3,nrow=ncol(occ_sub),dimnames=list(haplotypes=colnames(occ_sub),value=c('mean','conf_intervall_low_95','conf_intervall_up_95')))

res<-c()
for (i in 1:10)
{
  conf_int<-boot.ci(bootobj,type='basic',index=i)
  res<- c(res,mean(bootobj$t[,i]),conf_int$basic[4:5])
  
}
r= matrix(res,ncol=3,dimnames=list(haplotypes=colnames(occ_sub),value=c('mean','conf_intervall_low_95','conf_intervall_up_95')))

print(r)

#############
# 5
#############

## How reliable are results based on the bootstrapping?
png('boxplot_boot.png')
a=bootobj$t
colnames(a)=colnames(occ_sub)
boxplot(a,ylab='properbility',col=c('white','red','red','white','red','red','white','white','white','white'),las=2)
dev.off()


