library(seqinr)
reads<-read.fasta(file = 'fna-files/fna-files/3.TCA.454Reads.fna',as.string=TRUE)

read_Cat<- sapply(reads,function(x){if(substr(x[1],1,10)=='acgagtgcgt'){return (1)} else if(substr(x[1],1,10)=='aggctcgaca'){return (2)}else {return (3)}})
clone_mix_reads_raw <- reads[read_Cat==1]
clone_mix_reads_seq <- sapply(clone_mix_reads_raw,function(x){substring(x[1],11)})
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

occ_bool<-as.logical(apply(occ,1,sum))
occ_sub<-occ[occ_bool,]

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
  likelyhood <-c()
  
  
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

    likelyhood<-c(likelyhood,sum(Ur*log(Pr.R_r)))
    if(all(abs(ph-ph_head)<epsilon))
    {
      return(list(ph=ph_head,likelyhoods=likelyhood,epsilon=epsilon))
    }
    else
    {
      ph=ph_head
      }
  
    
    }
  
  }
  





u_sub <- U[names(U) %in% rownames(occ_sub)]
Ur <-u_sub[order(names(u_sub))]
occ_sub<-occ_sub[order(rownames(occ_sub)),]


a=our.EM(occ_sub,ph_init,Ur)


bootres<-list()


for (i in 1:1000)
{ print(i)
  samples<-sample(nrow(occ_sub),nrow(occ_sub),replace=T)
  bootstrap_sub<-occ_sub[sort(unique(samples)),]
  
  count<- table(samples)
  newUr<- as.vector(count)
  names(newUr)<-names(occ_sub)[as.numeric(names(count))]
  
  ph_i<- rep(1/10,10)
  
  r=our.EM(bootstrap_sub,ph_i,newUr,1e-5)

  bootres[i]<-list(r$ph)
}


