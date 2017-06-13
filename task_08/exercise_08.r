library(seqinr)
reads<-read.fasta(file = 'fna-files/fna-files/3.TCA.454Reads.fna',as.string=TRUE)

read_Cat<- sapply(reads,function(x){if(substr(x[1],1,10)=='acgagtgcgt'){return (1)} else if(substr(x[1],1,10)=='aggctcgaca'){return (2)}else {return (3)}})
clone_mix_reads_raw <- reads[read_Cat==1]
clone_mix_reads_seq <- sapply(clone_mix_reads_raw,function(x){substring(x[1],11)})

references<- read.fasta('./clonal_sequences.fasta',as.string=T)
referencestrings <-sapply(references,function(x){x[1]})


#######
# 2
#######

occ<-matrix(nrow=length(clone_mix_reads_seq),ncol=length(referencestrings),dimnames=(list(read=names(clone_mix_reads_seq),reference=names(referencestrings))))

for (i in 1:10)
{
  print(i)
  occ[,i]<-sapply(clone_mix_reads_seq,function(read){grepl(read,referencestrings[i])})
}

occ_bool<-as.logical(apply(occ,1,sum))
occ_sub<-occ[occ_bool,]

############
# 3
############



Pr.R_r<-function(read_matrix,ph)
{
  Pr.R_r.cond.H_h<-apply(read_matrix,2,function(col){1/sum(col)})
  
  return(
    apply(read_matrix,1,
          function(row)
          {
            sum(
              Pr.R_r.cond.H_h*row*ph
            )}
    )
  )
}


own.EM <- function(occurances)
  {
  ph=rep(1/ncol(occurances),ncol(occurances))
  
  Ur=apply(occ_sub,1,sum)
  Pr.R_r.cond.H_h<-apply(occ_sub,2,function(col){1/sum(col)})
  for (i in 1:100)
  {
  # Expectation Step
  Urh<-t(
                        Ur*((ph_i*Pr.R_r.cond.H_h2*t(occ_sub))/Pr.R_r(occ_sub,ph))
                        )
  #Maximisation Step
  ph<-apply(Urh,2,sum)/sum(Ur)
  }
  browser()
  return(ph*Pr.R_r.cond.H_h)
}