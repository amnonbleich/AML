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
