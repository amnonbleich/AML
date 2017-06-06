#source("https://bioconductor.org/biocLite.R")
#biocLite("nem")
library('nem')
data(BoutrosRNAi2002)

D <- BoutrosRNAiDiscrete[,9:16]

control = set.default.parameters(unique(colnames(D))) 
res <- nem(D,inference="search", control=control,verbose = T)

# Use the exhaustive search inference algorithm to find the best scoring model
# task 1.1
print(res)
#$mLL: -239.5778
plot.nem(res)


#To which S genes are the E genes most likely connected to?

best_model <- res$pos[[334]]
best <- apply(best_model,1,function(x){colnames(best_model)[which.max(x)]})
# task 1.2
print(best)


# task 1.3

par(mfrow=c(2,3))
models <- enumerate.models(c("rel","key", "tak", "mkk4hep"))


top5 <- order(res$mLL,decreasing = T)[2:5]
i=1
for (elem in top5)
{
  
  plot(graph.adjacency(models[[elem]]-diag(4)),main=paste('Top',i))
  i=i+1
}
#as seen in the plots 4th and 5th models are similar to the best one, 2nd and 3rd differ from the best 

# task 1.4
steps<- seq(0.01,0.99,0.03)

for (alpha in steps){}


#2.1a)
par(mfrow=c(1,3))
res2 <- nem(D,inference="pairwise", control=control)
plot(res2, main="pairwise")
print(res2)
#$mLL: -256.9889
#Summary of MAP estimates:
#all  ..  -> <-> 
#  6   2   3   1

#2.1b)
res3 <- nem(D,inference="triples", control=control)
plot(res3,main="triples")
print(res3)
# $mLL: -239.5778

#2.1c)
res3 <- nem(D,inference="nem.greedy", control=control)
plot(res3,  main="nem.greedy")
print(res3)
# $mLL: -239.5778

## The three algorithms end up with the same graph.
# The exaustive search, triples , nem.greedy have the same mLL of  -239.5778 the pairwise have a better mLL of -256

##2.2
# Still missing


#4.1

all_models<- enumerate.models(c("rel","key", "tak", "mkk4hep"))
statis<- sapply(all_models,function(v){v[4,3]&v[2,1]})
