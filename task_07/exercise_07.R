#AML Exercise 7 - 08.06.17 - Ben Wulf, Lie Hong, Amnon Bleich

#source("https://bioconductor.org/biocLite.R")
#biocLite("nem")
library('nem')
library(Rgraphviz)
library(gplots)

data(BoutrosRNAi2002)

# remove the control groups
D <- BoutrosRNAiDiscrete[,9:16]

# perform discretizazion to get alpha and beta errors with pseudocounts
res.disc = nem.discretize(BoutrosRNAiExpression,neg.control=1:4,pos.control=5:8)

#######################
# Task 1.1
# Use the exhaustive search inference algorithm to find the best scoring model
#######################
control = set.default.parameters(unique(colnames(D)),para=res.disc$para) 
res <- nem(D,inference="search", control=control,verbose = T)

print(res)
# best model =334
#$mLL: -239.5778
#67  selected E-genes:
#--> rel : 30  attached E-genes
#--> key : 30  attached E-genes
#--> tak : 9  attached E-genes
#--> mkk4hep : 28  attached E-genes
#--> null : 1  attached E-genes
plot.nem(res)

#Task 1.3
#To which S genes are the E genes most likely connected to?

best_model <- res$pos[[334]]
best <- apply(best_model,1,function(x){colnames(best_model)[which.max(x)]})
# task 1.2
print(best)


# task 1.3
# How do the top 5 exhaustive models differ?
Sgenes=unique(colnames(D))
models <- enumerate.models(Sgenes)
best5 <- order(res$mLL,decreasing = T)[1:5]
col<-c("yellow","yellow","green","blue")
names(col) = Sgenes
png(file='Top5Models.png',width = 600, height = 500)
par(mfrow=c(2,3),oma = c( 0, 0, 2, 0 ))

for (i in 1:5) {
  graph <- as(models[[best5[i]]]-diag(4),"graphNEL")
  print(best5[i])
  plot(graph,
       nodeAttrs=list(fillcolor=col),
       main=paste("-",i, "-"))        
}
title('Top 5 Models', outer = T)
dev.off()


# task 1.4
# Evaluate stability for different choices of alpha, beta.
m=0
n=0
mLLM<-matrix(nrow=10, ncol=10) 
diffM<-matrix(nrow=10, ncol=10)

selmodel<- models[best5[1]][[1]]

for (alpha in seq(0.01,0.99,0.1))
{
  m<- m+1
  n=0
  for (beta in seq(0.01,0.99,0.1))
  { n<-n+1
    print(paste(m,n))
    control = set.default.parameters(unique(colnames(D)),para=c(alpha,beta)) 
    resT <- nem(D,inference="search", control=control)
    maxi<-which.max(resT$mLL)
    mLLM[m,n]<-resT$mLL[maxi]
    diffM[m,n]<-sum(models[maxi][[1]]!=selmodel)
    
  }
}


heatmap.2(diffM,dendrogram='none',col=rev(redgreen(12)),Rowv=FALSE, Colv=FALSE,trace='none',labCol=seq(0.01,0.99,0.1),labRow=seq(0.01,0.99,0.1),main='# ofChanges in Models',xlab='beta',ylab='alpha')
heatmap.2(mLLM,dendrogram='none',col=rev(redgreen(256)),Rowv=FALSE, Colv=FALSE,trace='none',labCol=seq(0.01,0.99,0.1),labRow=seq(0.01,0.99,0.1),main='mLL of the best Models',xlab='beta',ylab='alpha')


#2.1a)
par(mfrow=c(1,3))
control = set.default.parameters(unique(colnames(D)),para=res.disc$para) 
res2 <- nem(D,inference="pairwise", control=control)
plot.nem(res2, main="pairwise")
print(res2)
#$mLL: -257.4055
#Summary of MAP estimates:
#all  ..  -> <-> 
#  6   2   3   1

#2.1b)
res3 <- nem(D,inference="triples", control=control)
plot.nem(res3,main="triples")
print(res3)
# $mLL: -239.5778

#2.1c)
res3 <- nem(D,inference="nem.greedy", control=control)
plot.nem(res3,  main="nem.greedy")
print(res3)
# $mLL: -239.5778

## The three algorithms end up with the same graph.
# The exaustive search, triples , nem.greedy have the same mLL of  -239.5778 the pairwise have a worser mLL of -256

##2.2
# Our own divede an conquer approach

edgecount<-matrix(rep(0,16),nrow=4,ncol=4,dimnames=list(col=unique(colnames(D)[1:8]),row=unique(colnames(D)[1:8])))

for(x in seq(1,7,2))
{
  base<- 1:8
  sel<-base[-c(x,x+1)]
  control = set.default.parameters(unique(colnames(D)[sel]),para=res.disc$para) 
  resDC <- nem(D[,sel],inference="search", control=control)
  mod<- as(resDC$graph,'matrix')
  for (row in rownames(mod))
  {
    for (col in colnames(mod))
    {
      edgecount[row,col]<-edgecount[row,col]+mod[row,col]
    }
  }
}

png(file='OwnDC.png',width = 700, height = 600)
par(mfrow=c(1,3),oma = c( 0, 0, 2, 0 ))
col<-c("yellow","yellow","green","blue")
names(col)<-colnames(edgecount)
for (i in 1:3) {
  graph <- as((edgecount>=i),"graphNEL")
  plot(graph,
       nodeAttrs=list(fillcolor=col),
       main=paste("value=>",i))        
}
title('divide and conquer algorithm for triple relations', outer = T)
dev.off()


#######################
#3. Structure Priors

prior1 <- matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1), nrow=4)
colnames(prior1) <- c( "rel", "key", "tak", "mkk4hep")
rownames(prior1) <- c( "rel", "key", "tak", "mkk4hep")
prior2 <- matrix(c(1,0,0,0,1,1,0,0,1,0,1,0,1,0,0,1), nrow=4)
colnames(prior2) <- c( "rel", "key", "tak", "mkk4hep")
rownames(prior2) <- c( "rel", "key", "tak", "mkk4hep")
prior3 <- matrix(c(1,1,1,1,0,1,1,1,0,0,1,1,0,0,0,1), nrow=4)
colnames(prior3) <- c( "rel", "key", "tak", "mkk4hep")
rownames(prior3) <- c( "rel", "key", "tak", "mkk4hep")

p1_mLL <- c()
p2_mLL <- c()
p3_mLL <- c()

for (l in 10^(0:20))
{ print(l)
  control.1 = set.default.parameters(unique(colnames(D)), lambda=l, Pm=prior1)
  control.2 = set.default.parameters(unique(colnames(D)), lambda=l, Pm=prior2)
  control.3 = set.default.parameters(unique(colnames(D)), lambda=l, Pm=prior3)
  resP1 <- nem(D, inference="search", control=control.1)
  resP2 <- nem(D, inference="search", control=control.2)
  resP3 <- nem(D, inference="search", control=control.3)
  
  p1_mLL <- c(p1_mLL, max(resP1$mLL))
  p2_mLL <- c(p2_mLL, max(resP2$mLL))
  p3_mLL <- c(p3_mLL, max(resP3$mLL))
}


plot(10^(0:20),p1_mLL,log='x',xlab='lambda',ylab='mLL',type='l',col='red', main="influence of Lamda")
lines((10^(0:20)),p2_mLL,log='x',col='blue')
lines((10^(0:20)),p3_mLL,log='x',col='green')
abline(v=100)
legend('bottomright',legend=c('Prior 1','Prior 2','Prior 3'),col=c('red','blue','green'),lty=1)



#4.1

###################
# We assume that an directed edge a->b will exclude an edge b->a, otherwise it would be also a directed edge
###################


all_models<- enumerate.models(c("rel","key", "tak", "mkk4hep"))
restricted<- sapply(all_models,function(v){v["mkk4hep","tak"]&v["key","rel"]&!v["tak","mkk4hep"]&!v['rel','key']})
sum(restricted) ## There are 29 Models. If reverse edges allowed there are 48

restricted_models<-all_models[restricted]
restricted_models_order<-order(res$mLL[restricted],decreasing = T)



############
# Plot Top 5 restricted models
############

col<-c("yellow","yellow","green","blue")
names(col) = c("rel","key", "tak", "mkk4hep")
png(file='Top5Models_restricted.png',width = 600, height = 500)
par(mfrow=c(2,3),oma = c( 0, 0, 2, 0 ))

for (i in 1:5) {
  print(i)
  graph <- as(restricted_models[[restricted_models_order[i]]]-diag(4),"graphNEL")

  plot(graph,
       nodeAttrs=list(fillcolor=col),
       main=paste("-",i, "-"))   
  print(i)
}
plot(as(all_models[[which.max(res$mLL)]]-diag(4),"graphNEL"),nodeAttrs=list(fillcolor=col), main="unrestricted")
title('Top 5 Models', outer = T)
dev.off()


###########
# 5
###########
control = set.default.parameters(unique(colnames(D)),para=res.disc$para) 

repeats<-1000
thresh<- 0.95
bt_default<-set.default.parameters(unique(colnames(D)),para=res.disc$para) 
ec<-matrix(rep(0,16),nrow=4,ncol=4,dimnames=list(source=c("rel","key", "tak", "mkk4hep"),target=c("rel","key", "tak", "mkk4hep"))) # stores the edgecount



for (i in 1:repeats)
{
  print(paste('Durchlauf',i))
  rows<- sample(nrow(D),round(nrow(D)*0.5),replace=T) # get genes
  bt_res <- nem(D[rows,],inference="nem.greedy", control=control,verbose = F)
  ec<-ec+as(bt_res$graph,'matrix')
}

edge_properbility<- ec/repeats # edge proberbilitys
likely_edges <-edge_properbility>=thresh # all edges with high proberbility

col<-c("yellow","yellow","green","blue")
names(col) = c("rel","key", "tak", "mkk4hep")
par(mfrow=c(1,2))
plot(as(likely_edges,"graphNEL"),nodeAttrs=list(fillcolor=col), main="own boostrap model")



buildin_bt_model<-nem.bootstrap(D, thresh=thresh, nboot=repeats,inference="nem.greedy",models=NULL,control=bt_default, verbose=TRUE)
plot(buildin_bt_model$graph,nodeAttrs=list(fillcolor=col), main="buildin boostrap model")
##

