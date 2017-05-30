## Task 6 AML - Ben Wulf, Amnon Bleich, Lie Hong - 2.6.17

# 4.a)
library(bnlearn)


#4.b)
data(learning.test)
constraint_based<-gs(learning.test) # Grow-Shrink
scoring_based <- hc(learning.test)  # Hill Climb

compare(constraint_based,scoring_based)
# scoring based
# arcs:                                5 
# undirected arcs:                     0 
# directed arcs:                       5 
# average branching factor:            0.83 
###################################################
# contraints based
# arcs:                                5 
# undirected arcs:                     1 
# directed arcs:                       4 
# average branching factor:            0.67 

# one undirected edge in the constrinats based leads to a lower branching factor 
##

# look on edge weigths ss ß

par(mfrow=c(1,2))

plot(constraint_based,main="constraint based")
plot(scoring_based,main="scoring based")



# 4.c)

# set, reverse, remove arcs with bnlearn.
par(mfrow=c(2,2))
plot(constraint_based,main='constraint based')
a<- set.arc(constraint_based, from='A', to='F', check.cycles = TRUE, check.illegal = TRUE, debug = FALSE)
plot(a,main='constraint based add edge from A to F')
b<-reverse.arc(a, from='F', to='A', check.cycles = TRUE, check.illegal = TRUE, debug = FALSE)
plot(b,main='constraint based reverse edge from A to F')
c<-drop.arc(b, from='A', to='F', debug = FALSE)
plot(c,main='constraint based remove edge from F to A')


# turn to dag

plot( set.arc(constraint_based, from='A', to='B', check.cycles = TRUE,check.illegal = TRUE, debug = FALSE),main='turned to DAG')



#4.D
# using some prior knowlage to improve graph

par(mfrow=c(1,2))
blacklisted_arcs<- data.frame(from=c('A','B'),to=c('B','A'))
cb2<-gs(learning.test, undirected=F, whitelist=c("E","F"), blacklist=blacklisted_arcs)
plot(cb2,main='Blacklisted arc A-B; whitelisted arc E->F')
cb3<-gs(learning.test, undirected=F, whitelist=c("A","F"), blacklist=c("C","D"))
plot(cb3,main='whitelist A->F; blacklist C->D')



par(mfrow=c(1,1))
# 4.E) Hill-climbing algorithm:

set.seed(3)
#Start with empty network G
nodes <- colnames(learning.test)
hc.graph <- empty.graph(nodes)
hc.score <- score(hc.graph, learning.test)

for (i in 1:1000) {
  manipulation.random <- sample(1:2,1)
  nodes.random <- sample(1:6,2,replace=F) # Get a random arc
  
  if (manipulation.random == 1)       #Try to add the choosen arc
    {
    prop.graph <- tryCatch(set.arc(hc.graph, nodes[nodes.random[1]], nodes[nodes.random[2]]), error=function(ex){hc.graph})
    }
  
  if (manipulation.random == 2) {     # Try to remove choosen arc
    prop.graph <- drop.arc(hc.graph, nodes[nodes.random[1]], nodes[nodes.random[2]])
  }
  
  #Compute the posterior of the neighborhood of G
  prop.score <- score(prop.graph, learning.test)
  
  #Select graph with highest score
  if (prop.score > hc.score) {
    hc.graph <- prop.graph
    hc.score <- prop.score
  }
  
  #Repeat
}

plot(hc.graph,main='own simple HC')
