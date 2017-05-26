## Task 6
library(bnlearn)
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

# arc operations.
par(mfrow=c(2,2))
plot(constraint_based,main='constraint based')
a<- set.arc(constraint_based, from='A', to='F', check.cycles = TRUE, check.illegal = TRUE, debug = FALSE)
plot(a,main='constraint based add edge from A to F')
b<-reverse.arc(a, from='F', to='A', check.cycles = TRUE, check.illegal = TRUE, debug = FALSE)
plot(b,main='constraint based reverse edge from A to F')
c<-drop.arc(b, from='A', to='F', debug = FALSE)
plot(c,main='constraint based remove edge from F to A')


# turn to dag

dag <- set.arc(constraint_based, from='A', to='B', check.cycles = TRUE,check.illegal = TRUE, debug = FALSE)
