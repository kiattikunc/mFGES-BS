
####### check ayclic####

remove_weak_cycle <- function(checked_pag,learned_mag,posterior_int_weak,n){

  
#  
#  
# 

find.cycles <- function(graph, k) {
  ring <- graph.ring(k, TRUE)
  subgraph_isomorphisms(ring, graph)
}

#find all cycles
g <- graph.adjacency(checked_pag)
l <- unlist(lapply(1L:n, find.cycles, graph=g), recursive=FALSE)

#extract the vertices in each cycle
x <- l[lapply(l, length) >1]



Matrix_x <- matrix(as.numeric(unlist(x)), ncol = 2, byrow = TRUE)
Matrix_x <- Matrix_x[!duplicated(Matrix_x[ ]),]
Matrix_x <- cbind(Matrix_x,matrix(0, nrow(Matrix_x), 2))

colnames(Matrix_x) <- c("from","to","weight","remove")

rows <- unique(which(Matrix_x[, 1] == Matrix_x[, 2], arr.ind = TRUE))
if (length(rows)!=0)
{
  Matrix_x <- Matrix_x[-rows, ]
}

for (index_weak in 1:nrow(Matrix_x))
{
  
  Matrix_x[index_weak,3] =posterior_int_weak[Matrix_x[index_weak,1] ,Matrix_x[index_weak,2]]
}

rows2 <- unique(which(Matrix_x[,3] >=0, arr.ind = TRUE))
if (length(rows2)==nrow(Matrix_x))
{
  learned_mag[Matrix_x[1,1] ,Matrix_x[1,2]]<-0 
}

#then you can do whatever you want
#remove rows

#remove columns

if (!is.null(nrow(Matrix_x)))
{
  Matrix_x <- Matrix_x[order(Matrix_x[,3],decreasing=TRUE),]

} else {
  Matrix_x <- t(Matrix_x)
}  

#colnames(learned_mag)[20]
#colnames(learned_mag)[33]
#colnames(learned_mag)[32]
#learned_mag[33,19]
#learned_mag[29,19]
if (is.null(nrow(Matrix_x)))
{
  
  
  n_matrix <-1
  
  while(learned_mag[Matrix_x[n_matrix,1],Matrix_x[n_matrix,2]]==0)
  {
    
    n_matrix <-n_matrix+1
  }
  
  
  learned_mag[Matrix_x[n_matrix,1] ,Matrix_x[n_matrix,2]]<-0
  
}




#g <- graph.adjacency(checked_pag)
#plot(g)


return_list <- list('remove_graph'=learned_mag,'cycle'=Matrix_x)



}
############ END CHECK CYCIC ################