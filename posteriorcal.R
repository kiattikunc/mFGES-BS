posteriorcal <- function(prior,graphs,Data,n){
#use obs to cal the score
#prior_uniform<- matrix(1/2, n, n)
prior_uniform<- prior
diag(prior_uniform)<-0
true_mag_amat <- amat(true_mag)
rownames(prior_uniform)<- rownames(prior)
colnames(prior_uniform)<- colnames(prior)
depence_posterior<- matrix(0, n, n)
indepence_posterior<- matrix(0, n, n)
#amat_round <- round(amat)

# cprior_uniformatch Bdeu score for constraints
for (run in length(Data))
{
  indepence_posterior_tmp<- matrix(0, n, n)
  for (nodeA in 1:n)
  {
   
    for( nodeB in 1:n)
    {
      if (nodeA != nodeB)
      {
        res = set.arc(indep_graphs,rownames(prior_uniform)[nodeA],rownames(prior_uniform)[nodeB])
        depend_score <- score(res,Data[[run]]$data, type = "bde")
        res3 <- drop.arc(res, rownames(prior_uniform)[nodeA],rownames(prior_uniform)[nodeB])
        independ_score <- score(res3, Data[[run]]$data, type = "bde")
    
      
        indepence_posterior_tmp[nodeA,nodeB] <- ((1-prior_uniform[nodeA,nodeB])*independ_score)/(((prior_uniform[nodeA,nodeB])*depend_score)+(1-prior_uniform[nodeA,nodeB])*independ_score)
        
       
        if (run==1)
        {
         indepence_posterior<- indepence_posterior_tmp
        }
        
      }
      
    }
    
  }
 
  indepence_posterior <-multiplymatrix(indepence_posterior,indepence_posterior_tmp)
  
  my_list <- list("indepence_posterior" = indepence_posterior, "indepence_posterior_obs" = indepence_posterior_tmp,'prior_uniform'=prior_uniform)
  return(my_list)
  
}

}

multiplymatrix <- function(A,B) {
  size <- nrow(A)
  C<- matrix(0, size, size)
  for (i in 1:size) {
    for (j in 1:size) {
      C[i, j] <- A[i, j]*B[i, j]
      C[j, i] <- C[i, j]
    }

   
  }
  return(C)
}