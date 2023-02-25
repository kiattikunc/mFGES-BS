posteriorcal <- function(prior,indep_graphs,Data,n,INT){
  
  prior_infor<- prior
  diag(prior_infor)<-0
  true_mag_amat <- amat(true_mag)
  rownames(prior_infor)<- rownames(prior)
  colnames(prior_infor)<- colnames(prior)
  depence_posterior<- matrix(0, n, n)
  indepence_posterior<- matrix(0, n, n)
  
  
library(bnlearn)
#use obs to cal the score
#prior_uniform<- matrix(1/2, n, n)



# cprior_uniformatch Bdeu score for constraints
for (run in 1:INT)
{

  for (nodeA in 1:n)
  {
   
    for( nodeB in 1:n)
    {
      if (nodeA != nodeB)
      {
        res = set.arc(indep_graphs,rownames(prior_infor)[nodeA],rownames(prior_infor)[nodeB])
        depend_score <- bnlearn::score(res,Data[[run]]$data, type = "bde")
        res3 <- drop.arc(res, rownames(prior_infor)[nodeA],rownames(prior_infor)[nodeB])
        independ_score <- bnlearn::score(res3, Data[[run]]$data, type = "bde")
        indepence_posterior[nodeA,nodeB] <- ((1-prior_infor[nodeA,nodeB])*depend_score)/(((prior_infor[nodeA,nodeB])*independ_score)+(1-prior_infor[nodeA,nodeB])*depend_score)
        
        
      }
      
    }
    
  }

 
    prior_infor <- pmax((1-indepence_posterior), Data[[run]]$prior)
    #prior_infor [ Data[[run]]$prior==0 & t(Data[[run]]$prior>0.5)] <- 0
    prior_infor [ Data[[run]]$prior==0 &  t(Data[[run]]$prior)  >0.5] <- 0
    prior_infor [ Data[[run]]$prior>0 &  t(int_rule)  ==1] <- 0

    
    rownames(prior_infor)<- rownames(prior)
    colnames(prior_infor)<- colnames(prior)
 

  
}
indepence_posterior<-(1-prior_infor)
diag(indepence_posterior)<-0
rownames(indepence_posterior)<- rownames(prior)
colnames(indepence_posterior)<- colnames(prior)

my_list <- list("indepence_posterior" = indepence_posterior, "indepence_posterior_INT" = indepence_posterior)
return(my_list)
}

multiplymatrix <- function(A,B) {
  size <- nrow(A)
  C<- matrix(0, size, size)
  for (i in 1:size) {
    for (j in 1:size) {
      C[i, j] <- A[i, j]*B[i, j]
     
    }

   
  }
  return(C)
}