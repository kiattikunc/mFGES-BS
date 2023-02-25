checkDependencies <- function( amat_obs, suffStat,  alpha)
{

  verbose =FALSE
  amat_obs_depend <- amat_obs
  amat_pvalue<- matrix(0, n, n)
  rownames(amat_pvalue)<- colnames(dataset)
  colnames(amat_pvalue)<- colnames(dataset)
  
  # 
  # which( colnames(dataset_pcalg)=="DISCONNECT" )
  # which( colnames(dataset_pcalg)=="VENTALV" )
  
  # colnames(dataset_pcalg)[nbrsA]
   # colnames(dataset_pcalg)[23]
  # KINKEDTUBE - VENTTUBE
  # which( colnames(dataset_pcalg)=="KINKEDTUBE" )
  # which( colnames(dataset_pcalg)=="CATECHOL" )
  # 
  # LVEDVOLUME- STROKEVOLUME
  # disCItest(which( colnames(dataset_pcalg)=="LVEDVOLUME" ),which( colnames(dataset_pcalg)=="STROKEVOLUME" ),which( colnames(dataset_pcalg)=="CO" ),suffStat)
  # 
  # disCItest(which( colnames(dataset_pcalg)=="INTUBATION" ),which( colnames(dataset_pcalg)=="SAO2" ), which( colnames(dataset_pcalg)=="CATECHOL" ),suffStat)
  # 
  # 
  # disCItest(which( colnames(dataset_pcalg)=="KINKEDTUBE" ),which( colnames(dataset_pcalg)=="VENTALV" ),which( colnames(dataset_pcalg)=="VENTLUNG" ),suffStat)
  # disCItest(which( colnames(dataset_pcalg)=="KINKEDTUBE" ),which( colnames(dataset_pcalg)=="VENTTUBE" ),NULL,suffStat)
  # disCItest(which( colnames(dataset_pcalg)=="INTUBATION" ),which( colnames(dataset_pcalg)=="VENTTUBE" ),NULL,suffStat)
  # disCItest(which( colnames(dataset_pcalg)=="KINKEDTUBE" ),which( colnames(dataset_pcalg)=="VENTALV" ),NULL,suffStat)
  # disCItest(35,34,NULL,suffStat)
  # 
  # disCItest(34,7,c(35,22),suffStat)
  # disCItest(34,7,NULL,suffStat)
  # disCItest(1,2,NULL,suffStat)
  # nbrsA <- which(amat[,34] != 0)
  # nbrsA <-  which(amat[34,] != 0)
  # S <- nbrsA[which(amat_skeleton[34,] != 0)]
  # 
  for (i in 1:n)
  {
         A <-i
         listB <- which(amat_obs_depend[i,] != 0)
         
         if (length(listB)>0)
         {
           for (j in 1:length(listB))
           {
             B <- listB[j]
             nbrsA <- which(amat_obs_depend[,i] != 0)
             nbrsB <- which(amat_obs_depend[,B] != 0)
             AB <- c(A, B)
             C <- c(nbrsA, nbrsB)
             C <- C[!duplicated(C)]
             C <- C[ !C %in% AB]
             
             
             # ## start with the neighbours of a
             if ((nn <- length(C)) > 0) {
               allComb <- expand.grid(lapply(integer(nn), function(.) 0:1))
               ## loop through all subsets of neighbours
               for (k in 1:nrow(allComb)) { ## == 1:(2^nn)
                 S <- C[which(allComb[k,] != 0)]
                 pval <- disCItest(A, B, S, suffStat)
                 ## save the pval and the set that produced this pval
                 if (verbose) cat("a:",A,"b:",B," S =",S," - pval =",pval,"\n")
                 if (pval >= alpha ) {
                   
                   amat_pvalue[A,B] <-pval
                   amat_pvalue[B,A] <-pval
                 }
               }
               
               
             }
             
           } 
         }
           
              
              
              
              
            
         
        }
   
  
  
  my_list <- list("amat_pvalue" = amat_pvalue)
  return(my_list)

 
} 
