###############################################################################
ugraphToTetradGraph <- function(ugmat, node_list){
  numNodes <- ncol(ugmat)
  varnames <- node_list
  edgelist <- c()

  for (i in 2:numNodes)
  {
    for (j in 1:(i-1))
    {
      if (ugmat[i,j]==3 & ugmat[j,i]==3) 
      {
        edgelist <- c(edgelist,paste(varnames[j],"->",varnames[i],"[arrowtail=none, arrowhead=none];"))
      }
      
      else if (ugmat[i,j]==1 & ugmat[j,i]==0)
      {edgelist <- c(edgelist,
                     paste(varnames[i],
                           "->",
                           varnames[j],"[arrowtail=none, arrowhead=normal];"))
      }
      
      else if (ugmat[i,j]==0 & ugmat[j,i]==1)
      {edgelist <- c(edgelist,
                     paste(varnames[j],
                           "->",
                           varnames[i],"[arrowtail=none, arrowhead=normal];"))
      }
      
      
      else if (ugmat[i,j]==1 & ugmat[j,i]==1) 
      {edgelist <- c(edgelist,
                     paste(varnames[i],
                           "->",
                           varnames[j],"[dir=both] ;"))
      }
      
      
    }
    
  }
  

    return( paste0("digraph g {",gsub(",","",toString(edgelist)),"}"))

}
