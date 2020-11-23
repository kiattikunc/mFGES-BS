#alarm n data sets with input targets
# version 5 to add conservative rule to add v-structure to add bi-direction
sink(file='myoutput.txt')
library(bnlearn)
library(rcausal)
library(DOT)
library(tictoc)
library(pcalg)
library(pheatmap)
library(tictoc)
library(Hmisc)

#hyper parameters
iterFGES <- list()
thred <- 0.5
target <-5
set.seed(5)
n_iter <-2
#load dataset
set <-'trainingData_ALARM_L5_10k.csv'
dataset <- read.csv(set,header = TRUE,na.strings=c(""))
dataset$DISCONNECT <- as.factor(dataset$DISCONNECT)
dataset$ERRCAUTER <- as.factor(dataset$ERRCAUTER)
dataset$ERRLOWOUTPUT <- as.factor(dataset$ERRLOWOUTPUT)
dataset$KINKEDTUBE <- as.factor(dataset$KINKEDTUBE)
dataset$HYPOVOLEMIA <- as.factor(dataset$HYPOVOLEMIA)
dataset$HISTORY <- as.factor(dataset$HISTORY)
dataset$ANAPHYLAXIS <- as.factor(dataset$ANAPHYLAXIS)
dataset$HISTORY <- as.factor(dataset$HISTORY)
dataset$INSUFFANESTH <- as.factor(dataset$INSUFFANESTH)
dataset$PULMEMBOLUS <- as.factor(dataset$PULMEMBOLUS)
dataset <- dataset[,order(colnames(dataset))]
dataset2 <- dataset[,order(colnames(dataset))]

#load interventional dataset
set <-'trainingData_ALARM_N_10k.csv'
truedag_dataset <- read.csv(set,header = TRUE,na.strings=c(""))
truedag_dataset$DISCONNECT <- as.factor(truedag_dataset$DISCONNECT)
truedag_dataset$ERRCAUTER <- as.factor(truedag_dataset$ERRCAUTER)
truedag_dataset$ERRLOWOUTPUT <- as.factor(truedag_dataset$ERRLOWOUTPUT)
truedag_dataset$KINKEDTUBE <- as.factor(truedag_dataset$KINKEDTUBE)
truedag_dataset$LVFAILURE <- as.factor(truedag_dataset$LVFAILURE)
truedag_dataset$HYPOVOLEMIA <- as.factor(truedag_dataset$HYPOVOLEMIA)
truedag_dataset$HISTORY <- as.factor(truedag_dataset$HISTORY)
truedag_dataset$ANAPHYLAXIS <- as.factor(truedag_dataset$ANAPHYLAXIS)
truedag_dataset$HISTORY <- as.factor(truedag_dataset$HISTORY)
truedag_dataset$INSUFFANESTH <- as.factor(truedag_dataset$INSUFFANESTH)
truedag_dataset$PULMEMBOLUS <- as.factor(truedag_dataset$PULMEMBOLUS)
truedag_dataset <- truedag_dataset[,order(colnames(truedag_dataset))]



modelstring = paste0("[HISTORY|LVFAILURE][CVP|LVEDVOLUME][PCWP|LVEDVOLUME][HYPOVOLEMIA][LVEDVOLUME|HYPOVOLEMIA:LVFAILURE][LVFAILURE]",
                     "[STROKEVOLUME|HYPOVOLEMIA:LVFAILURE][ERRLOWOUTPUT][HRBP|ERRLOWOUTPUT:HR][HREKG|ERRCAUTER:HR][ERRCAUTER][HRSAT|ERRCAUTER:HR][INSUFFANESTH]",
                     "[ANAPHYLAXIS][TPR|ANAPHYLAXIS][EXPCO2|ARTCO2:VENTLUNG][KINKEDTUBE][MINVOL|INTUBATION:VENTLUNG][FIO2][PVSAT|FIO2:VENTALV]",
                     "[SAO2|PVSAT:SHUNT][PAP|PULMEMBOLUS][PULMEMBOLUS][SHUNT|INTUBATION:PULMEMBOLUS][INTUBATION][PRESS|INTUBATION:KINKEDTUBE:VENTTUBE][DISCONNECT]",
                     "[MINVOLSET][VENTMACH|MINVOLSET][VENTTUBE|DISCONNECT:VENTMACH][VENTLUNG|INTUBATION:KINKEDTUBE:VENTTUBE][VENTALV|INTUBATION:VENTLUNG]",
                     "[ARTCO2|VENTALV][CATECHOL|ARTCO2:INSUFFANESTH:SAO2:TPR][HR|CATECHOL][CO|HR:STROKEVOLUME][BP|CO:TPR]")
true_graph = model2network(modelstring)


true_mag = empty.graph(c("VENTLUNG","MINVOLSET","DISCONNECT","VENTMACH","PRESS","FIO2","MINVOL","KINKEDTUBE","STROKEVOLUME","LVEDVOLUME","HYPOVOLEMIA","VENTALV","CATECHOL"
                         ,"ARTCO2","HR","CO","BP","PVSAT","EXPCO2","HREKG","ERRCAUTER","HRBP","ERRLOWOUTPUT","HRSAT"
                         ,"HISTORY","ANAPHYLAXIS","PAP","INTUBATION","VENTTUBE","TPR","PCWP","INSUFFANESTH","PULMEMBOLUS","CVP","SAO2"))


arc.set = matrix(c("ANAPHYLAXIS","TPR","ARTCO2","CATECHOL","ARTCO2","EXPCO2","CATECHOL","HR","CO","BP","DISCONNECT","VENTTUBE","ERRCAUTER","HREKG"
                   ,"ERRCAUTER","HRSAT","ERRLOWOUTPUT","HRBP","FIO2","PVSAT","HR","CO","HR","HRBP","HR","HREKG","HR","HRSAT"
                   ,"HYPOVOLEMIA","LVEDVOLUME","HYPOVOLEMIA","STROKEVOLUME","INSUFFANESTH","CATECHOL","INTUBATION","MINVOL","INTUBATION","PRESS"
                   ,"INTUBATION","VENTALV","INTUBATION","VENTLUNG","KINKEDTUBE","PRESS","KINKEDTUBE","VENTLUNG","LVEDVOLUME","CVP"
                   ,"LVEDVOLUME","PCWP","MINVOLSET","VENTMACH","PULMEMBOLUS","PAP","PVSAT","SAO2","SAO2","CATECHOL","STROKEVOLUME","CO","TPR","BP","TPR","CATECHOL"
                   ,"VENTALV","ARTCO2","VENTALV","PVSAT","VENTLUNG","EXPCO2","VENTLUNG","MINVOL","VENTLUNG","VENTALV"
                   ,"VENTMACH","VENTTUBE","VENTTUBE","PRESS","VENTTUBE","VENTLUNG","HISTORY","LVEDVOLUME"
                   ,"HISTORY","STROKEVOLUME","LVEDVOLUME","STROKEVOLUME","INTUBATION","SAO2","PULMEMBOLUS","SAO2","LVEDVOLUME","HISTORY","STROKEVOLUME","HISTORY","STROKEVOLUME","LVEDVOLUME"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))

TetradGetAdjmat <- function(tetradrunner) {
  p <- length(tetradrunner$nodes)
  adjmat <- matrix(0, p, p)
  for (e in edges) {
    edgevec <- unlist(strsplit(e, " "))
    i <- match(edgevec[1], tetradrunner$nodes)
    j <- match(edgevec[3], tetradrunner$nodes)
    if (edgevec[2] == '---') {
      edge <- c(edgevec[1], edgevec[3])
      adjmat[i, j] <- 1
      adjmat[j, i] <- 1
    } else if (edgevec[2] == "-->") {
      adjmat[i, j] <- 1
    }
  }
  return(adjmat)
}



arcs(true_mag) = arc.set
n <- ncol(dataset)
true_mag_amat <- amat(true_mag)

local_BIC <- matrix(0, nrow = n+1, ncol = n_iter+2)


for(iteration in n_iter)
{
  #type == 'single'
  # In each experiment there is n (target) intervention on each variable,
  # except the last experiment which has no interventions.
  
  
  #nexp <- sample(1:n, 1)+1
  nexp <-iteration
  m <- matrix(0, nexp , n)
  round <- list()
  tic.clearlog()
  
  for(i in 1:nexp)
  {
    round[[i]]<- cbind(rep(i, each=target), sample(1:n, target))
    #round < rbind(round,round[i])
  }
  
  target_matrix <- round[[1]]
  for(i in 2:nexp)
  {
    
    target_matrix <- rbind(target_matrix,round[[i]])
  }
  m[target_matrix] <- 1
  E <- rbind(m,rep(0,n))
  
  #####################################################################
  
  Data<- list()  
  
  tic("Tetradtime fges")
  for ( i in 1:nrow(E)) {
    # Create the tuple that will be stored in D, storing the vector e.
    Data[[i]]<-list(e=E[i,])
    # Get indexes of intervened variables in this particular experimental setting.
    J<-which( Data[[i]]$e==1 )
    # The data consist of manipulated graphs where
    # Edge heads into the intervened variables are cut.
    
    Data[[i]]$data <- truedag_dataset
    interventiontarget_dataset <- dataset[,order(colnames(dataset))]
    #skeleton = hc(Data[[i]]$data,debug=FALSE)
    
    if (length(J) !=0)
    {
      list_to_remove <- character(0)
      for(index in J)
      {
        list_to_remove_tempt <- incoming.arcs(true_graph, toString(data.frame(colnames(interventiontarget_dataset))[index,1]))
        list_to_remove <- rbind(list_to_remove_tempt,list_to_remove)
      }
      
    }
    
    
    if (length(J) !=0 & nrow(list_to_remove) !=0 )
    {
      skeleton <- true_graph
      for ( incoming_index in 1: nrow(list_to_remove)) {
        skeleton <- drop.arc(skeleton, list_to_remove[incoming_index,1], list_to_remove[incoming_index,2])
      }
      
    }
    else{
      invertentional_dataset <- dataset
      param_skeleton <- NULL
    }
    
    if (length(J) !=0)
    {
      
      param_skeleton <- bn.fit(skeleton, Data[[i]]$data, method = "bayes", keep.fitted = TRUE, debug = TRUE)
      
      
      nms <- names(interventiontarget_dataset[J])
      for(index_target in 1:length(nms))
      {
        T=param_skeleton[[nms[index_target]]]$prob
        
        uniform <- 1/length(T)
        for (index in 1:length(T))
        {
          
          T[[index]]=uniform
          
          # if (index ==2)
          # {
          #   T[[index]]=1
          # }
          # else
          # {
          #   T[[index]]=0
          # }
        }
        
        
        param_skeleton[[nms[index_target]]] =T
        
        
      }
      
      
      
      interventional_dataset <- rbn(param_skeleton, n = 10000, debug = FALSE)
      interventional_dataset$LVFAILURE <- NULL
      interventional_dataset$SHUNT <- NULL
      
      Data[[i]]$data <- interventional_dataset
      Data[[i]]$param <- param_skeleton
   
      tetradrunner <- tetradrunner(algoId = 'fges',
                                   df = Data[[i]]$data,
                                   scoreId = 'bic',
                                   dataType = 'discrete',
                                   faithfulnessAssumed = TRUE,
                                   maxDegree = -1, #argv$maxDegree,
                                   verbose = FALSE)



      Data[[i]]$pag <-tetradrunner
    }
    else
    {
      Data[[i]]$data <- dataset
    
      tetradrunner <- tetradrunner(algoId = 'fges',
                                   df = Data[[i]]$data,
                                   scoreId = 'bic',
                                   dataType = 'discrete',
                                   faithfulnessAssumed = TRUE,
                                   maxDegree = -1, #argv$maxDegree,
                                   verbose = FALSE)

      Data[[i]]$pag <-tetradrunner
      
      
     
    }
    # index the local BIC for node
    print(paste("iteration",iteration))
    nodesx <-tetradrunner$graph
    nodesx$getAttribute('BIC')
    
    nodes <- nodesx$getNodes()
    
    local_BIC[n+1,1]<- 'model score'
    for(i_bic in 0:as.integer(nodes$size()-1)){
      
      local_BIC[i_bic+1,1] <-nodes$get(as.integer(i_bic))$getName()
      local_BIC[i_bic+1,i+1] <-nodes$get(as.integer(i_bic))$getAttribute('BIC')
      local_BIC[n+1,i+1]<- nodesx$getAttribute('BIC')
      
    
    
    
  
  }

    
    
    
  
   }
  
  toc(log = TRUE, quiet = TRUE)
  log.txt <- tic.log(format = TRUE)
  
  
    amat<- matrix(0, n, n)
    for ( i in 1:nrow(E)) {
      
      #if (i==2||i==3 ||i==4||i==5||i==9)
      
      {
        edges <- Data[[i]]$pag$edges
        amat_temp <- TetradGetAdjmat(Data[[i]]$pag)
        
        

        rownames(amat_temp) <- colnames(dataset)
        colnames(amat_temp) <- colnames(dataset)
        
        
        amat_temp <- amat_temp[order(rownames(amat_temp)),order(colnames(amat_temp))]
        amat_temp [ amat_temp ==t(amat_temp)  ] <- 0
        #amat <- round(amat/nrow(E), digits=2)
        #amat_temp <- Data[[i]]$pag@amat
        # if (i==nrow(E))
        # {
        #   amat <- amat+nrow(E)*amat_temp+amat_temp
        # }
        # else{
        amat <- amat+amat_temp
        # }
        
        
        
      }
      
      
    }  
  
  
    amat <- (amat/(nrow(E)))
    pheatmap(t(amat), cluster_rows = FALSE, cluster_cols = FALSE,color = colorRampPalette(c("white","grey", "blue"))(100))
    
    #amat <- round(amat)
    amat [ amat <thred ] <- 0
    amat [ amat >thred] <- 1
    #print(mmhc(dataset))
    
    #pheatmap(t(true_mag_amat), cluster_rows = FALSE, cluster_cols = FALSE,color = colorRampPalette(c("white","green"))(100))
    
    
   # set skeleton symmetry
    ske_amat <- amat
    ske_amat [ ske_amat != t(ske_amat)] <- 1

    skeleton <- new("pcAlgo")
    skeleton1 <- new("graphAM",ske_amat)
    skeleton2 <- as(skeleton1, "graphNEL")
    skeleton@graph <-skeleton2
    dataset_pcalg<- Data[[i]]$data
    cols <- colnames(dataset_pcalg)
    for(j in cols){
      dataset_pcalg[[j]]<- as.numeric(dataset_pcalg[[j]])
      dataset_pcalg[[j]]<- dataset_pcalg[[j]]-1
    }
    datamatrix <- data.matrix(na.omit(dataset_pcalg))
    i <- 3
    varNames <- colnames(Data[[i]]$data)
    suffStat = list(dm = datamatrix, nlev = sapply(Data[[i]]$data , nlevels),adaptDF = FALSE)

    skel.fit <- pcalg::skeleton(suffStat, indepTest=disCItest, label=c(varNames),alpha=0.01,  verbose = FALSE)

    p <-list()
    q <-list()
    for(k in 1:n)
    {
      p[[k]] = integer()


    }
      for(l in 1:n)
      {
        q[[l]] <- p
      }
    
    skeleton@sepset <- q
    source("pc.cons.intern2.R")
    rule0_graph <- pc.cons.intern2(skeleton, suffStat, indepTest=disCItest, alpha=0.01, version.unf = c(1, 1),maj.rule = TRUE, verbose = TRUE)
    #varNames <- colnames(Data[[i]]$data)
    #pag <- fci(suffStat, indepTest=disCItest,rules = rep(TRUE,10),alpha = 0.01,label=c(varNames))

   unsh <- find.unsh.triple(ske_amat, check=TRUE)
   vstruct <- unsh$unshTripl
  
   true_mag_amat <- amat(true_mag) 
   true_mag_amat<- true_mag_amat[order(as.character(rownames(true_mag_amat))), order(as.character(colnames(true_mag_amat)))]
  #add bi-directed edge
   true_mag_amat[21,12] <-1
   true_mag_amat[30,12] <-1
   true_mag_amat[30,21] <-1
   header <- 1:n
   rownames(true_mag_amat) <-c(header)
   colnames(true_mag_amat) <-c(header)

   x <- as.bn(as(true_mag_amat, "graphNEL"))
   v <- vstructs(x)
   v <- t(v)
   matrix(logical(), nrow = ncol(vstruct), ncol = 1)
   vstrucr_true <- rbind(vstruct,  matrix(logical(), nrow = 1, ncol = ncol(vstruct)))



   for(index_sepset_i in 1:n)
     {
     for(index_sepset_j in 1:n)
     {
       size <- length(rule0_graph[["sk"]]@sepset[[index_sepset_i]][[index_sepset_j]])

      

               for(vstruct_true_index in 1:ncol(vstrucr_true))
                 {
                 
                         if (index_sepset_i==vstrucr_true[1,vstruct_true_index] && index_sepset_j==vstrucr_true[3,vstruct_true_index] )
                         {
                           
                           if (size !=0)
                           {
                                 count_maj <-0
                                       for(sepset_index in 1:size)
                                       {
                                         if(rule0_graph[["sk"]]@sepset[[index_sepset_i]][[index_sepset_j]][sepset_index]!=vstrucr_true[2,vstruct_true_index] )
                                         {
                                           count_maj <- count_maj+1
                                           print(rule0_graph[["sk"]]@sepset[[index_sepset_i]][[index_sepset_j]][sepset_index])
                          
                                         }
                          
                                       }
                                 if (count_maj==size)
                                 {
                                   vstrucr_true[4,vstruct_true_index] <-TRUE
                                 }
                                 else if (count_maj==0)
                                 {
                                   vstrucr_true[4,vstruct_true_index] <-FALSE
                                 }
              
                           }
                           else if (size ==0)
                           {
                             vstrucr_true[4,vstruct_true_index] <-TRUE
                           }
      
                }

                
       }
       
       
     }

}










   tp <- 0 
   ans <- 0
   for(index_truev in 1:ncol(v))
   {
     for(index_v in 1:ncol(vstrucr_true))
     {
       if (is.na(vstrucr_true[4,index_v])==FALSE && vstrucr_true[4,index_v]==TRUE&&v[1,index_truev]==vstrucr_true[1,index_v] &&v[2,index_truev]==vstrucr_true[2,index_v]&&v[3,index_truev]==vstrucr_true[3,index_v])
       {
         tp <- tp+1
       }
      
       if(is.na(vstrucr_true[4,index_v]) ==FALSE && vstrucr_true[4,index_v]==TRUE)
       {
         ans <- ans+1
       }
     }
     
   }
   precition_v <-0
   recall_v <-0
   ans <- ans/ncol(v)
   precition_v <-  tp/ans
   recall_v <-  tp/ncol(v)






   
  # which(length(rule0_graph$sk@sepset[[2]]!=0))
  
  source("precisionrecall.R")
  true_mag_amat <- amat(true_mag)
  true_mag_amat<- true_mag_amat[order(rownames(true_mag_amat)),order(colnames(true_mag_amat))]
  #amat_round <- round(amat)
  k <-precisionrecall(amat, true_mag_amat)
  iterFGES$nExp[[iteration]] <- nexp
  iterFGES$precision[[iteration]] <- k[["precision"]] 
  iterFGES$recall[[iteration]] <- k[["recall"]]
    
  iterFGES$BSF[[iteration]] <- k[["BSF"]]
  iterFGES$time[[iteration]] <- log.txt[[1]]
  C <- matrix(0, nrow=1, ncol=35)
  
  E <- rbind(C, E)
  E <- E[-nrow(E),] 
  
  L1 <- matrix(0, nrow=iteration+1, ncol=1)
  E1 <- cbind(E[,1:21],L1)
  E2 <- cbind(E1,E[,22:35])
  
  
  E3 <- cbind(E2[,1:30],L1)
  E4 <- cbind(E3,E2[,31:36])
  
  
  
  set4 <- "nExp.csv"
  path <-"TestINT/"
  write.csv(E4,paste0(path,iteration,set4))
  
}
# library(DOT)
# 
# graph_dot <- tetradrunner.tetradGraphToDot(Data[[1]][["pag"]][["graph"]])
# dot(graph_dot)

  # set1 <- "precision.csv"
  # set2 <- "recall.csv"
  # set3 <- "BSF.csv"
  # set5 <- "time.csv"
  # write.csv(iterFGES$precision,paste0(path,set1))
  # write.csv(iterFGES$recall,paste0(path,set2))
  # write.csv(iterFGES$BSF,paste0(path,set3))
  # write.csv(iterFGES$time,paste0(path,set5))
  


sink()
  # write.csv(true_mag_amat,"mag.csv")