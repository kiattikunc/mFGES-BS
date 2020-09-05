#alarm n data sets with input targets

library(rcausal)
library(DOT)
library(tictoc)
library(pcalg)
library(bnlearn)
library(pheatmap)
library(tictoc)


#load dataset

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
                   ,"HISTORY","STROKEVOLUME","LVEDVOLUME","STROKEVOLUME","INTUBATION","SAO2","PULMEMBOLUS","SAO2"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
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

true_mag_amat <- amat(true_mag)

iterFGES <- list()
thred <- 0.2
target <-5
set.seed(5)
for(iteration in 8)
{
  #type == 'single'
  # In each experiment there is n (target) intervention on each variable,
  # except the last experiment which has no interventions.
  
  n <- ncol(dataset)
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
    
    print(paste("iteration",iteration))
    
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
  path <-"/Users/kiattikunchobtham/Documents/GitHub/R-structure-learning-with-latent/TestINT/"
  write.csv(E4,paste0(path,iteration,set4))
}
  
  

# set1 <- "precision.csv"
# set2 <- "recall.csv"
# set3 <- "BSF.csv"
# set5 <- "time.csv"
# write.csv(iterFGES$precision,paste0(path,set1))
# write.csv(iterFGES$recall,paste0(path,set2))
# write.csv(iterFGES$BSF,paste0(path,set3))
# write.csv(iterFGES$time,paste0(path,set5))


# write.csv(true_mag_amat,"/Users/kiattikunchobtham/Documents/GitHub/R-structure-learning-with-latent/mag.csv")
