#alarm n data sets with input targets
# version 7 to integrate with bayesian prob with depenc or indepen with uniform
#sink(file='myoutput.txt')
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
n_iter <-35
n_intervent <- 10000


set <-'alarm_L5_full_10k.csv'
dataset2 <- read.csv(set,header = FALSE,na.strings=c(""))

#load observational dataset
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



#load N dataset for CPT under true DAG.
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


#construct true DAG
modelstring = paste0("[HISTORY|LVFAILURE][CVP|LVEDVOLUME][PCWP|LVEDVOLUME][HYPOVOLEMIA][LVEDVOLUME|HYPOVOLEMIA:LVFAILURE][LVFAILURE]",
                     "[STROKEVOLUME|HYPOVOLEMIA:LVFAILURE][ERRLOWOUTPUT][HRBP|ERRLOWOUTPUT:HR][HREKG|ERRCAUTER:HR][ERRCAUTER][HRSAT|ERRCAUTER:HR][INSUFFANESTH]",
                     "[ANAPHYLAXIS][TPR|ANAPHYLAXIS][EXPCO2|ARTCO2:VENTLUNG][KINKEDTUBE][MINVOL|INTUBATION:VENTLUNG][FIO2][PVSAT|FIO2:VENTALV]",
                     "[SAO2|PVSAT:SHUNT][PAP|PULMEMBOLUS][PULMEMBOLUS][SHUNT|INTUBATION:PULMEMBOLUS][INTUBATION][PRESS|INTUBATION:KINKEDTUBE:VENTTUBE][DISCONNECT]",
                     "[MINVOLSET][VENTMACH|MINVOLSET][VENTTUBE|DISCONNECT:VENTMACH][VENTLUNG|INTUBATION:KINKEDTUBE:VENTTUBE][VENTALV|INTUBATION:VENTLUNG]",
                     "[ARTCO2|VENTALV][CATECHOL|ARTCO2:INSUFFANESTH:SAO2:TPR][HR|CATECHOL][CO|HR:STROKEVOLUME][BP|CO:TPR]")
true_graph = model2network(modelstring)

#construct true MAG
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

#function to convert tetrad format to Matrix format 
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

#add bi-directed edge to true MAG
true_mag_amat <- amat(true_mag) 
true_mag_amat<- true_mag_amat[order(as.character(rownames(true_mag_amat))), order(as.character(colnames(true_mag_amat)))]
true_mag_amat[21,12] <-1
true_mag_amat[30,12] <-1
true_mag_amat[30,21] <-1
####



local_BIC <- matrix(0, nrow = n+1, ncol = n_iter+2)




#type == 'single'
# In each experiment there is n (target) intervention on each variable,
# except the last experiment which has no interventions.
for(iteration in 2:n_iter)
{
  #nexp <- sample(1:n, 1)+1
  nexp <-iteration
  m <- matrix(0, nexp , n)
  round <- list()
  tic.clearlog()
  
  for(ir in 1:nexp)
  {
    round[[ir]]<- cbind(rep(ir, each=target), sample(1:n, target))
    #round < rbind(round,round[i])
  }
  
  target_matrix <- round[[1]]
  for(is in 1:nexp)
  {
    
    target_matrix <- rbind(target_matrix,round[[is]])
  }
  m[target_matrix] <- 1
  E <- rbind(m,rep(0,n))
  
  #####################################################################
  
  
  
  Data<- list()  
  
  #generate n interventional dataset
  
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
      inverventional_dataset <- dataset
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
      
      
      #generate interventioanl dataset
      interventional_dataset <- rbn(param_skeleton, n = n_intervent, debug = FALSE)
      
      #export int-data
      cols_int_all <- colnames(interventional_dataset)
      for(cols_int in cols_int_all){
        interventional_dataset[[cols_int]]<- as.numeric(interventional_dataset[[cols_int]])
        interventional_dataset[[cols_int]]<- interventional_dataset[[cols_int]]-1
      }
      set_int <- paste0('Data/',iteration,"-",i,"-alarm_N_10k.csv")
      
      write.table(interventional_dataset,paste0(set_int),row.names = FALSE, col.names=FALSE,sep=",")
  
      
    }
    else
    {
    
      Data[[i]]$data <- dataset
      set_int <- paste0('Data/',iteration,"-",i,"-alarm_N_10k.csv")
      write.table(dataset2,paste0(set_int),row.names = FALSE, col.names=FALSE,sep=",")
      
      
      
    }

    
  }

}

