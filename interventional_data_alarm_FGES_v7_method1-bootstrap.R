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
dataset2 <- dataset[,order(colnames(dataset))]

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
      adjmat[i, j] <- 3
      adjmat[j, i] <- 3
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



bn.boot(data, statistic, R = 200, m = nrow(data), algorithm,
        algorithm.args = list(), statistic.args = list(), cluster = NULL,
        debug = FALSE)





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
     
      # #export int-data
      # cols_int_all <- colnames(interventional_dataset)
      # for(cols_int in cols_int_all){
      #   interventional_dataset[[cols_int]]<- as.numeric(interventional_dataset[[cols_int]])
      #   interventional_dataset[[cols_int]]<- interventional_dataset[[cols_int]]-1
      # }
      # set_int <- paste0('Data/',iteration,"-",i,"-alarm_N_10k.csv")
      # 
      # write.table(interventional_dataset,paste0(set_int),row.names = FALSE, col.names=FALSE,sep=",")
      # 
   
      #remove latent variables
      interventional_dataset$LVFAILURE <- NULL
      interventional_dataset$SHUNT <- NULL

      
      Data[[i]]$data <- interventional_dataset
      Data[[i]]$param <- param_skeleton
   
     
    }
    else
    {
      Data[[i]]$data <- dataset
    }
   
    print(paste("iteration",iteration))
  #   nodesx <-tetradrunner$graph
  #  # nodesx$getAttribute('BIC')
  #  
  # 
  #   nodes <- nodesx$getNodes()
  #   
  #   local_BIC[n+1,1]<- 'model score'
  #   for(i_bic in 0:as.integer(nodes$size()-1)){
  #     
  #     local_BIC[i_bic+1,1] <-nodes$get(as.integer(i_bic))$getName()
  #     local_BIC[i_bic+1,i+1] <-nodes$get(as.integer(i_bic))$getAttribute('BIC')
  #     local_BIC[n+1,i+1]<- nodesx$getAttribute('BIC')
  #     
  #   
  # 
  # }

  
   }
  
  ############## START #########################
  tic("Tetradtime fges")
  
  for ( data_j in 1:nrow(E)) 
    {
        tetradrunner <- tetradrunner(algoId = 'fges',
                               df = Data[[data_j]]$data,
                               scoreId = 'bdeu',
                               dataType = 'discrete',
                               faithfulnessAssumed = TRUE,
                               maxDegree = -1, #argv$maxDegree,
                               verbose = FALSE)

  
        Data[[data_j]]$pag <-tetradrunner
    }
  

  # plain model averaging concept from frequentist approach
  amat<- matrix(0, n, n)
  amat_undirect <- matrix(0, n, n)
  
  for ( i in 1:nrow(E)) {
    
    #if (i==2||i==3 ||i==4||i==5||i==9)
    
    { 
      edges <- Data[[i]]$pag$edges
     
      amat_temp <- TetradGetAdjmat(Data[[i]]$pag)
      amat_temp2 <- matrix(0, n, n) # for undirected edge for intv dataset
      rownames(amat_temp) <- colnames(dataset)
      colnames(amat_temp) <- colnames(dataset)
      
      ##split undirect from matrix
      #amat_temp <- amat_temp[order(rownames(amat_temp)),order(colnames(amat_temp))]
      amat_temp2 [ amat_temp ==t(amat_temp) & amat_temp ==3  ] <- 0.5
      amat_temp [ amat_temp ==t(amat_temp)  ] <- 0
      
      #amat <- round(amat/nrow(E), digits=2)
      #amat_temp <- Data[[i]]$pag@amat
      # if (i==nrow(E))
      # {
      #   amat <- amat+nrow(E)*amat_temp+amat_temp
      # }
      # else{
      amat <- amat+amat_temp + amat_temp2
      
  
    }
    
  }  
  amat <- (amat/(nrow(E)))
  amat_temp <- TetradGetAdjmat(Data[[i]]$pag)
  amat_undirect [ amat_temp ==t(amat_temp) &  amat_temp ==3 ] <- 3
  
    #pheatmap(t(amat), cluster_rows = FALSE, cluster_cols = FALSE,color = colorRampPalette(c("white","grey", "blue"))(100))
    
    #amat <- round(amat)
    # amat [ amat <thred ] <- 0
    # amat [ amat >thred] <- 1
    #print(mmhc(dataset))
    
    #pheatmap(t(true_mag_amat), cluster_rows = FALSE, cluster_cols = FALSE,color = colorRampPalette(c("white","green"))(100))
    
    
    

    
    
    
    
## BOX2: use the skeleton from essntial graph from BOX1 and perform the statistic test to identify v-structure using conservative rule with majority rule

#     #edges <- Data[[iteration+1]]$pag$edges
#     #essential_amat <- TetradGetAdjmat(Data[[iteration+1]]$pag)
#     essential_amat <- amat
#     essential_amat [ essential_amat <thred ] <- 0
#     essential_amat [ essential_amat >thred] <- 1
# 
# 
#     rownames(essential_amat) <- colnames(dataset)
#     colnames(essential_amat) <- colnames(dataset)
# 
#     ske_amat <-essential_amat
#     ske_amat [ ske_amat != t(ske_amat)] <- 1  # create skeleton as undirected graph
#     #amat [ amat != t(amat)] <- 1
#     skeleton <- new("pcAlgo")
#     skeleton1 <- new("graphAM",ske_amat)
#     skeleton2 <- as(skeleton1, "graphNEL") # convert amat format to pcalg graph format
#     skeleton@graph <-skeleton2
# 
#     obs_index <- i # use observational dataset to identify skeleton and v-structure
#     varNames <- colnames(Data[[obs_index]]$data)
# 
#    #df<-Data[[1]]$data
#     df<-Data[[obs_index]]$dat
#     #### combine all dataset to identify v-strucure
#     # for(data_index in 2:obs_index)
#     # {
#     #   df <- rbind(df, Data[[data_index]]$data)
#     # }
#     ##########
#     dataset_pcalg<- df
#     cols <- colnames(dataset_pcalg)
#     for(j in cols){
#       dataset_pcalg[[j]]<- as.numeric(dataset_pcalg[[j]])
#       dataset_pcalg[[j]]<- dataset_pcalg[[j]]-1
#     }
#     datamatrix <- data.matrix(na.omit(dataset_pcalg)) # convert rawdata from category to level 0, 1,2,..
# 
#     suffStat = list(dm = dataset_pcalg, nlev = sapply(Data[[obs_index]]$data , nlevels),adaptDF = FALSE)
#     #skel.fit <- pcalg::skeleton(suffStat, indepTest=disCItest, label=c(varNames),alpha=0.01,  verbose = FALSE)
# 
#     p <-list()
#     q <-list()
#     for(k in 1:n)
#     {
#       p[[k]] = integer()
#     }
#       for(l in 1:n)
#       {
#         q[[l]] <- p
#       }
# 
# 
#     skeleton@sepset <- q
#     source("pc.cons.intern2.R")
# 
# 
# 
#     rule0_graph <- pc.cons.intern2(skeleton, suffStat, indepTest=disCItest, alpha=0.01, version.unf = c(1, 1),maj.rule = TRUE, verbose = FALSE)
#     #varNames <- colnames(Data[[i]]$data)
#     #pag <- fci(suffStat, indepTest=disCItest,rules = rep(TRUE,10),alpha = 0.01,label=c(varNames))
# 
# 
# 
#     unsh <- find.unsh.triple(ske_amat, check=TRUE)
#     # unsh <- find.unsh.triple(amat, check=TRUE)
#     vstruct <- unsh$unshTripl
# 
# 
# 
# 
# # 
#    header <- 1:n
#    rownames(ske_amat) <-c(header)
#    colnames(ske_amat) <-c(header)
# 
# 
#    matrix(logical(), nrow = ncol(vstruct), ncol = 1)
#    vstrucr_true <- rbind(vstruct,  matrix(logical(), nrow = 1, ncol = ncol(vstruct)))
# 
#    # mark which unshield triple to v-structue using conservative rule with majority rule
#    for(index_sepset_i in 1:n)
#      {
#      for(index_sepset_j in 1:n)
#      {
#        size <- length(rule0_graph[["sk"]]@sepset[[index_sepset_i]][[index_sepset_j]])
# 
#                for(vstruct_true_index in 1:ncol(vstrucr_true))
#                  {
# 
#                          if (index_sepset_i==vstrucr_true[1,vstruct_true_index] && index_sepset_j==vstrucr_true[3,vstruct_true_index] )
#                          {
# 
#                            if (size !=0)
#                            {
#                              ratio <- rule0_graph[["sepsetratio"]][[index_sepset_i]][[index_sepset_j]]
#                              if(length(ratio)>0)
#                              {
#                                if(ratio>0.5)
#                                {
#                                  vstrucr_true[4,vstruct_true_index] <-FALSE
#                                }
# 
# 
#                                else if (ratio<0.5)
#                                {
#                                  vstrucr_true[4,vstruct_true_index] <-TRUE
#                                }
#                                else
#                                {
#                                  vstrucr_true[4,vstruct_true_index] <- 'N/A'
#                                }
#                              }
# 
#                            }
#                            else if (size ==0)
#                            {
# 
#                                 vstrucr_true[4,vstruct_true_index] <- TRUE
# 
# 
#                            }
# 
#                 }
# 
# 
#        }
# 
# 
#      }
# 
# }



#soft contrain v-structue integrated to prior

#     
# for(index_sepset_i in 1:n)
#      {
#       for(index_sepset_j in 1:n)
#        {
#           size <- length(rule0_graph[["sk"]]@sepset[[index_sepset_i]][[index_sepset_j]])
# 
#           if (size ==0)
#           {
#             a<- which(vstruct[1,] == index_sepset_i & vstruct[3,] == index_sepset_j , arr.ind = TRUE)
#             if (length(a) !=0)
#             {
# 
#               if (length(a) >1)
#               {
#                 for (a_index in 1: length(a))
#                 {
# 
#                   collider<- vstruct[[2,a[a_index]]]
#                    #ratio <-1-rule0_graph[["sepsetratio"]][[index_sepset_i]][[index_sepset_j]]
# 
#                   ratio <-max(1-rule0_graph[["sepsetratio"]][[index_sepset_i]][[index_sepset_j]],0.5)
#                   if (ratio>0.5)
#                   {
#                     print(paste0(colnames(amat)[index_sepset_i],'-',colnames(amat)[collider],'-',colnames(amat)[index_sepset_j], ' ratio=',ratio))
#                     
#                   }
#                   amat[index_sepset_i,collider] <- max(amat[index_sepset_i,collider],ratio)
#                   amat[index_sepset_j,collider] <- max(amat[index_sepset_j,collider],ratio)
# 
#                 }
#               }
# 
# 
#             }
# 
#           }
#           else
#           {
#             a<- which(vstruct[1,] == index_sepset_i & vstruct[3,] == index_sepset_j , arr.ind = TRUE)
#             if (length(a) !=0)
#             {
#               if (length(a) >1)
#               {
#                 for (a_index in 1: length(a))
#                 {
# 
#                   collider<- vstruct[[2,a[a_index]]]
#                    #ratio <-1-rule0_graph[["sepsetratio"]][[index_sepset_i]][[index_sepset_j]]
#                   ratio <-max(1-rule0_graph[["sepsetratio"]][[index_sepset_i]][[index_sepset_j]],0.5)
#                   if (ratio>0.5)
#                   {
#                     print(paste0(colnames(amat)[index_sepset_i],'-',colnames(amat)[collider],'-',colnames(amat)[index_sepset_j], ' ratio=',ratio))
#                     
#                   }
#                   amat[index_sepset_i,collider] <- max(amat[index_sepset_i,collider],ratio)
#                   amat[index_sepset_j,collider] <- max(amat[index_sepset_j,collider],ratio)
# 
#                 }
#               }
#               else
#               {
#                 collider<- vstruct[[2,a]]
#                 #ratio <-1-rule0_graph[["sepsetratio"]][[index_sepset_i]][[index_sepset_j]]
#                 ratio <-max(1-rule0_graph[["sepsetratio"]][[index_sepset_i]][[index_sepset_j]],0.5)
#                 if (ratio>0.5)
#                 {
#                   print(paste0(colnames(amat)[index_sepset_i],'-',colnames(amat)[collider],'-',colnames(amat)[index_sepset_j], ' ratio=',ratio))
#                   
#                 }
#                 amat[index_sepset_i,collider] <- max(amat[index_sepset_i,collider],ratio)
#                 amat[index_sepset_j,collider] <- max(amat[index_sepset_j,collider],ratio)
#               }
# 
# 
#             }
# 
# 
#           }
# 
# 
#         }
#       }




   #check the performance of conservative-rule
   # true v-structure 
   # x <- as.bn(as(true_mag_amat, "graphNEL"))
   # v <- vstructs(x)
   # v <- t(v)
   # 
   # 
   # tp <- 0
   # ans <- 0
   # for(index_truev in 1:ncol(v))
   # {
   #   for(index_v in 1:ncol(vstrucr_true))
   #   {
   #     if (is.na(vstrucr_true[4,index_v])==FALSE && vstrucr_true[4,index_v]==TRUE&&v[1,index_truev]==vstrucr_true[1,index_v] &&v[2,index_truev]==vstrucr_true[2,index_v]&&v[3,index_truev]==vstrucr_true[3,index_v])
   #     {
   #       tp <- tp+1
   #     }
   # 
   #     if(is.na(vstrucr_true[4,index_v]) ==FALSE && vstrucr_true[4,index_v]==TRUE)
   #     {
   #       ans <- ans+1
   #     }
   #   }
   # 
   # }
   # 
   # 
   # v_precition <-0
   # v_recall <-0
   # ans <- ans/ncol(v)
   # v_precition <-  tp/ans
   # v_recall <-  tp/ncol(v)

#############################################################################

# catch Bdeu score for constraints
  
   indep_graphs = empty.graph(c("VENTLUNG","MINVOLSET","DISCONNECT","VENTMACH","PRESS","FIO2","MINVOL","KINKEDTUBE","STROKEVOLUME","LVEDVOLUME","HYPOVOLEMIA","VENTALV","CATECHOL"
                                ,"ARTCO2","HR","CO","BP","PVSAT","EXPCO2","HREKG","ERRCAUTER","HRBP","ERRLOWOUTPUT","HRSAT"
                                ,"HISTORY","ANAPHYLAXIS","PAP","INTUBATION","VENTTUBE","TPR","PCWP","INSUFFANESTH","PULMEMBOLUS","CVP","SAO2"))
   
   source("posteriorcal.R")
    #method1 - prior/N_data
   
    #method2 - prior/N_data*3
  
   posterior<- posteriorcal(amat,indep_graphs,Data,n)
   posterior_obs <-posterior$indepence_posterior_obs
   
   
   
posterior_obs<- 1-posterior_obs
posterior_obs [ posterior_obs >thred] <- 1
posterior_obs [ posterior_obs <thred ] <- 0
diag(posterior_obs) <- 0

learned_mag <- matrix(0, n, n)
learned_mag[posterior_obs==1] <- 1
#learned_mag[posterior_obs==1 & amat_undirect==0 ] <- 1
#learned_mag[posterior_obs==1 & amat_undirect==3 ] <- 1
learned_mag[posterior_obs==0 & t(posterior_obs)==0 & amat_undirect==3 ] <- 3


############## STOP #########################
toc(log = TRUE, quiet = TRUE)
log.txt <- tic.log(format = TRUE)
#############################################

  # which(length(rule0_graph$sk@sepset[[2]]!=0))
  
  source("precisionrecall_undirected_edge.R")
  true_mag_amat <- amat(true_mag)
  
  
  
  
 
  true_mag_amat<- true_mag_amat[order(rownames(true_mag_amat)),order(colnames(true_mag_amat))]

  k <-precisionrecall(learned_mag, true_mag_amat)
 

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
  


# library(DOT)
# 
# graph_dot <- tetradrunner.tetradGraphToDot(Data[[1]][["pag"]][["graph"]])
# dot(graph_dot)

  set1 <- paste0(n_iter,"-precision.csv")
  set2 <- paste0(n_iter,"-recall.csv")
  set3 <- paste0(n_iter,"-BSF.csv")
  set5 <- paste0(n_iter,"-time.csv")
  write.csv(iterFGES$precision,paste0(path,set1))
  write.csv(iterFGES$recall,paste0(path,set2))
  write.csv(iterFGES$BSF,paste0(path,set3))
  write.csv(iterFGES$time,paste0(path,set5))
  
}

#sink()
   #write.csv(true_mag_amat,"mag.csv")