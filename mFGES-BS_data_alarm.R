#v9 to use obs dataset for searching by FGS
# bias the interventional learned output
# sub version 9.2 to integrate with prior from CPDAG from obs +INT i with INT i+1 dataset
# sub version 9.3 to integrate with skeleton from statistics
# sub version 9.4 to No Prior in structure learning FGS
# sub version 9.7 to Max P(obs) + P(INT) 9.7.1 use min
# sub version 9.8 to skeleton from pcalg
# sub version 9.20.6 added to the last iteration
#sink(file='myoutput.txt')

######SET UP##########
library(bnlearn)
library(rcausal)
library(DOT)
library(tictoc)
library(pcalg)
library(pheatmap)
library(tictoc)
library(Hmisc)
library(igraph)
library(stringr) 
source("code V 9/precisionrecall_undirected_edge_combine.R")
#hyper parameters

thred <- 0.50  #ratio of FGS
thred2 <- 0.5  # ratio of v-structure
thred3 <- 0.001  # the threshold of latent confounders
target <-1
n_average <-1
n_iter <-10
n_intervent <- 1000
stat_skeleton = TRUE


#load data sets
MAG_set <-'Input/ALARM_MAGtrue_L5_cML_cIL_cSL_cMISL.csv'
DAG_set <-'Input/DAGtrue_ALARM.csv'
training_set <-'Input/trainingData_ALARM_L5_1k.csv'
param_set <-'Input/trainingData_ALARM_N_10k.csv'
latent1 <- c("LVFAILURE", "SHUNT")
iterFGES <- list()
setname <- "/ALARM/"
nexp_input <-'C:/Users/kiattikun/Documents/GitHub/COmbINE-master/Data/Alarm Target=1 N=1000/'




######PRE RUN FUNCTION##########

#load observational dataset
dataset<- read.csv(training_set,header = TRUE,na.strings=c(""),check.names = FALSE)

for(i in 1:ncol(dataset)){
  
  dataset[,i] <- as.factor(dataset[,i])
  
}


dataset <- dataset[,order(colnames(dataset))]

#resample <- resample(dataset, R = 1000)


#load N dataset for CPT under true DAG.
truedag_dataset <- read.csv(param_set,header = TRUE,na.strings=c(""),check.names = FALSE)

for(i in 1:ncol(truedag_dataset)){
  
  truedag_dataset[,i] <- as.factor(truedag_dataset[,i])
  
}

truedag_dataset <- truedag_dataset[,order(colnames(truedag_dataset))]



n <- ncol(dataset)
n_obs  <- ncol(truedag_dataset)

#construct true DAG#####
DAG_RAW <- read.csv(DAG_set,header = TRUE,na.strings=c(""),check.name=FALSE)
data_matrix <- as.matrix(DAG_RAW)
data_matrix <- data_matrix[,-1]
data_matrix <- data_matrix[,-2]

# Convert friends matrix to an igraph object
g <- graph.edgelist(data_matrix, directed = TRUE)
NEL <- igraph.to.graphNEL(g)
true_graph <- as.bn(NEL)
####################



#construct true MAG#####
MAG_RAW <- read.csv(MAG_set,header = TRUE,na.strings=c(""),check.name=FALSE)

data_matrix <- as.matrix(MAG_RAW)
data_matrix_bidirect <-as.matrix(MAG_RAW)
data_matrix <- data_matrix[,-1]
data_matrix <- data_matrix[,-2]
#add bi-directed edge to true MAG (if possible)
for (row_index in 1:nrow(data_matrix_bidirect))
{
  if (data_matrix_bidirect [row_index,3]=="<->")
  {
    bidirect <- c(data_matrix_bidirect [row_index,2],data_matrix_bidirect [row_index,4])
    data_matrix <- rbind(data_matrix , bidirect)
    bidirect <- c(data_matrix_bidirect [row_index,4],data_matrix_bidirect [row_index,2])
    data_matrix <- rbind(data_matrix , bidirect)
  }
  
}



##################################


# Convert friends matrix to an igraph object
g <- graph.edgelist(data_matrix, directed = TRUE)
NEL <- igraph.to.graphNEL(g)
true_mag <- as.bn(NEL)
true_mag_amat <- amat(true_mag) 
true_mag_amat<- true_mag_amat[order(as.character(rownames(true_mag_amat))), order(as.character(colnames(true_mag_amat)))]
####################




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
      adjmat[j, i] <- 0
    }
    else if (edgevec[2] == "o-o") {
      edge <- c(edgevec[1], edgevec[3])
      adjmat[i, j] <- 3
      adjmat[j, i] <- 3
    }
    else if (edgevec[2] == "o->") {
      adjmat[i, j] <- 1
      adjmat[j, i] <- 3
    }
    else if (edgevec[2] == "<->") {
      edge <- c(edgevec[1], edgevec[3])
      adjmat[i, j] <- 1
      adjmat[j, i] <- 1
    }
    
  }
  return(adjmat)
}


#type == 'single'
# In each experiment there is n (target) intervention on each variable,
# except the last experiment which has no interventions.

#type == 'single'
# In each experiment there is n (target) intervention on each variable,
# except the last experiment which has no interventions.
for(average in 1:n_average)
{
  set.seed(5*average)
  print(paste0("######  iteration : ",average))
  for(iteration in 1:n_iter)
  {
    
    nexp <-iteration
    tic.clearlog()
    Data<-list()
    total <- n_iter+1
    round <- list()
    nexp_path <- paste0(nexp_input,nexp,average,"-nExp.csv")
    obd_nexp <- matrix(0,1,n_obs)
    m <- read.csv(nexp_path,header = TRUE,na.strings=c(""),check.names = FALSE)
    x<- nrow(m)
    m <- m[-1,-1] # remove 1st column and 1st row
    rest<- nexp+1
    for ( i_index in rest:total) {
      #read targets
      
      m <- rbind(m,obd_nexp)
      
    }
    colnames(m) <- colnames(truedag_dataset)
    m <- m[ , !(names(m) %in% latent1)]
    E <- m 
    
    
    
    
    #####################################################################
    
    #generate n interventional dataset
    
    for ( i_index in 1:total) {
      
      # Create the tuple that will be stored in D, storing the vector e.
      Data[[i_index]]<-list(e=E[i_index,])
      # Get indexes of intervened variables in this particular experimental setting.
      J<-which( Data[[i_index]]$e==1 )
      # The data consist of manipulated graphs where
      # Edge heads into the intervened variables are cut.
      
      Data[[i_index]]$data <- truedag_dataset
      interventiontarget_dataset <- dataset[,order(colnames(dataset))]
      #skeleton = hc(Data[[i_index]]$data,debug=FALSE)
      skeleton <- true_graph
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
        
        for ( incoming_index in 1: nrow(list_to_remove)) {
          skeleton <- drop.arc(skeleton, list_to_remove[incoming_index,1], list_to_remove[incoming_index,2])
        }
        
      }
      
      
      if (length(J) !=0 & nrow(list_to_remove) !=0 )
      {
        
        param_skeleton <- bn.fit(skeleton, Data[[i_index]]$data, method = "bayes", keep.fitted = TRUE, debug = TRUE)
        
        
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
        interventional_dataset <- interventional_dataset[,order(colnames(interventional_dataset))]
        
        
        
        #export int-data
        # interventional_dataset2 <-interventional_dataset
        # cols_int_all <- colnames(interventional_dataset2)
        # for(cols_int in cols_int_all){
        #   interventional_dataset2[[cols_int]]<- as.numeric(interventional_dataset2[[cols_int]])
        #   interventional_dataset2[[cols_int]]<- interventional_dataset2[[cols_int]]-1
        # }
        # set_int <- paste0('Data/average-',average,"-",iteration,"-",i,"-alarm_N_10k.csv")
        # 
        # write.table(interventional_dataset2,paste0(set_int),row.names = FALSE, col.names=FALSE,sep=",")
        # 
        
        #remove latent variables
        interventional_dataset <- interventional_dataset[ , !(names( interventional_dataset) %in% latent1)]
        
        Data[[i_index]]$data <- interventional_dataset
        Data[[i_index]]$param <- param_skeleton
        
        
      }
      
      else if (length(J) !=0 & nrow(list_to_remove) ==0 )
      {
        param_skeleton <- bn.fit(skeleton, Data[[i_index]]$data, method = "bayes", keep.fitted = TRUE, debug = TRUE)
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
        interventional_dataset <- interventional_dataset[,order(colnames(interventional_dataset))]
        #remove latent variables
        
        interventional_dataset <- interventional_dataset[ , !(names( interventional_dataset) %in% latent1)]
        
        
        Data[[i_index]]$data <- interventional_dataset
        Data[[i_index]]$param <- param_skeleton
        
        
      }
      
      else if (length(J) ==0)
      {
        interventional_dataset <- dataset[ , !(names( dataset) %in% latent1)]
        Data[[i_index]]$data <- interventional_dataset
        
      }
      
      
    }
    
    
    
    
    
   
    
    ############## START ############################################
    tic("Tetradtime fges")
    
    ###################### 1 run FGS on obs data #####################
    
    tetradrunner <- tetradrunner(algoId = 'fges',
                                 df = Data[[total]]$data,
                                 scoreId = 'bdeu',
                                 dataType = 'discrete',
                                 faithfulnessAssumed = TRUE,
                                 maxDegree = -1, #argv$maxDegree,
                                 verbose = FALSE)
    
    
    Data[[total]]$pag <-tetradrunner
    
    #################################################################
    #### 1.1 BOX1 run frequrentist method from FGS
    
    maj_amat <- matrix(0, n, n)    # majority Rule
    rownames(maj_amat) <- colnames(dataset)
    colnames(maj_amat) <- colnames(dataset)
    
    int_rule <- matrix(0, n, n)   # intervnetion Rule
    rownames(int_rule) <- colnames(dataset)
    colnames(int_rule) <- colnames(dataset)
    
    
    amat<- matrix(0, n, n)
    amat_undirect <- matrix(0, n, n)
    rownames(amat) <- colnames(dataset)
    colnames(amat) <- colnames(dataset)
    
    edges <- Data[[total]]$pag$edges
    
    amat <- TetradGetAdjmat(Data[[total]]$pag)
    rownames(amat) <- colnames(dataset)
    colnames(amat) <- colnames(dataset)
    
    # for undirected edge for intv dataset
    amat_obs <- amat
    ##split undirect from matrix
    amat_undirect [ amat ==t(amat) &  amat ==3 ] <- 3
    amat[ amat ==t(amat) &  amat ==3 ] <- 0
    Data[[total]]$amat <-amat
    
    ################################################################
    #### 1.2  BOX2: use the skeleton from essential graph from BOX1 and perform the statistic test to identify v-structure using conservative rule with majority rule
    
    amat_stat <- matrix(0, n, n)
    

    obs_index <- total # use observational dataset to identify skeleton and v-structure
    varNames <- colnames(Data[[obs_index]]$data)
    
    
    df<-Data[[obs_index]]$data
    
    
    dataset_pcalg <-Data[[obs_index]]$data
    cols <- colnames(dataset_pcalg)
    for(j in cols){
      dataset_pcalg[[j]]<- as.integer(dataset_pcalg[[j]])
      dataset_pcalg[[j]]<- dataset_pcalg[[j]]-1
    }
    datamatrix <- data.matrix(na.omit(dataset_pcalg)) # convert rawdata from category to level 0, 1,2,..
    
    suffStat = list(dm = dataset_pcalg, nlev = sapply(Data[[obs_index]]$data , nlevels),adaptDF = FALSE)
   
    amat_skeleton <- matrix(0, n, n)
    amat_skeleton [ amat_obs != t(amat_obs)] <- 1
    amat_skeleton [ amat_obs == t(amat_obs) & amat_obs!=0] <- 1
    
    if(stat_skeleton)
    {
      source("code V 9/checkDependencies.R")
      print("calcualting skeleton")
      tetradrunner_fci <- tetradrunner(algoId = 'fci',
                                       df = df,
                                       dataType = 'discrete',
                                       faithfulnessAssumed = TRUE,
                                       maxDegree = -1, #argv$maxDegree,
                                       verbose = FALSE)
      
      edges <- tetradrunner_fci$edges
      amat_stat <- TetradGetAdjmat(tetradrunner_fci)
      # # 
      # skel.fit <- pcalg::skeleton(suffStat , n = n,
      #                             indepTest = disCItest, ## (partial correlations)
      #                             alpha = 0.5, labels = c(varNames), verbose = FALSE)
      # amat_stat <- as(skel.fit, "amat")
      
      amat_skeleton <- matrix(0, n, n)
      amat_skeleton [ amat_obs != t(amat_obs)] <- 1
      amat_skeleton [ amat_obs == t(amat_obs) & amat_obs!=0] <- 1
      amat_skeleton [ amat_stat == t(amat_stat) & amat_stat!=0] <- 1
      
    }
    
    
   
    
    ##################################################
    
    #### 1.3  BOX3  run the local Bdeu on obs data
    nodesx <-tetradrunner$graph
    nodes <- nodesx$getNodes()
    
    local_BIC <- matrix(0, nrow = n+1, ncol = nexp+2)
    local_BIC[n+1,1]<- 'model score'  # create local Bdeu table 

    for(i_bic in 0:as.integer(nodes$size()-1))
    {
      local_BIC[i_bic+1,1] <-nodes$get(as.integer(i_bic))$getName()
      local_BIC[i_bic+1,nexp+2] <-nodes$get(as.integer(i_bic))$getAttribute('BIC')
      local_BIC[n+1,nexp+2]<- nodesx$getAttribute('BIC')
      
    }
    
    ###################### 2 run FGS on INT data #####################
   
    # start to run the experiment on INT data 
    if (nexp >0)
    {
      
      for ( n_row in 1:(nexp) ) # run INT-1 experiment and use the prior to cacualte posetia for the last INT
      {
        
        
        print(paste0("Perform the eperiment:",n_row))
        
        #pre_amat <- Data[[total]]$amat
        #pre_amat [ pre_amat ==1  ] <- 0.5
       # pre_amat [ pre_amat ==t(pre_amat) & pre_amat ==3  ] <- 0.5
        
    
       
        tetradrunner <- tetradrunner(algoId = 'fges',
                                     df = Data[[n_row]]$data,
                                     scoreId = 'bdeu',
                                     dataType = 'discrete',
                                     faithfulnessAssumed = TRUE,
                                     maxDegree = -1, #argv$maxDegree,
                                     verbose = FALSE,
                                    
        )
        
        
        nodesx <-tetradrunner$graph
        nodes <- nodesx$getNodes()
        
        for(i_bic in 0:as.integer(nodes$size()-1))
        {
          
          local_BIC[i_bic+1,n_row+1] <-nodes$get(as.integer(i_bic))$getAttribute('BIC')
          local_BIC[n+1,n_row+1]<- nodesx$getAttribute('BIC')
          
        }
        
        
        ############
        xk <- list()
        edges <- tetradrunner$edges
        X <-  TetradGetAdjmat(tetradrunner)
        xk <-precisionrecall(X, true_mag_amat)
        path <-paste0("code V 9/TestINT/Target=",target," N=",n_intervent,setname)
        
        setxk <- paste0(average,"-",n_iter,"-",n_row,"-xk.csv")
        write.csv(xk,paste0(path,setxk))
        
        
        #####################
        
        
        Data[[n_row]]$pag <-tetradrunner
        
        #### 2.1 BOX1 run frequentest method from FGS
       
        amat<- matrix(0, n, n)
        amat_undirect <- matrix(0, n, n) # for 3 value undirected edge for intv dataset
        amat_fgs  <- matrix(0, n, n)
        
        rownames(amat) <- colnames(dataset)
        colnames(amat) <- colnames(dataset)
        
        edges <- Data[[n_row]]$pag$edges
        
        amat<- TetradGetAdjmat(Data[[n_row]]$pag)
        rownames(amat) <- colnames(dataset)
        colnames(amat) <- colnames(dataset)
        
        
        ##split undirect from matrix
        Data[[n_row]]$amat <-amat
        
        ### calculate prob of FGS
        edge_count <- matrix(0, n, n)
        rownames(edge_count) <- colnames(dataset)
        colnames(edge_count) <- colnames(dataset)
        
        for (i_nexp in 1:(n_row))  
        {

            for (i_n_row in 1:(i_nexp))
            {
              
              
              # reduce the weight for intervention target
              for ( var_index in 1: length(Data[[i_n_row]][["e"]]))
              {
                # check if the variable is interventional target
                if( Data[[i_n_row]][["e"]][var_index] == 1 )
                  
                { 
                 
                  
                  parent <- which(Data[[i_n_row]]$amat[,var_index]  >0, arr.ind = TRUE)
                  if (length(parent)==1)
                  {
                    print ('update rule')
                    Data[[i_n_row]]$amat[parent,var_index]  <- 0
                    if(Data[[i_n_row]]$amat[var_index,parent]==3)
                    {
                      Data[[i_n_row]]$amat[var_index,parent]  <- 0.5
                    }
                    if(Data[[i_n_row]]$amat[var_index,parent]==1)
                    {
                      Data[[i_n_row]]$amat[var_index,parent]  <- 1
                    }
                    
                   
                   
                    int_rule[var_index,parent]  <- 1
                    
                  }
                  else if(length(parent)>0)
                  {
                    for (parent_i in 1: length(parent))
                    {
                      Data[[i_n_row]]$amat[parent[parent_i],var_index]  <- 0
                      if(  Data[[i_n_row]]$amat[var_index,parent[parent_i]] ==1)
                      {
                        Data[[i_n_row]]$amat[var_index,parent[parent_i]]  <- 1
                      }
                      if(  Data[[i_n_row]]$amat[var_index,parent[parent_i]] ==3)
                      {
                        Data[[i_n_row]]$amat[var_index,parent[parent_i]]  <- 0.5
                      }
                      
                  
                      int_rule[var_index,parent[parent_i]]  <- 1
                    }
                    
                  }
                  
                  
                }
              }
              
            }
            edge_count[ Data[[total]]$amat ==1 | t(Data[[total]]$amat) ==1] <- 1
    
            for (i_n_row in 1:(i_nexp))
            {
              edge_count_tmp <- matrix(0, n, n)
              edge_count_tmp [ Data[[i_n_row]]$amat ==1 | t(Data[[i_n_row]]$amat) ==1] <- 1
              edge_count <- edge_count+ edge_count_tmp
              
              amat_int_run <- Data[[i_n_row]]$amat
              amat_int_run [ amat_int_run ==t(amat_int_run)  ] <- 0
              amat_fgs= amat_fgs+amat_int_run
     
            }
        
            amat_direct_obs<- Data[[total]]$amat
            amat_direct_obs [ amat_direct_obs ==t(amat_direct_obs) & amat_direct_obs ==3] <- 0
            edge_count [ edge_count < 0 ] <- 0
            amat <-(amat_fgs+amat_direct_obs)/edge_count
            amat[is.nan(amat)] <-0
            amat [ t(amat)!=0.5  &  amat ==0.5 ] <- 0  
            amat [ amat >0  &  t(int_rule)  ==1 ] <- 0  #int_rule  # Intervention rule
            amat [int_rule  ==1 ] <- 0.5
            amat[ amat >1] <- 1
            
      
            ##  skeleton from FGES##
            edges <- Data[[total]]$pag$edges
            amat_obs  <- TetradGetAdjmat(Data[[total]]$pag)
            
            amat_undirect [ amat_stat == t(amat_stat) & amat_stat==1] <- 3
            amat_undirect [ amat_stat == t(amat_stat) & amat_stat==3] <- 3
            amat_undirect [ amat_stat == 1 & t(amat_stat)==0] <- 3
            amat_undirect [ amat_stat == 1 & t(amat_stat)==3] <- 3
            amat[amat<=0.5 & t(amat)<=0.5 & amat_undirect==3 ] <- 0.5
           
            
            amat[ amat ==t(amat) & amat ==0 &  amat_undirect ==3 ] <- 0.5
            #amat[ amat_pvalue >0.05 ] <- 0 ## remove dependencies
            amat[amat<=0.5 & t(amat)<=0.5 & amat_stat != t(amat_stat) & amat_stat !=0] <- 0.5
            amat[amat<=0.5 & t(amat)<=0.5 & amat_stat == t(amat_stat) & amat_stat ==3] <- 0.5
            amat[amat<=0.5 & t(amat)<=0.5 & amat_stat == t(amat_stat) & amat_stat ==1] <- 0.5
            
            Data[[i_nexp]]$prior <- amat
           
        }   
        
        
        print('calculate prob of FGS completed')
        
        
      }

    }
    
    
    
    ####### 2.2 Box  added majority rule

    #amat_undirect_obs[amat_skeleton_stat==0 & t(amat_skeleton_stat)==0] <-0
  
    amat[amat<=0.5 & t(amat)<=0.5 & amat_stat==1 ] <- 0.5
    amat[amat<=0.5 & t(amat)<=0.5 & amat_stat != t(amat_stat) & amat_stat !=0] <- 0.5
    amat[amat>0  & amat_skeleton==0 ] <- 0  #spurouse edge

    amat_undirect [ amat_stat == t(amat_stat) & amat_stat==1] <- 3
    amat_undirect [ amat_stat == t(amat_stat) & amat_stat==3] <- 3
    amat_undirect [ amat_stat == 1 & t(amat_stat)==0] <- 3
    amat_undirect [ amat_stat == 1 & t(amat_stat)==3] <- 3
    amat[amat<=0.5 & t(amat)<=0.5 & amat_undirect==3 ] <- 0.5

    
    skeleton <- new("pcAlgo")
    skeleton1 <- new("graphAM",amat_skeleton)
    skeleton2 <- as(skeleton1, "graphNEL") # convert amat format to pcalg graph format
    skeleton@graph <-skeleton2
    
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
    source("code V 9/pc.cons.intern2.R")
    
    print('calculate v-structure')
    
    rule0_graph <- pc.cons.intern2(skeleton, suffStat, indepTest=disCItest, alpha=0.01, version.unf = c(1, 1),maj.rule = TRUE, verbose = FALSE)
    
    #pag <- fci(suffStat, indepTest=disCItest,rules = rep(TRUE,10),alpha = 0.01,label=c(varNames))
    
    unsh <- find.unsh.triple(amat_skeleton, check=TRUE)
    # unsh <- find.unsh.triple(amat, check=TRUE)
    vstruct <- unsh$unshTripl
    
    
    ###################soft contrain v-structue integrated to prior
  
    amat[ amat ==t(amat) & amat ==0 &  amat_undirect ==3 ] <- 0.5
    #amat[ amat_pvalue >0.05 ] <- 0 ## remove dependencies
    amat[amat<=0.5 & t(amat)<=0.5 & amat_stat != t(amat_stat) & amat_stat !=0] <- 0.5
    amat[amat<=0.5 & t(amat)<=0.5 & amat_stat == t(amat_stat) & amat_stat ==3] <- 0.5
    amat[amat<=0.5 & t(amat)<=0.5 & amat_stat == t(amat_stat) & amat_stat ==1] <- 0.5
    
    for(index_sepset_a in 1:n)
    {
      for(index_sepset_c in 1:n)
      {
        for(index_sepset_b in 1:n)
        {
          maj = rule0_graph[["sepsetratio"]][[index_sepset_a]][[index_sepset_c]][[index_sepset_b]]

          if ( maj >= 0 & maj < 0.5 & !is.na(maj))
          {
            a<- which(vstruct[1,] == index_sepset_a & vstruct[3,] == index_sepset_c & vstruct[2,] == index_sepset_b , arr.ind = TRUE)
            b<- which(vstruct[1,] == index_sepset_c & vstruct[3,] == index_sepset_a  & vstruct[2,] == index_sepset_b , arr.ind = TRUE)

            if (length(a) ==1 || length(b) ==1)

            {
              collider<- index_sepset_b
              # ratio <-1-rule0_graph[["sepsetratio"]][[index_sepset_i]][[index_sepset_j]]
              ratio <-max(1-maj,0.5)
              if (ratio>thred2)
              {
                print(paste0(colnames(amat)[index_sepset_a],'-',colnames(amat)[collider],'-',colnames(amat)[index_sepset_c], ' ratio=',ratio))

                amat[index_sepset_a,collider] <- max(amat[index_sepset_a,collider],ratio)
                amat[index_sepset_c,collider] <- max(amat[index_sepset_c,collider],ratio)


              #  maj_amat[index_sepset_a,collider] <- max(amat[index_sepset_a,collider],ratio)
               # maj_amat[index_sepset_c,collider] <- max(amat[index_sepset_c,collider],ratio)

                
                maj_amat[index_sepset_a,collider] <-max(maj_amat[index_sepset_a,collider],ratio)
                maj_amat[index_sepset_c,collider] <-max(maj_amat[index_sepset_c,collider],ratio)

              }

            }

          }

        }



      }

    }


    amat_int<- amat
    prior_int_0 <-amat
    prior_int_0 <- pmax(prior_int_0,maj_amat)
  
    for (i_nexp in 1:(nexp))  
    {
      Data[[i_nexp]]$prior<- pmax(Data[[i_nexp]]$prior,maj_amat)
     
    }
   
   # amat_int<- pmax(amat,maj_amat)
    
    
    ####### 2.3  BOX3  run the local Bdeu on obs data
   
    learned_latent_conf <- matrix(0, n, n)
    amat_interv_eff <- matrix(0, n, n)
   
    if (nexp >0) 
    {
      for ( n_row in 1:(nexp) ) # run INT-1 experiment and use the prior to cacualte posetia for the last INT
      {
        
        ###################### 2.3 BOX 3

        
        ratio_BIC <- local_BIC
        tempt <- ncol(ratio_BIC)
        ratio_BIC[,2:tempt] <-  abs(1-(as.double(ratio_BIC[,2:tempt]) /  as.double(ratio_BIC[, ncol(ratio_BIC)])))
        
        # 
        for ( nexp_index in 1:(tempt-1))
        {
          
          for ( var_index in 1: length(Data[[nexp_index]][["e"]]))
          {
            
            
            if(as.double(ratio_BIC[var_index,nexp_index+1]) >1 )
            {
              ratio_BIC[var_index,nexp_index+1] <- 1
            }
            
          }
          
        }
        
        ############# 2.3.1 check the local Bdeu 

        for ( var_index in 1: length(Data[[n_row]][["e"]]))
        {
          
          if( Data[[n_row]][["e"]][var_index] == 1 )
            
          { 
            
            print(paste0("Experiment:", n_row ," has ",colnames(amat_int)[var_index]))
            for (row_index in 1:length(Data[[n_row]][["e"]]) )
            {
              
              #list the variable which is not interventional target
              if (row_index !=var_index)
              {
                
                
                if (as.double(ratio_BIC[row_index,n_row+1]) >thred3 & Data[[n_row]]$amat[var_index,row_index]>0)
                  # else if (amat_undirect[var_index,row_index]==1)
                  
                {
                  amat_interv_eff[var_index,row_index]<-min(amat_interv_eff[var_index,row_index]+as.double(ratio_BIC[row_index,n_row+1]),1)
                  amat_int[var_index,row_index] <- min(amat_int[var_index,row_index]+as.double(ratio_BIC[row_index,n_row+1]),1)
                  print(paste0(" enhance ",colnames(amat_int)[var_index]," -> ",  colnames(amat_int)[row_index],' by ', as.double(ratio_BIC[row_index,n_row+1])))
                  
                }
                else if (as.double(ratio_BIC[row_index,n_row+1]) <thred3 & amat[var_index,row_index]>0)
                  # else if (amat_undirect[var_index,row_index]==1)
                  
                {
                  learned_latent_conf[var_index,row_index]<-max(learned_latent_conf[var_index,row_index]+as.double(ratio_BIC[row_index,n_row+1]),0)
                  print(paste0(" found the possiblity of latent confounders ",colnames(amat_int)[var_index]," <-> ",  colnames(amat_int)[row_index],' with the score ', as.double(ratio_BIC[row_index,n_row+1])))
                  
                }
                
                
              }
            }
            
          }
        }
        
        Data[[n_row]]$prior <- Data[[n_row]]$prior+amat_interv_eff
        Data[[n_row]]$prior <- pmin(Data[[n_row]]$prior,1) 
       
        }
      
      
    
      
    }
   
   
     
    
 
    
    ###################### catch Bdeu score for constraints
    source("code V 9/posteriorcal3.R")
    
    print(paste0('calcualting the posterior data:',nexp))
    indep_graphs <- empty.graph(varNames)
    posterior<- posteriorcal(prior_int_0,indep_graphs,Data,n,nexp)

    posterior_int <-posterior$indepence_posterior_INT
    posterior_int_weak <-posterior$indepence_posterior_INT
   
    
    rownames(posterior_int) <- colnames(dataset)
    colnames(posterior_int) <- colnames(dataset)
    
    
    posterior_int<- 1-posterior_int
    diag(posterior_int) <- 0
    
    #pick the higest posteria prob when less than threshold
    #posterior_obs[posterior_obs< thred & t(posterior_obs) <thred & amat_undirect==3 & (posterior_obs> t(posterior_obs)) ] <- 1

    learned_mag <- matrix(0, n, n)
    rownames(learned_mag) <- colnames(dataset)
    colnames(learned_mag) <- colnames(dataset)

    rownames(amat_undirect) <- colnames(dataset)
    colnames(amat_undirect) <- colnames(dataset)
    
    learned_mag[ posterior_int > thred ] <- 1
    learned_mag [ posterior_int < thred ] <- 0
    learned_mag [ posterior_int > thred  & posterior_int > t(posterior_int)] <- 1
    learned_mag [ posterior_int > thred  & posterior_int < t(posterior_int)] <- 0
    learned_mag[posterior_int==0 & t(posterior_int)==0 & amat_undirect==3 ] <- 3
    learned_mag [ posterior_int > thred  & round(posterior_int, digits = 4) ==  round(t(posterior_int), digits = 4)] <- 3
   
   
    learned_latent_conf[learned_mag==t(learned_mag) & learned_mag==1]
    
    
    ####### check_cyclic #######
    
    
    checked_pag <- matrix(0, n, n)
    rownames(checked_pag) <- c(1:n)
    colnames(checked_pag) <- c(1:n)
    checked_pag[learned_mag==1 & t(learned_mag)==0 ] <- 1
    checked_pag[learned_mag==1 & t(learned_mag)==1 ] <- 0 # remove bi-directed edges
    checked_pag[learned_mag==3 & t(learned_mag)==3 ] <- 0 # remove undirected edges
    check_cyclic <- pcalg::isValidGraph(checked_pag, type = "pdag", verbose = FALSE)
    
    
    # g <- graph.adjacency(checked_pag)
    # plot(g)
    
    source("code V 9/remove_weak_cycle.R")
    while(!check_cyclic)
    {
      
      removed_graph <- remove_weak_cycle(checked_pag,learned_mag,posterior_int_weak,n)
      learned_mag <- removed_graph[["remove_graph"]] 
      
      
      checked_pag <- matrix(0, n, n)
      checked_pag[learned_mag==1 & t(learned_mag)==0 ] <- 1
      checked_pag[learned_mag==1 & t(learned_mag)==1 ] <- 0 # remove bi-directed edges
      checked_pag[learned_mag==3 & t(learned_mag)==3 ] <- 0 # remove undirected edges
      check_cyclic <- pcalg::isValidGraph(checked_pag, type = "pdag", verbose = TRUE)
      print(paste0("cyclic=",check_cyclic))
      
    }
    
    
    
    
    
    ################################################################
    
    
    
    
    
    
    
    
    
    
    ############## STOP THE mFGS-BS #########################
    toc(log = TRUE, quiet = TRUE)
    log.txt <- tic.log(format = TRUE)
    #############################################
    
    
    
    
    ### Evaluation section  #################################
   
    true_mag_amat <- amat(true_mag)
    
    true_mag_amat<- true_mag_amat[order(rownames(true_mag_amat)),order(colnames(true_mag_amat))]
    
    
    k <- list()
# edges <- Data[[11]]$pag$edges
   #rownames(X) <- rownames(true_mag_amat)
#colnames(X) <- colnames(true_mag_amat)

  # k <-precisionrecall(amat, true_mag_amat)
    k <-precisionrecall(learned_mag, true_mag_amat)
    
    ### END Evaluation section #########
    
    
    
    
    path <-paste0("code V 9/TestINT/Target=",target," N=",n_intervent,setname)
    iterFGES$nExp[[iteration]] <- nexp
    iterFGES$precision[[iteration]] <- k[["precision"]] 
    iterFGES$recall[[iteration]] <- k[["recall"]]
    
    iterFGES$BSF[[iteration]] <- k[["BSF"]]
    iterFGES$time[[iteration]] <- log.txt[[1]]
    #iterFGES$cyclic[[iteration]] <-check_cyclic
    iterFGES$nEdge[[iteration]] <-k[["nEdgesPag"]]
    
    
    
    set4 <- paste0(average,"-nExp.csv")
    set5 <- paste0(average,"-learn_mag.csv")
    set6 <- paste0(average,"-learn_latent_conf.csv")
    
    
    write.csv(learned_mag,paste0(path,iteration,set5))
    
    write.csv(learned_latent_conf,paste0(path,iteration,set6))
    
   
    ##################### END EXPORT ####################### 
    
    
    
    source("graphutils.R")
    #  library(DOT)
    #
    # graph_dot <- tetradrunner.tetradGraphToDot(Data[[1]][["pag"]][["graph"]])
    # dot(graph_dot)
    # graph_dot <- tetradrunner.tetradGraphToDot(Data[[3]][["pag"]][["graph"]])
    # 
    #  dot(graph_dot)
    #varNames <- colnames(dataset)
    #x <- ugraphToTetradGraph(learned_mag_int1, varNames)
    # dot(x)
    
    set1 <- paste0(average,"-",n_iter,"-precision.csv")
    set2 <- paste0(average,"-",n_iter,"-recall.csv")
    set3 <- paste0(average,"-",n_iter,"-BSF.csv")
    #set4 <- paste0(average,"-",n_iter,"-cyclic.csv")
    set5 <- paste0(average,"-",n_iter,"-time.csv")
    set6 <- paste0(average,"-",n_iter,"-nedge.csv")
    write.csv(iterFGES$precision,paste0(path,set1))
    write.csv(iterFGES$recall,paste0(path,set2))
    write.csv(iterFGES$BSF,paste0(path,set3))
    #write.csv(iterFGES$cyclic,paste0(path,set4))
    write.csv(iterFGES$time,paste0(path,set5))
    write.csv(iterFGES$nEdge,paste0(path,set6))
    
  }
}



  source("graphutils.R")
 # library(DOT)
#graph_dot <- tetradrunner.tetradGraphToDot(Data[[3]][["pag"]][["graph"]])
# dot(graph_dot)


#  library(DOT)
# dot("digraph g {\n \"asia\" -> \"tub\" [arrowtail=none, arrowhead=normal]; 
# \"asia\" -> \"either\" [arrowtail=none, arrowhead=normal];
# \"either\" -> \"xray\" [arrowtail=none, arrowhead=normal];  
# \"lung\" -> \"either\" [arrowtail=none, arrowhead=normal];
# \"smoke\" -> \"lung\" [arrowtail=none, arrowhead=normal]; 
# \"smoke\" -> \"bronc\" [arrowtail=none, arrowhead=normal];
# \"bronc\" -> \"dysp\" [arrowtail=none, arrowhead=normal];
# \"either\" -> \"dysp\" [arrowtail=none, arrowhead=normal];dysp [style=filled, color=lightgrey];}")
# #either [shape=square]; 
# graph_dot <- tetradrunner.tetradGraphToDot(Data[[10]][["pag"]][["graph"]])
# dot(graph_dot)

#dot("digraph g {\n \"asia\" -> \"tub\" [arrowtail=none, arrowhead=none]; asia [shape=square];\n}")


#graph_dot <- tetradrunner.tetradGraphToDot(tetradrunner_fci$graph)
#dot(graph_dot)
#learned_mag <- read.csv(paste0(path,'41-learn_mag.csv'),header = FALSE,na.strings=c(""),check.name=FALSE)
 # colnames(dataset) <- str_replace_all(colnames(dataset), "[^[:alnum:]]", "")
#   varNames <- colnames(dataset)
#x <- ugraphToTetradGraph(amat(indep_graphs), varNames)
#     dot(x)
