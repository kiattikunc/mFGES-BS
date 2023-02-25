
library(DOT)
graph_dot <- tetradrunner.tetradGraphToDot(tetradrunner$graph)
library(rcausal)
source("graphutils.R")


 graph_dot <- tetradrunner.tetradGraphToDot(Data[[1]][["pag"]][["graph"]])
 
 varNames <- colnames(dataset)
x <- ugraphToTetradGraph(learned_mag, varNames)
dot(x)

varstring <- paste(varnames, collapse=" ")
graphfilename <- "impossibly_long_graph_file_name_temporary.txt"

  graphfile <- .jnew("java/io/File", graphfilename)
  newug_tetrad <- .jcall("edu/cmu/tetrad/graph/GraphUtils",
                         "Ledu/cmu/tetrad/graph/Graph;",
                         "loadGraphTxt", graphfile)
numNodes <- ncol(learned_mag)
varnames <- varNames
edgelist <- c()

  newug_tetrad <- .jcast(newug_tetrad, "edu/cmu/tetrad/graph/Graph",
                         check=TRUE)
  

  
source("exportamat2toCSV.R")

set4 <- "dot.csv"

rownames(learned_mag) <- colnames(dataset)
colnames(learned_mag) <- colnames(dataset)

exportamattoCSV(learned_mag,set4)

dot(graph_dot)


dot("digraph {A -> B;  A [shape=recangle,color=red,style=filled];}")


dot("graph { ANAPHYLAXIS -- TPR;
CATECHOL	->	ARTCO2;
EXPCO2	->	ARTCO2;
ARTCO2	->	VENTALV;
BP	->	CO;
BP	->	TPR;
HR	->	CATECHOL;
CATECHOL	->	SAO2;
CATECHOL	->	TPR;
CO	->	HR;
CO	->	STROKEVOLUME;
CVP	->	LVEDVOLUME;
DISCONNECT	--	VENTMACH;
VENTTUBE	->	DISCONNECT;
HREKG	->	ERRCAUTER;
HRSAT	->	ERRCAUTER;
HRBP	->	ERRLOWOUTPUT;
EXPCO2	->	VENTLUNG;
PVSAT	->	FIO2;
LVEDVOLUME	->	HISTORY;
STROKEVOLUME	->	HISTORY;
HRBP	->	HR;
HREKG	->	HR;
HRSAT	->	HR;
LVEDVOLUME	->	HYPOVOLEMIA;
STROKEVOLUME	->	HYPOVOLEMIA;
MINVOL	->	INTUBATION;
PRESS	->	INTUBATION;
SAO2	->	INTUBATION;
VENTALV	->	INTUBATION;
VENTLUNG	->	INTUBATION;
PRESS	->	KINKEDTUBE;
VENTLUNG	->	KINKEDTUBE;
PCWP	->	LVEDVOLUME;
LVEDVOLUME	--	STROKEVOLUME;
MINVOL	->	VENTLUNG;
MINVOLSET	->	VENTMACH;
PAP	-	PULMEMBOLUS;
PRESS	-	VENTALV;
PRESS	->	VENTTUBE;
SAO2	->	PULMEMBOLUS;
SAO2	->	PVSAT;
PVSAT	->	VENTALV;
VENTALV	->	VENTLUNG;
VENTLUNG	->	VENTTUBE;
VENTTUBE	->	VENTMACH }")


dot("graph c{2--9; 9 -> 3;}")

