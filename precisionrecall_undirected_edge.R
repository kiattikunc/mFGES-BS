precisionrecall <- function(pag,gtpag){
  nVars<-nrow(gtpag)
 tps <- (((pag-gtpag) ==0) & (t(pag)-t(gtpag) ==0)) & (pag!=0)
 tp <- length(tps[tps==TRUE])
 
 tps_undirect <- ((pag ==3) & (pag-gtpag ==2)& (t(pag)-t(gtpag) !=2))
 tps_bidirect <- ((pag ==3) & (gtpag ==1) & (pag-gtpag ==2)& (t(pag)-t(gtpag) ==2))
 tp_bidirect <- length(tps_bidirect[tps_bidirect==TRUE])
 tp_undirect <- length(tps_undirect[tps_undirect==TRUE])



 
 fps <-  (pag!=0)  & ((pag-gtpag) != (t(pag)-t(gtpag)))
 fps2 <-  (pag==0 & t(pag)==1)  & ((pag-gtpag) != (t(pag)-t(gtpag)))
 
 fns <-   (gtpag==0 & t(gtpag)==1)   & (((pag-gtpag) !=0) | (t(pag)-t(gtpag) !=0)) & (pag!=3 & t(pag)!=3)
 fn <- length(fns[fns==TRUE])
 

 tns <-  gtpag==0 & t(gtpag)==0 & pag==0 & t(pag)==0
 tn<- 0.5*(length(tns[tns==TRUE]) -length(diag(tns)))
 

 temp2 <- gtpag!=0
 temp3 <- (gtpag ==1)& t(gtpag== 1)
 nEdgesGtPag <- length(temp2[temp2==TRUE])
 nBiGtEdges <- 0.5*length(temp3[temp3==TRUE])
 nEdgesGtPag <- nEdgesGtPag-nBiGtEdges
 
 
 temp2 <- learned_mag!=0
 nEdgesMag <- length(temp2[temp2==TRUE])
 
 temp3 <- (learned_mag ==1)& t(learned_mag== 1)
 nBiEdges <- length(temp3[temp3==TRUE])
 
 temp3 <- (learned_mag ==3)& t(learned_mag== 3)
 nUnEdges <- 0.5*length(temp3[temp3==TRUE])
 nEdgesMag<- nEdgesMag-nBiEdges-nUnEdges
 
 

 temp4 <- pag!=0
 temp5 <- (pag ==1)& t(pag== 1)
 temp6 <- (pag ==3)& t(pag== 3)
 nEdgesPag <- length(temp4[temp4==TRUE])
 nBiEdges <- 0.5*length(temp5[temp5==TRUE])
 nUnEdges <- 0.5*length(temp6[temp6==TRUE])
 nEdgesPag <- nEdgesPag-nBiEdges-nUnEdges
 
 tp <- length(tps[tps==TRUE])

 tp<-  tp+0.5*(tp_undirect+tp_bidirect)
 fp <- 0.5*length(fps[fps==TRUE]) + 0.5*length(fps2[fps2==TRUE])

 a <- nEdgesGtPag

i <-(nVars*(nVars-1)/2)-a

BSF <- 0.5*(tp/a+tn/i-fp/i-fn/a)
  #nCorrectEdges = nrow(((pag-gtpag)==0)& ((pag'-gtpag')==0) & ~~pag);


  precision <- tp/(tp+fp)
  recall <- tp/(tp+fn)

  return_list <- list('precision'=precision,'recall'=recall,'BSF'=BSF,'tp'=tp,'fp'=fp,'fn'=fn,'tn'=tn,'nEdgesGtPag'=nEdgesGtPag,'a'=a,'i'=i,'nEdgesPag'=nEdgesPag,'nUndirectEdges'=nUnEdges)

}
