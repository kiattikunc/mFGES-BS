precisionrecall <- function(pag,gtpag){
  nVars<-nrow(gtpag)
 tps <- ((pag-gtpag) ==0) & (t(pag)-t(gtpag) ==0) & (pag!=0)
 temp2 <- gtpag!=0
 temp3 <- pag!=0
 fps <-  (pag!=0)  & ((pag-gtpag) != (t(pag)-t(gtpag)))
 fns <-  (gtpag!=0) & (((pag-gtpag) !=0) | (t(pag)-t(gtpag) !=0))
 tns <-  (gtpag==0) & (pag==0)

 nCorrectEdges <- length(tps[tps==TRUE])
 nEdgesGtPag <- length(temp2[temp2==TRUE])
 nEdgesPag <- length(temp3[temp3==TRUE])

 tp <- length(tps[tps==TRUE])
 fp <- length(fps[fps==TRUE])
 fn<- length(fns[fns==TRUE])
 tn<- length(tns[tns==TRUE]) -length(diag(tns))
a <- 0.5*nEdgesGtPag
print(nVars*(nVars-1)/2)
print(a)
i <-(nVars*(nVars-1)/2)-a
print(i)
BSF <- 0.5*(0.5*tp/a+0.5*tn/i-0.5*fp/i-0.5*fn/a)
  #nCorrectEdges = nrow(((pag-gtpag)==0)& ((pag'-gtpag')==0) & ~~pag);


  precision <- nCorrectEdges/nEdgesPag
  recall <- nCorrectEdges/nEdgesGtPag

  return_list <- list('precision'=precision,'recall'=recall,'BSF'=BSF,'tp'=0.5*tp,'fp'=0.5*fp,'fn'=0.5*fn,'tn'=0.5*tn,'nEdgesGtPag'=0.5*nEdgesGtPag,'a'=a,'i'=i,'nEdgesPag'=0.5*nEdgesPag)

}
