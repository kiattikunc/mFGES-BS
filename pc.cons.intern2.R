##########################################################
pc.cons.intern2 <- function(sk, suffStat, indepTest, alpha,
                           version.unf = c(NA,NA), maj.rule = FALSE, verbose = FALSE)
{
  ## Purpose:  For any unshielded triple A-B-C, consider all subsets D of
  ## the neighbors of A and of the neighbors of C, and record the sets
  ## D for which A and C are conditionally independent given D. If B
  ## is in none of these sets, do nothing (it is a
  ## v-structure) and also delete B from sepset(A,C) if present (so we are
  ## sure that a v-structure will be created). If B is in all sets, do nothing
  ## (it is not a v-structure) and also add B to sepset(A,C) if not present
  ## (so we are sure that a v-structure will not be created). If maj.rule=FALSE
  ## the normal conservative version is applied, hence if B is in
  ## some but not all sets, mark the triple as "ambiguous". If maj.rule=TRUE
  ## we mark the triple as "ambiguous" if B is in exactly 50% of the cases,
  ## if less than 50% define it as a v-structure, and if in more than 50%
  ## no v-structure.
  ## ----------------------------------------------------------------------
  ## Arguments: - sk: output returned by function "skeleton"
  ##            - suffStat: Sufficient statistics for independent tests
  ##            - indepTest: Function for independence test
  ##            - alpha: Significance level of test
  ##            - version.unf[1]: 1 it checks if b is in some sepsets,
  ##                              2 it also checks if there exists a sepset
  ##                              which is a subset of the neighbours.
  ##            - version.unf[2]: 1 same as in Tetrad (do not consider
  ##                              the initial sepset), 2 it also considers
  ##                              the initial sepset
  ##            - maj.rule: FALSE/TRUE if the majority rule idea is applied
  ## ----------------------------------------------------------------------
  ## Value: - unfTripl: Triple that were marked as unfaithful
  ##        - vers: vector containing the version (1 or 2) of the
  ##                corresponding triple saved in unfTripl (1=normal
  ##                unfaithful triple that is B is in some sepsets;
  ##                2=triple coming from version.unf[1]==2
  ##                that is a and c are indep given the initial sepset
  ##                but there doesn't exist a subset of the neighbours
  ##                that d-separates them)
  ##        - sk: updated skelet object, sepsets might have been updated
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 12 Feb 2010, 10:43
  ## Modifications: Diego Colombo
  
  g <- as(sk@graph,"matrix")
  stopifnot(all(g == t(g))) ## g is guaranteed to be symmetric
  p <- as.numeric(dim(g)[1]) ## p is number of node
  
  unfTripl <- vers <- rep(NA,min(p*p,100000)) ## index p x p
  counter <- 0
  
  
  sepsetratio <- sk@sepset
  if (sum(g) > 0) {
    ind <- which(g == 1, arr.ind = TRUE) # logical; should array indices be returned when x is an array
    tripleMatrix <- NULL
    ## Go through all edges
    for (i in seq_len(nrow(ind))) {
      a <- ind[i,1]  # which one is a
      b <- ind[i,2] # which one is b
      allC <- setdiff(which(g[b,] == 1),a) ## a-b-c
      newC <- allC[g[a,allC] == 0]
      tmpMatrix <- cbind(rep(a,length(newC)),rep(b,length(newC)),newC)
      tripleMatrix <- rbind(tripleMatrix,tmpMatrix)
      colnames(tripleMatrix) <- c("","","")
    }
    if ((m <- nrow(tripleMatrix)) > 0) {
      deleteDupl <- logical(m)# all FALSE
      for (i in seq_len(m))
        if (tripleMatrix[i,1] > tripleMatrix[i,3])
          deleteDupl[i] <- TRUE
        if(any(deleteDupl))
          tripleMatrix <- tripleMatrix[!deleteDupl,, drop = FALSE]
        
        for (i in seq_len(nrow(tripleMatrix))) {
          ## pay attention to the size of counter
          if (counter+1L == length(unfTripl)) {
            n.xtra <- min(p*p, 100000)
            new.len <- counter+1L + n.xtra
            length(unfTripl) <- new.len
            length(vers)     <- new.len
          }
          a <- tripleMatrix[i,1]
          b <- tripleMatrix[i,2]
          c <- tripleMatrix[i,3]
          nbrsA <- which(g[,a] != 0) ## G symm; c no nbr of a
          nbrsC <- which(g[,c] != 0)
          if (verbose) {
            cat("\nTriple:", a,b,c,"and sepset by skelet:",
                unique(sk@sepset[[a]][[c]],sk@sepset[[c]][[a]]),"\n")
          }
          r.abc <- checkTriple(a, b, c, nbrsA, nbrsC,
                               sk@sepset[[a]][[c]], sk@sepset[[c]][[a]],
                               suffStat = suffStat, indepTest = indepTest, alpha = alpha,
                               version.unf = version.unf, maj.rule = maj.rule, verbose = verbose)
          
          sepsetratio[[a]][[c]]<- r.abc$ratio
          sepsetratio[[c]][[a]]<- r.abc$ratio
         
          
          if (verbose) {
            cat("\ndecision:", r.abc$decision,"\n")
          }
           # r.abc  = results of a-b-c
          ## 1: in NO set; 2: in ALL sets; 3: in SOME but not all
          ## Take action only if case "3"
          if (r.abc$decision == 3) {
            ## record ambiguous triple
            counter <- counter + 1
            unfTripl[counter] <- triple2numb(p,a,b,c)
            vers[counter] <- r.abc$version
          }
          ## can happen the case in Tetrad, so we must save the triple
          ## as ambiguous:
          ## a and c independent given S but not given subsets of the
          ## adj(a) or adj(c)
          if ((version.unf[1] == 2) && (r.abc$version == 2) && (r.abc$decision != 3)) {
            counter <- counter + 1
            unfTripl[counter] <- triple2numb(p,a,b,c)
            vers[counter] <- r.abc$version
          }
          sk@sepset[[a]][[c]] <- r.abc$SepsetA
       
          sk@sepset[[c]][[a]] <- r.abc$SepsetC
          
          
         
          #print(r.abc$ratio)
          
        }
    }
  }
  length(unfTripl) <- length(vers) <- counter
  list(unfTripl = unfTripl, vers = vers, sk = sk, sepsetratio=sepsetratio)
}


## Called both from pc.cons.intern() and rfci.vStruc() :
checkTriple <- function(a, b, c, nbrsA, nbrsC, sepsetA, sepsetC,
                        suffStat, indepTest, alpha, version.unf = c(NA,NA),
                        maj.rule = FALSE, verbose = FALSE)
{
  ## Purpose: For each subset of nbrsA and nbrsC where a and c are cond.
  ## independent, it is checked if b is in the conditioning set.
  ## ----------------------------------------------------------------------
  ## Arguments: - a,b,c: Nodes (positions in adjacency matrix)
  ##            - nbrsA: Neighbors of a
  ##            - nbrsC: Neighbors of c
  ##            - sepsetA: sepset(a,c)
  ##            - sepsetC: sepset(c,a)
  ##            - suffStat: Sufficient statistics for independent tests
  ##            - indepTest: Function for independence test
  ##            - alpha: Significance level of test
  ##            - version.unf[1]: 1 it checks if b is in some sepsets,
  ##                              2 it also checks if there exists a sepset
  ##                              which is a subset of the neighbours.
  ##            - version.unf[2]: 1 same as Tetrad (do not consider the initial
  ##                              sepset), 2 consider if b is in sepsetA
  ##                              or sepsetC
  ##            - maj.rule: FALSE/TRUE if the majority rule idea is applied
  ## ----------------------------------------------------------------------
  ## Value: - decision: res
  ##          res = 1: b is in NO sepset (-> v-structure)
  ##          res = 2: b is in ALL sepsets (-> no v-structure)
  ##          res = 3: b is in SOME but not all sepsets (-> ambiguous triple)
  ##        - version: version (1 or 2) of the ambiguous triple
  ##                (1=normal ambiguous triple that is b is in some sepsets;
  ##                2=triple coming from version.unf[1]==2 that is a and c are
  ##                indep given the initial sepset but there doesn't exist a
  ##                subset of the neighbours that d-separates them)
  ##        - sepsetA and sepsetC: updated separation sets
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 12 Feb 2010, 12:13
  ## Modifications: Diego Colombo, Martin Maechler
  
  
  ## loop through all subsets of parents
 
  
  
  nr.indep <- 0
  stopifnot(length(version.unf) == 2, version.unf %in% 1:2)
  ## Tetrad
  tmp <- if (version.unf[2] == 2) ## our version
    (b %in% sepsetA || b %in% sepsetC) ## else NULL = Tetrad version
  version <- 0
  ## start with the neighbours of a
  if ((nn <- length(nbrsA)) > 0) {
    allComb <- expand.grid(lapply(integer(nn), function(.) 0:1))
    ## loop through all subsets of neighbours
    for (i in 1:nrow(allComb)) { ## == 1:(2^nn)
      S <- nbrsA[which(allComb[i,] != 0)]
      pval <- indepTest(a, c, S, suffStat)
      ## save the pval and the set that produced this pval
      if (verbose) cat("a: S =",S," - pval =",pval,"\n")
      if (pval >= alpha) {
        nr.indep <- nr.indep + 1
        ## is b in set?
        tmp <- c(tmp, b %in% S)
        version <- 1
      }
    }
  }
  ## now with the neighbours of c
  if ((nn <- length(nbrsC)) > 0) {
    allComb <- expand.grid(lapply(integer(nn), function(.) 0:1))
    ## loop through all subsets of neighbours
    for (i in 1:nrow(allComb)) { ## == 1:(2^nn)
      S <- nbrsC[which(allComb[i,] != 0)]
      pval <- indepTest(a, c, S, suffStat)
      ## save the pval and the set that produced this pval
      if (verbose) cat("c: S =",S," - pval =",pval,"\n")
      if (pval >= alpha) {
        nr.indep <- nr.indep + 1
        ## is b in set?
        tmp <- c(tmp, b %in% S)
        version <- 1
      }
    }
  }
  if (version.unf[1] == 2  && nr.indep == 0) {
    version <- 2
  }
  if (is.null(tmp)) tmp <- FALSE
  ratio<- sum(tmp)/length(tmp)
  
  if (verbose) {
    cat("\nCount:",sum(tmp)," N-",length(tmp), "\nprob:",sum(tmp)/length(tmp),"\ratio=:",ratio,"\n")
  }
  
  if (all(tmp)) {
    res <- 2 ## in ALL sets
    ## therefore a - b - c is not a v-structure, hence add b to sepset(a,c)
    ## and sepset(c,a)
    ## for example it can happen that b is not in sepset(a,c) or sepset(c,a)
    ## but now it is in each set that separates a and c given the neighbours
    if (b %nin% sepsetA) sepsetA <- c(sepsetA, b)
    if (b %nin% sepsetC) sepsetC <- c(sepsetC, b)
    
  } else {
    if (all(!tmp)) {
      res <- 1 ## in NO set
      ## therefore a - b - c is a v-structure, hence delete b from
      ## sepset(a,c) and sepset(c,a)
      ## for example it can happen that b is in sepset(a,c) or sepset(c,a)
      ## but now it is in no set that separates a and c given the neighbours
      sepsetA <- setdiff(sepsetA,b)
      sepsetC <- setdiff(sepsetC,b)
     
     
      

      
    } else {
      ## normal conservative PC, b is in some sets hence the triple
      ## is unfaithful
      if (!maj.rule) {
        res <- 3 ## in SOME sets
      } else {
  
        
        ## use the majority rule to test if the triple is faithful
        ## or not and then decide if it is a v-structure or not accordigly
        ## NEW check the percentage of b in the conditioning sets that
        ## make a and c independent
        if (sum(tmp)/length(tmp) < 0.5) {
          ## we accept that b is in NO set
          res <- 1 ## in NO set
          ## therefore a - b - c is a v-structure, hence delete b
          ## from sepset(a,c) and sepset(c,a)
          ## for example it can happen that b is in sepset(a,c) or
          ## sepset(c,a) but now it is in no set that
          ## separates a and c given the neighbours
          sepsetA <- setdiff(sepsetA,b)
          sepsetC <- setdiff(sepsetC,b)
        } else if (sum(tmp)/length(tmp) > 0.5) {
          ## we accept that b is in ALL set
          res <- 2 ## in ALL sets
          ## therefore a - b - c is not a v-structure, hence add b
          ## to sepset(a,c) and sepset(c,a)
          ## for example it can happen that b is not in sepset(a,c)
          ## or sepset(c,a) but now it is in each set that
          ## separates a and c given the neighbours
          if (b %nin% sepsetA) sepsetA <- c(sepsetA,b)
          if (b %nin% sepsetC) sepsetC <- c(sepsetC,b)
        } else if (sum(tmp)/length(tmp) < 0.6 || sum(tmp)/length(tmp) > 0.4) {
          ## define the triple as unfaithful, because half of the
          ## times b is in the set and half of them in not in
          res <- 3 ## in SOME sets
        }
      }
    }
  }
  if (verbose && res == 3) cat("Triple ambiguous\n")
  
  ## if you save a variable <- NULL into a list it will delete this element!
  ## The following also transforms NULL sepset* to integer(0):
  lapply(list(decision = res,ratio = ratio, version = version, SepsetA = sepsetA, SepsetC = sepsetC),
         as.double)
} ## {checkTriple}