####file containing  pc algorithm

###first version is following the description in the pcalg package documentation
##G_0 is the initial adjacency matrix - usually a complete graph but made an input
##for added flexibility
pcalg <- function(obsDat,
                  alpha,
                  test="RCIT",
                  G_0 = NULL,
                  supprMessages = FALSE,
                  n_cl=NULL,
                  run_parallel=TRUE,
                  ...){
  ##messages
  print(paste0("Running PC algorithm with ", test))
  if(is.null(n_cl) & run_parallel){
    warning("run_parallel is TRUE but number of cores it not set; running sequentially")
    run_parallel <- FALSE
  }else{
    if(run_parallel) print(paste0("Running in parallel with ", n_cl, " cores")) else print("Running sequentially")
  }

  #choose independence test
  if(test == "RCIT"){
    independence_test <- RIT
    CI_test <- RCIT
  }
  if(test == "partial correlation"){
    independence_test <- partial_correlations
    CI_test <- partial_correlations_CI
  }

  nvar <- length(obsDat)

  #if G_0 is not supplied, make it a complete graph
  if(is.null(G_0)){
    G_0 <- matrix(1, nrow=nvar, ncol=nvar) - diag(nvar)
  }

  #run step 1 of the algorithm to get the skeleton
  print("Starting skeleton algorithm")
  skeleton <- pcskel(G_0, obsDat, alpha,
                     independence_test, CI_test,
                     supprMessages,
                     n_cl=n_cl,
                     run_parallel=run_parallel, ...)
  print("Skeleton found")
  #skeleton matrix
  G_skel <- skeleton[[1]]
  G_out <- G_skel
  G_temp <- G_out
  #separating sets
  S <- skeleton[[2]]

  #get all unshielded triples i-j-k, for each orient them if
  #j is not in the sepset of i,k
  #actually keep track of any that should be oriented and orient them
  #afterwards if there are no contradictions
  triples <- getUnTrpls(G_skel)

  #tripleMat is a matrix for keeping track of which edges should be oriented
  tripleMat <- matrix(0, ncol=nvar, nrow=nvar)
  for(triple in triples){
    i <- triple[1]
    j <- triple[2]
    k <- triple[3]
    #if j is not in sepset orient
    if(!(j %in% S[[deparse(c(i, k))]])){
      tripleMat[i, j] <- 1
      tripleMat[k, j] <- 1
    }
  }

  #go through edges in tripleMat, if i-j is > 0 but j-i is 0 there is no disagreement
  #and we can orient
  for(i in 1:nvar){
    for(j in 1:nvar){
      if(tripleMat[i, j] > 0 & tripleMat[j, i] == 0){
        G_temp[j, i] <- 0
      }
    }
  }

  ##apply rules 1-3 repeatedly until nothing changes
  changed <- TRUE
  while(changed){
    rule1.out <- .rule1(G_temp)
    rule2.out <- .rule2(rule1.out[[1]])
    rule3.out <- .rule3(rule2.out[[1]])
    G_temp <- rule3.out[[1]]
    changed <- any(c(rule1.out[[2]], rule2.out[[2]], rule3.out[[2]]))
  }

  ##if at any time we deleted both directions of an edge, this is wrong
  ##so put them back (unless they were not there in the original
  ##graph - these are the extra if statements inside)
  for(i in 1:nvar){
    for(j in 1:i){
      if(anyEdge(G_out, i, j) && !anyEdge(G_temp, i, j)){
        if(G_0[i, j] == 1) G_temp[i, j] <- 1
        if(G_0[j, i] == 1) G_temp[j, i] <- 1
      }
    }
  }

  G_out <- G_temp

  print("CPDAG found")
  return(list(G_out, G_skel, S))
}



##function that produces the skeleton (i.e. undirected acyclic graph) using order-independent
##modifications to the algorithm.
pcskel <- function(G_0, obsDat, alpha, independence_test, CI_test, supprMessages = FALSE,
                   n_rff=5, n_rffz=5, n_cl=NULL, run_parallel=TRUE, ...){
  #G_old is the graph we will use for the adj sets during each loop
  #G_new is the graph with edges deleted as we go
  #need both since G_1 is only updated to G_2 each time the size of the
  #conditioning set grows
  G_old <- G_0
  G_new <- G_0

  `%do_switch%` <- ifelse(run_parallel, `%dopar%`, `%do%`)

  #S is the list of seperation sets
  S <- list()

  ##Step 1:Test for unconditional independence
  #get list of edges from adjacency matrix
  edges <- getEdges(G_old)

  if(run_parallel){
    cl<-makeCluster(n_cl)
    registerDoParallel(cl)
  }

  comb <- function(x, ...) {
    lapply(seq_along(x),
           function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
  }

  z <- foreach(i = 1:length(edges), .combine = comb,
               .init = list(list(), list()),
               .packages=c('RCITcpp'),
               .export=c('independence_test_wrapper')) %do_switch%{
                 edge <- edges[[i]]
                 if(independence_test_wrapper(obsDat[[edge[1]]], obsDat[[edge[2]]], alpha=alpha, independence_test, ...)){
                   # if(fstest(obsDat[[edge[1]]], obsDat[[edge[2]]], alpha=alpha, test=TRUE,
                   #           n_rff=n_rff)){
                   edgeremove <- NULL
                   sepset <- NULL
                 }else{
                   edgeremove <- edge
                   sepset <- NA
                 }
                 list(edgeremove, sepset)
               }

  edgeremove <- z[[1]]
  sepsets <- z[[2]]
  keep <- !unlist(lapply(edgeremove, is.null))
  edgeremove <- edgeremove[keep]
  sepsets <- sepsets[keep]
  if(length(edgeremove) > 0){
    for(k in 1:length(edgeremove)){
      edge <- edgeremove[[k]]
      sepset <- sepsets[[k]]
      i <- edge[1]
      j <- edge[2]
      G_new[i, j] <- 0
      G_new[j, i] <- 0
      S[[deparse(c(i, j))]] <- NA
      S[[deparse(c(j, i))]] <- NA
    }
  }

  #update adjacency matrix (actually could have just updated G_old for this
  #whole step, done to make it conceptually similar to later)
  G_old <- G_new

  if(length(getEdges(G_old)) == 0){
    if(run_parallel) stopCluster(cl)
    return(list(G_old, S))
  }

  ##get size of the maximum conditioning set...how long does it take??
  ptm <- proc.time()
  edges <- getEdges(G_old)
  maxsize <- 0
  belowMaxCondSize <- TRUE
  for(edge in edges){
    i <- edge[1]
    j <- edge[2]
    adjSet <- getAdj(G_old, i, j)
    maxsize <- max(maxsize, length(adjSet))
  }

  ##Step 2: Test for conditional independence
  condSetSize <- 1
  edges <- getEdges(G_old)
  #for each edge, check if the adjacency set is big enough (both directions)
  #to condition on
  #variable maxCondSize checks if the biggest possible size has been reached for
  #the conditioning sets
  belowMaxCondSize <- TRUE
  while(belowMaxCondSize){
    belowMaxCondSize <- FALSE
    if(!supprMessages) print(paste0("condSetSize ", condSetSize))
    #go through each edge


    z <- foreach(i = 1:length(edges), .combine = comb,
                 .init = list(list(), list()),
                 .packages=c('RCITcpp'),
                 .export=c('getAdj', 'getSubsets', 'CI_test_wrapper')) %do_switch% {
                   edge <- edges[[i]]
                   if(!supprMessages) print(paste0("testing edge ", edge[1], ",", edge[2], " with size ", condSetSize))
                   i <- edge[1]
                   j <- edge[2]
                   brokenEdge <- FALSE
                   #do the i, j direction
                   #get set adj(G,i)\j
                   adjSet <- getAdj(G_old, i, j)

                   #if the adjacency set is bigger than conditioning set size, test
                   #conditional independence
                   if(length(adjSet) >= condSetSize){
                     belowMaxCondSize <- TRUE
                     #test independence conditioning on  all subsets of right size
                     condSets <- getSubsets(adjSet, condSetSize)

                     for(condSet in condSets){
                       #make a matrix of conditioning data
                       condData <- do.call(cbind, obsDat[condSet])
                       #if they are conditionally independent, remove edge, add condSet to
                       #separation set list and stop trying this edge
                       if(!CI_test_wrapper(obsDat[[i]], obsDat[[j]], condData, alpha=alpha, CI_test, ...)){
                         # if(!condtest(obsDat[[i]], obsDat[[j]], condData, alpha=alpha,
                         #              n_rff=n_rff, n_rffz=n_rffz)){
                         edgeremove <- edge
                         sepset <- condSet
                         brokenEdge <- TRUE
                         break
                       }
                     }
                   }

                   if(!brokenEdge){
                     #do j, i direction
                     i <- edge[1]
                     j <- edge[2]

                     #do the i, j direction
                     #get set adj(G,i)\j
                     adjSet <- getAdj(G_old, j, i)

                     #if the adjacency set is bigger than conditioning set size, test
                     #conditional independence
                     if(length(adjSet) >= condSetSize){
                       #test independence conditioning on  all subsets of right size
                       condSets <- getSubsets(adjSet, condSetSize)

                       for(condSet in condSets){
                         #make a matrix of conditioning data
                         condData <- do.call(cbind, obsDat[condSet])
                         #if they are conditionally independent, remove edge, add condSet to
                         #separation set list and stop trying this edge
                         if(!CI_test_wrapper(obsDat[[i]], obsDat[[j]], condData, alpha=alpha, CI_test, ...)){
                           # if(!condtest(obsDat[[i]], obsDat[[j]], condData, alpha=alpha)){
                           edgeremove <- edge
                           sepset <- condSet
                           brokenEdge <- TRUE
                           break
                         }
                       }
                     }
                   }
                   if(!brokenEdge){
                     edgeremove <- NULL
                     sepset <- NULL
                   }
                   list(edgeremove, sepset)
                 }
    #update the adjacency matrix with edges remove

    edgeremove <- z[[1]]
    # print(length(z[[1]]))
    sepsets <- z[[2]]
    # print(length(z[[2]]))
    # print(edgeremove)s
    keep <- !unlist(lapply(edgeremove, is.null))
    edgeremove <- edgeremove[keep]
    sepsets <- sepsets[keep]
    if(length(edgeremove) > 0){
      for(k in 1:length(edgeremove)){
        edge <- edgeremove[[k]]
        sepset <- sepsets[[k]]
        i <- edge[1]
        j <- edge[2]
        G_new[i, j] <- 0
        G_new[j, i] <- 0
        S[[deparse(c(i, j))]] <- sepset
        S[[deparse(c(j, i))]] <- sepset
      }
    }
    G_old <- G_new
    condSetSize <- condSetSize + 1
    edges <- getEdges(G_old)
    maxsize <- 0
    for(edge in edges){
      i <- edge[1]
      j <- edge[2]
      adjSet <- getAdj(G_old, i, j)
      maxsize <- max(maxsize, length(adjSet))
    }
    if(!supprMessages) print(paste0("maxsize ", maxsize))
    belowMaxCondSize <- maxsize >= condSetSize
  }
  if(run_parallel) stopCluster(cl)
  colnames(G_old) <- names(obsDat)
  return(list(G_old, S))
}

