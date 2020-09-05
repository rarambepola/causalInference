pcalg_last_RCIT <- function(obsDat,
                            agg_info,
                            test_opts,
                            alpha,
                            last.index,
                            G_0 = NULL,
                            supprMessages = FALSE,
                            n_cl=NULL,
                            run_parallel=TRUE,
                            aggregated=FALSE,
                            n_subset=100){
  ##messages
  # print(paste0("Running PC algorithm with ", test))
  if(aggregated) print("aggregated")
  if(is.null(n_cl) & run_parallel){
    warning("run_parallel is TRUE but number of cores it not set; running sequentially")
    run_parallel <- FALSE
  }else{
    if(run_parallel) print(paste0("Running in parallel with ", n_cl, " cores")) else print("Running sequentially")
  }

  #choose independence test

  nvar <- length(obsDat)
  #if G_0 is not supplied, make it a complete graph
  if(is.null(G_0)){
    G_0 <- matrix(1, nrow=nvar, ncol=nvar) - diag(nvar)
    warnings("G_0 not supplied")
  }

  #run step 1 of the algorithm to get the skeleton
  if(!supprMessages) print("Starting skeleton algorithm")


  skeleton <- pcskel_last_RCIT(G_0=G_0,
                               obsDat=obsDat,
                               agg_info=agg_info,
                               test_opts=test_opts,
                               alpha=alpha,
                               last.index=last.index,
                               supprMessages=supprMessages,
                               n_cl=n_cl,
                               run_parallel=run_parallel,
                               n_subset=n_subset)
  if(!supprMessages) print("Skeleton found")
  #skeleton matrix
  #print(skeleton)

  colnames(skeleton[[1]]) <- names(obsDat)

  return(list(skeleton[[1]]))
}


##function that produces the skeleton (i.e. undirected acyclic graph) using order-independent
##modifications to the algorithm.
pcskel_last_RCIT <- function(G_0,
                             obsDat,
                             agg_info,
                             test_opts,
                             alpha,
                             last.index,
                             supprMessages,
                             n_cl,
                             run_parallel,
                             n_subset){
  `%do_switch%` <- ifelse(run_parallel, `%dopar%`, `%do%`)
  #G_old is the graph we will use for the adj sets during each loop
  #G_new is the graph with edges deleted as we go
  #need both since G_1 is only updated to G_2 each time the size of the
  #conditioning set grows
  G_old <- G_0
  G_new <- G_0

  #S is the list of seperation sets
  S <- list()

  ##Step 1:Test for unconditional independence
  #get list of edges only connected to last element

  getEdgesLast <- function(G_0, last.index){
    adj.nodes <- unique(c(which(G_0[last.index, ] > 0), which(G_0[, last.index] > 0)))
    #print(G_0)
    #print(adj.nodes)
    edges <- list()
    if(length(adj.nodes) == 0) return(NULL)
    for(i in 1:length(adj.nodes)){
      edges[[i]] <- c(last.index, adj.nodes[i])
    }
    return(edges)
  }

  edges <- getEdgesLast(G_0, last.index)

  #print(edges)

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
               .packages=c('RCITcpp')) %do_switch%{
                 edge <- edges[[i]]

                 testDat <- RCITcpp::format_vars(obsDat[edge],
                                        agg_info$agg_levels[edge],
                                        agg_info$poly_static,
                                        agg_info$poly_dynamic,
                                        agg_info$pixel_static,
                                        agg_info$pixel_dynamic,
                                        agg_info$times_dynamic,
                                        agg_info$inc_times,
                                        agg_info$inc_poly,
                                        agg_info$pop_dynamic,
                                        agg_info$pop_static,
                                        n_subset)

                 if(RCITcpp::RIT_disag_wrapper_v2(testDat$obs[[1]],
                                         testDat$obs[[2]],
                                         testDat$start_index,
                                         testDat$end_index,
                                         testDat$pop,
                                         test_opts$n_rff,
                                         test_opts$n_bs) < alpha){

                 # if(independence_test_wrapper(obsDat[[edge[1]]], obsDat[[edge[2]]], alpha=alpha, independence_test, ...)){
                   edgeremove <- NULL
                   sepset <- NULL
                 }else{
                   edgeremove <- edge
                   sepset <- NA
                 }
                 list(edgeremove, sepset)
               }


  edgeremove <- z[[1]]
  #print(edgeremove)
  sepsets <- z[[2]]
  keep <- !unlist(lapply(edgeremove, is.null))
  edgeremove <- edgeremove[keep]
  sepsets <- sepsets[keep]
  #print(edgeremove)
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

  if(length(getEdgesLast(G_old, last.index)) == 0){
    if(run_parallel) stopCluster(cl)
    return(list(G_old, S))
  }


  ##get size of the maximum conditioning set...how long does it take??
  ptm <- proc.time()
  edges <- getEdgesLast(G_old, last.index)

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
  edges <- getEdgesLast(G_old, last.index)
  #for each edge, check if the adjacency set is big enough (both directions)
  #to condition on
  #variable maxCondSize checks if the biggest possible size has been reached for
  #the conditioning sets
  belowMaxCondSize <- TRUE
  while(belowMaxCondSize){
    belowMaxCondSize <- FALSE
    print(paste0("condSetSize ", condSetSize))
    #go through each edge


    z <- foreach(i = 1:length(edges), .combine = comb,
                 .init = list(list(), list()),
                 .packages=c('RCITcpp'),
                 .export=c('getAdj', 'getSubsets')) %do_switch%{
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
                       testDat <- RCITcpp::format_vars(obsDat[c(i, j, condSet)],
                                              agg_info$agg_levels[c(i, j, condSet)],
                                              agg_info$poly_static,
                                              agg_info$poly_dynamic,
                                              agg_info$pixel_static,
                                              agg_info$pixel_dynamic,
                                              agg_info$times_dynamic,
                                              agg_info$inc_times,
                                              agg_info$inc_poly,
                                              agg_info$pop_dynamic,
                                              agg_info$pop_static,
                                              n_subset)

                       condData <- do.call(cbind, testDat$obs[-(1:2)])
                       #if they are conditionally independent, remove edge, add condSet to
                       #separation set list and stop trying this edge

                       if(RCITcpp::RCIT_disag_wrapper_v2(testDat$obs[[1]], testDat$obs[[2]], condData,
                                                testDat$start_index,
                                                testDat$end_index,
                                                testDat$pop,
                                                test_opts$n_rff,
                                                test_opts$n_rff_z,
                                                test_opts$n_bs) > alpha){


                       # #make a matrix of conditioning data
                       # condData <- do.call(cbind, obsDat[condSet])
                       # #if they are conditionally independent, remove edge, add condSet to
                       # #separation set list and stop trying this edge
                       # if(!CI_test_wrapper(obsDat[[i]], obsDat[[j]], condData, alpha=alpha, CI_test, ...)){

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
                         testDat <- RCITcpp::format_vars(obsDat[c(i, j, condSet)],
                                                         agg_info$agg_levels[c(i, j, condSet)],
                                                         agg_info$poly_static,
                                                         agg_info$poly_dynamic,
                                                         agg_info$pixel_static,
                                                         agg_info$pixel_dynamic,
                                                         agg_info$times_dynamic,
                                                         agg_info$inc_times,
                                                         agg_info$inc_poly,
                                                         agg_info$pop_dynamic,
                                                         agg_info$pop_static,
                                                         n_subset)

                         condData <- do.call(cbind, testDat$obs[-(1:2)])
                         #if they are conditionally independent, remove edge, add condSet to
                         #separation set list and stop trying this edge

                         if(RCITcpp::RCIT_disag_wrapper_v2(testDat$obs[[1]], testDat$obs[[2]], condData,
                                                           testDat$start_index,
                                                           testDat$end_index,
                                                           testDat$pop,
                                                           test_opts$n_rff,
                                                           test_opts$n_rff_z,
                                                           test_opts$n_bs) > alpha){

                         #if they are conditionally independent, remove edge, add condSet to
                         #separation set list and stop trying this edge
                         # if(!CI_test_wrapper(obsDat[[i]], obsDat[[j]], condData, alpha=alpha, CI_test, ...)){
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
    #print("z[[1]]")
    #print(z[[1]])
    sepsets <- z[[2]]
    keep <- !unlist(lapply(edgeremove, is.null))
    edgeremove <- edgeremove[keep]
    sepsets <- sepsets[keep]
    #print(edgeremove)
    if(length(edgeremove) > 0){
      #print(edgeremove)
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
    edges <- getEdgesLast(G_old, last.index)
    if(length(edges) == 0) return(list(G_old, S))

    if(!is.null(edges)){
      maxsize <- 0
      for(edge in edges){
        i <- edge[1]
        j <- edge[2]
        adjSet <- getAdj(G_old, i, j)
        maxsize <- max(maxsize, length(adjSet))
      }
      print(paste0("maxsize ", maxsize))
      belowMaxCondSize <- maxsize >= condSetSize
    }
  }
  if(run_parallel) stopCluster(cl)
  colnames(G_old) <- names(obsDat)
  return(list(G_old, S))
}
