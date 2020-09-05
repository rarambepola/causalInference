reaggregate <- function(cov_list, polygon_start, polygon_end){
  n_obs <- length(polygon_start)
  cov_list_out <- lapply(cov_list, function(l) lapply(1:n_obs, function(i) return(l[polygon_start[i]:polygon_end[i]])))
  return(cov_list_out)
}

pcalg_both_steps <- function(obsDat,
                             alpha,
                             last.index,
                             test="RCIT",
                             G_0 = NULL,
                             supprMessages = FALSE,
                             n_cl=NULL,
                             run_parallel=TRUE,
                             aggregated=FALSE,
                             alpha_step_2=NULL,
                             subset_step_1=NULL,
                             subset_step_2=NULL,
                             run_parallel_step_2=NULL,
                             pop_weighting=FALSE,
                             ...){
  if(is.null(alpha_step_2)) alpha_step_2 <- alpha
  nvar <- length(obsDat)
  if(is.null(G_0)){
    G_0 <- matrix(1, nrow=nvar, ncol=nvar) - diag(nvar)
  }
  if(is.null(run_parallel_step_2)) run_parallel_step_2 <- run_parallel
  obsDat_step1 <- obsDat[-last.index]
  if(!is.null(subset_step_1)){
    sample_use <- sample.int(length(obsDat_step1[[1]]), subset_step_1)
    obsDat_step1 <- lapply(obsDat[-last.index], function(l) l[sample_use])
  }

  pc_step_1 <- pcalg(obsDat = obsDat_step1,
                     alpha = alpha,
                     test = test,
                     G_0 = G_0[-last.index, -last.index],
                     supprMessages=supprMessages,
                     n_cl = n_cl,
                     run_parallel = run_parallel,
                     ...)

  G_1 <- pc_step_1[[1]]
  G_0[-last.index, -last.index] <- G_1

  #deal with subsets in step 2
  n_obs <- length(obsDat[[1]])
  if(!is.null(subset_step_2)){
    n_sample <- min(n_obs, subset_step_2)
    sample_use <- sample.int(n_obs, n_sample)
    if(aggregated){
      #undo aggregation
      obsDat[-last.index] <- reaggregate(obsDat[-last.index], polygon_start_index + 1, polygon_end_index + 1)
      if(pop_weighting) population <- reaggregate(population, polygon_start_index + 1, polygon_end_index + 1)
      obsDat <- lapply(obsDat, function(l) l[sample_use])
      polygon_sizes <- cumsum(sapply((obsDat[-last.index])[[1]], length))
      polygon_start_index <- c(0, polygon_sizes)
      polygon_start_index <- polygon_start_index[-length(polygon_start_index)]
      polygon_end_index <- polygon_sizes - 1
      obsDat <- lapply(obsDat, unlist)
      if(pop_weighting) population <- unlist(population)
    }else{
      obsDat <- lapply(obsDat, function(l) l[sample_use])
    }
  }

  pc_step_2 <- pcalg_last(obsDat=obsDat,
                          alpha=alpha_step_2,
                          last.index=last.index,
                          test=test,
                          G_0=G_0,
                          supprMessages=supprMessages,
                          n_cl=n_cl,
                          run_parallel=run_parallel_step_2,
                          aggregated=aggregated,
                          ...)
  if(!is.null(names(obsDat))){
    colnames(pc_step_2[[1]]) <- names(obsDat)
  }
  return(pc_step_2)
}
