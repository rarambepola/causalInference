

pcalg_both_steps_RCIT <- function(obsDat,
                                  alpha,
                                  last.index,
                                  agg_info,
                                  test_opts,
                                  G_0 = NULL,
                                  supprMessages = FALSE,
                                  n_cl=NULL,
                                  run_parallel=TRUE,
                                  subset_step_1=100,
                                  subset_step_2=100,
                                  run_parallel_step_2=NULL,
                                  alpha_step_2=NULL
){
  if(is.null(alpha_step_2)) alpha_step_2 <- alpha
  nvar <- length(obsDat)
  if(is.null(G_0)){
    G_0 <- matrix(1, nrow=nvar, ncol=nvar) - diag(nvar)
  }
  if(is.null(run_parallel_step_2)) run_parallel_step_2 <- run_parallel

  agg_info1 <- agg_info
  agg_info1$agg_levels <- agg_info1$agg_levels[-last.index]

  pc_step_1 <- pcalg_RCIT(obsDat = obsDat[-last.index],
                     agg_info = agg_info1,
                     test_opts = test_opts,
                     alpha = alpha,
                     G_0 = G_0[-last.index, -last.index],
                     supprMessages=supprMessages,
                     n_cl = n_cl,
                     run_parallel = run_parallel,
                     n_subset=subset_step_1)

  G_1 <- pc_step_1[[1]]
  G_0[-last.index, -last.index] <- G_1
  #deal with subsets in step 2
  pc_step_2 <- pcalg_last_RCIT(obsDat = obsDat,
                          agg_info = agg_info,
                          alpha = alpha_step_2,
                          last.index = last.index,
                          G_0=G_0,
                          supprMessages=supprMessages,
                          n_cl=n_cl,
                          run_parallel=run_parallel_step_2,
                          n_subset = subset_step_2,
                          test_opts = test_opts)
  if(!is.null(names(obsDat))){
    colnames(pc_step_2[[1]]) <- names(obsDat)
  }
  return(pc_step_2)
}
