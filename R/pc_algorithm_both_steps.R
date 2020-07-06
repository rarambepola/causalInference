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
                             ...){
  if(is.null(alpha_step_2)) alpha_step_2 <- alpha
  nvar <- length(obsDat)
  if(is.null(G_0)){
    G_0 <- matrix(1, nrow=nvar, ncol=nvar) - diag(nvar)
  }

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

  pc_step_2 <- pcalg_last(obsDat=obsDat,
                          alpha=alpha_step_2,
                          last.index=last.index,
                          test=test,
                          G_0=G_0,
                          supprMessages=supprMessages,
                          n_cl=n_cl,
                          run_parallel=run_parallel,
                          aggregated=aggregated,
                          ...)
  if(!is.null(names(obsDat))){
    colnames(pc_step_2[[1]]) <- names(obsDat)
  }
  return(pc_step_2)
}
