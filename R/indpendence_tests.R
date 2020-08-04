partial_correlations <- function(x, y, N_bootstrap=100){
  cor_ts <- cor(x, y, use="complete")
  N_obs <- length(x)
  cor_vec <- c()
  for(i in 1:N_bootstrap){
    cor_vec[i] <- cor(x[sample(N_obs)], y, use="complete")
  }

  p_val <- 1 - mean(cor_ts > cor_vec)
  return(p_val)
}

partial_correlations_CI <- function(x, y, z, N_bootstrap=100){
  fitxz <- lm(x ~ z)
  fityz <- lm(y ~ z)
  x <- fitxz$residuals
  y <- fityz$residuals
  cor_ts <- cor(x, y, use="complete")
  N_obs <- length(x)
  cor_vec <- c()
  for(i in 1:N_bootstrap){
    cor_vec[i] <- cor(x[sample(N_obs)], y, use="complete")
  }
  p_val <- 1 - mean(cor_ts > cor_vec)
  return(p_val)
}


##RCIT
RIT <- function(x, y,
                n_rff=25,
                n_bs=100,
                n_rff_z=NULL,
                polygon_start_index=NULL,
                polygon_end_index=NULL){
  return(RCITcpp::RIT_wrapper(x=x,
                              y=y,
                              n_rff=n_rff,
                              n_bs=100,
                              get_ts=FALSE))
}

RCIT <- function(x, y, z,
                 n_rff=5,
                 n_bs=100,
                 n_rff_z=50,
                 polygon_start_index=NULL,
                 polygon_end_index=NULL){
  return(RCITcpp::RCIT_wrapper(x=x,
                               y=y,
                               z=z,
                               n_rff=n_rff,
                               n_rffz=n_rff_z,
                               n_bs=n_bs,
                               get_ts=FALSE
                               ))
}


###RCIT disaggregation
RIT_disag <- function(x_poly, y, polygon_start_index, polygon_end_index,
                      n_rff=5, n_bs=100, n_rff_z=100){
  return(RCITcpp::RIT_disag_wrapper(x_poly=x_poly, y=y,
                                    polygon_start_index=polygon_start_index,
                                    polygon_end_index=polygon_end_index,
                                    get_ts=FALSE,
                                    n_rff=n_rff))
}

RCIT_disag <- function(x_poly, y, z,
                       polygon_start_index,
                       polygon_end_index,
                       n_rff=5,
                       n_bs=100,
                       n_rff_z=100){
  return(RCITcpp::RCIT_disag_wrapper(x_poly=x_poly, y=y, z=z,
                                     polygon_start_index = polygon_start_index,
                                     polygon_end_index = polygon_end_index,
                                     n_rff=n_rff,
                                     n_rffz=n_rff_z,
                                     n_bs=n_bs,
                                     get_ts=FALSE)
  )
}

