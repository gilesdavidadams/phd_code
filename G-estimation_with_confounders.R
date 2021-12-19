rm(list=ls())
gc()
library(matlib)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(foreach)
library(doParallel)
library(tidyverse)


create_sample <- function(prm){
  
  psi_star_0 <- prm$psi_star_0
  psi_star_x <- prm$psi_star_x
  psi_star_vec <- prm$psi_star_vec
  
  prb_x = prm$prb_x
  prb_a_x0 = prm$prb_a_x0
  prb_a_x1 = prm$prb_a_x1
  
  n <- prm$n
  t0_min <- prm$t0_min
  t0_max <- prm$t0_max
  
  # creates treatment variables
  x0_size <- floor(n*prb_x)
  x1_size <- n - floor(n*prb_x)
  x <- c(rep(0, x0_size), 
         rep(1, x1_size))
  
  a0_xeq0 <-  sample(c(rep(0, x0_size - floor(x0_size*prb_a_x0)),
              rep(1, floor(x0_size*prb_a_x0))),
              size=x0_size)
  a0_xeq1 <-  sample(c(rep(0, x1_size - floor(x1_size*prb_a_x1)),
                     rep(1, floor(x1_size*prb_a_x1))),
                   size=x1_size)
  
  df <- tibble(x=x, a_0 = c(a0_xeq0, a0_xeq1))
 
  df <- df %>% add_column(t0  = runif(n, min=t0_min, max=t0_max )) %>%
    #mutate(ti = 0, t0_rsd = t0) %>%
    mutate(ti = exp((psi_star_vec[[1]] + psi_star_vec[[2]]*x)*a_0)*t0)
  
  return(df)
  #return(select(df, -c(t0_rsd, temp, w_TOT, g_psi)))
}

fit_treatment_models <- function(df, prm){
  
  df <- df %>% mutate(fit_a0 = glm(a_0 ~ x, family=binomial)$fitted.values)
  return(df)
}

calculate_tau_k <- function(df, prm, psi_hat_vec,  a_k=0, ...){
  
  psi_hat_a0 <- psi_hat_vec[[1]]
  psi_hat_x <- psi_hat_vec[[2]]
  
  df <- df %>% mutate(t0_hat_k = exp(-(psi_hat_a0 + psi_hat_x*x)*a_0)*ti)
  
  return(df)
}

calculate_score <- function(df, prm, psi_hat_vec){
  
  df_temp <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec)
  
  S_a0 <- with(df_temp, sum((a_0 - fit_a0)*t0_hat_k))
  S_x <- with(df_temp, sum((a_0 - fit_a0)*t0_hat_k*x))
  
  S_vec <- c(S_a0, S_x)
  
  return(S_vec)
}

calculate_jacobian <- function(df, prm, psi_hat_vec){
  
  df_temp <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec)
  
  jacobi_vec <- c(
    with(df_temp, sum((a_0 - fit_a0)*t0_hat_k*a_0)),
    with(df_temp, sum((a_0 - fit_a0)*t0_hat_k*x*a_0)),
    with(df_temp, sum((a_0 - fit_a0)*t0_hat_k*x*a_0)),
    with(df_temp, sum((a_0 - fit_a0)*t0_hat_k*x*x*a_0))
  )
  
  return(matrix(jacobi_vec, byrow=T, nrow=2))
  #return(matrix(jacobi_vec, byrow=T, nrow= (length(psi_hat_CT) + length(psi_hat_TOT) - 1)))
}

calculate_variance <- function(df, psi_hat_vec, prm){
  
  D_vec <- c(rep(1, prm$n), df$x)
  D <- matrix(D_vec, ncol=2, byrow=F)
  
  df_temp <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec)
  
  D_psi_vec <- c(df_temp$t0_hat_k, df_temp$t0_hat_k*df_temp$x)
  D_psi <- matrix(D_psi_vec, byrow=F, ncol=2)
  
  fit <- df$fit_a0
  
  #df_temp <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec)
  #tau_vec <- df_temp$t0_hat_k
  
  Jbb_vec <- c()
  for(i in 1:dim(D)[2]){ 
    for(j in 1:dim(D)[2]){
      Jbb_vec <- c(Jbb_vec, sum(D[,i]*D[,j]*fit*(1-fit)))
    }
  }
  Jbb <- matrix(Jbb_vec, nrow=sqrt(length(Jbb_vec)), byrow=T)
  
  
  Jtt <- t(sapply(1:dim(D_psi)[[2]], function(i){
    sapply(1:dim(D_psi)[[2]], function(j){
      sum(D_psi[,i]*D_psi[,j]*fit*(1-fit))
    })
  }))
  
  Jtb <- t(sapply(1:dim(D_psi)[[2]], function(i){
    sapply(1:dim(D)[[2]], function(j){
      sum(D_psi[,i]*D[,j]*fit*(1-fit))
    })
  }))
  
  
  JTT <- Jtt - Jtb%*%solve(Jbb)%*%t(Jtb)
  
  Jacobian <- df %>% calculate_jacobian(prm=prm, psi_hat_vec=psi_hat_vec)
  J_inv <- solve(Jacobian)
  
  return(J_inv %*% JTT %*% t(J_inv))
}





newton_raphson_grad <- function(df, prm,
                                psi_start_vec=c(0,0), 
                                tol=0.001, max_iter=20, psi_max = 10){

  psi_hat_vec <- psi_start_vec
  psi_old_vec <- psi_hat_vec
  steps <- 1L
  fail_flag <- F
  
  while(((sum(abs(psi_hat_vec - psi_old_vec)) > tol) | steps == 1) &
        (steps <= max_iter) & (max(abs(psi_hat_vec)) < psi_max)){
    
    psi_old_vec <- psi_hat_vec
    
    (S_vec <- df %>% calculate_score(prm=prm, psi_hat_vec=psi_hat_vec))
    (J <- df %>% calculate_jacobian(prm=prm, psi_hat_vec=psi_hat_vec))
    
    (psi_hat_vec <- solve(J, S_vec) + psi_hat_vec)
    
    steps <- steps + 1
  }
  
  if (max(abs(psi_hat_vec)) > psi_max){
    fail_flag <- T
  }
  return(list(psi_hat_vec, steps, fail_flag))
}




nr_run <- function(psi_star_0, 
                   psi_star_x,
                   prb_x = 0.5,
                   prb_a_x0 = 0.4,
                   prb_a_x1 = 0.7,
                   n = 1000,
                   t0_min = 10,
                   t0_max = 100,
                   sims = 10){
  
  
  prm <- list(psi_star_0=psi_star_0,
              psi_star_x=psi_star_x,
              psi_star_vec=c(psi_star_0, psi_star_x),
              prb_x = prb_x,
              prb_a_x0 = prb_a_x0,
              prb_a_x1 = prb_a_x1,
              #psi_star_list=list(psi_star_CT, psi_star_TOT),
              #t_a_vec=t_a_vec,
              n=n,
              t0_min=t0_min,
              t0_max=t0_max,
              sims=sims)
  
  #prm$t_df <- get_time_df(prm=prm)
  
  cl <- makePSOCKcluster(28)
  registerDoParallel(cl)
  #registerDoParallel()
  
  nr_out <- foreach(icount(sims), .combine=rbind,
                    .export=c("create_sample", "fit_treatment_models",
                              "calculate_tau_k",
                              "calculate_score", "calculate_jacobian", 
                             "newton_raphson_grad",
                              "calculate_variance"
                    ), .packages="tidyverse") %dopar%
    {
      df <- create_sample(prm=prm)
      df <- df %>% fit_treatment_models(prm=prm)
      
      nri_out <- df %>% newton_raphson_grad(prm=prm)
      psi_hat_vec <- nri_out[[1]]
      
      var_hat <- df %>% calculate_variance(psi_hat_vec=psi_hat_vec, prm=prm)
      
      c(psi_hat_vec, diag(var_hat))
      #c(unlist(psi_hat_list), diag(var_hat))
    }
  
  stopCluster(cl)
  
  for (j in 1:2){
    if (j == 1){ psi_lab <- "psi_a0"}
    else {psi_lab <- "psi_x"}

    cat(
      paste(psi_lab, "\t", "MEAN", "\t\t", "VAR", collapse="\t"),
      "\n",
      paste(
        c("TRU",
        format(c(prm$psi_star_vec[[j]],
                 var(nr_out[,j], na.rm=T))
               , nsmall=5)),
        collapse = "\t"),
      "\n",
      paste(
        c("EST",
        format(c(mean(nr_out[,j] , na.rm=T),
                 mean(nr_out[,j+2], na.rm=T))
               , nsmall=5)),
        collapse = "\t"),
      "\n")
  }
  
  #return(NULL)
}

nr_run(psi_star_0 = c(log(2)), 
       psi_star_x = c(log(1.5)),
       prb_x = 0.5,
       prb_a_x0 = 0.4,
       prb_a_x1 = 0.6,
       #t_a_vec  = c(0, 30),
       n = 40,
       sims=1000
)
 
