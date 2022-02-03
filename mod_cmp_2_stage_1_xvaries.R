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
  with(prm, {
    
    a0_size <- max(floor(n/2),1)
    a1_size <- max(floor(n/4),1)
    x0_size <-  max(floor(n/8),1)
    x1_size <- max(floor(n/16),1)
    
    df <- tibble(a_0 = c(rep(0, a0_size), rep(1, n - a0_size))) %>%
      
      add_column(a_1 = c(rep(0, a1_size), rep(1, a0_size - a1_size),
                         rep(0, a1_size), rep(1, n - (a0_size + a1_size)))) %>%
      
      add_column(x_0 = c(rep(0, x0_size), rep(1, a1_size - x0_size),
                       rep(0, x0_size), rep(1, (a0_size - a1_size) - x0_size),
                       rep(0, x0_size), rep(1, a1_size - x0_size),
                       rep(0, x0_size), rep(1, n - (a0_size + a1_size + x0_size)))) %>%
      
      add_column(x_1 = c(rep(0, x1_size), rep(1, x0_size - x1_size),
                        rep(0, x1_size), rep(1, a1_size - x0_size - x1_size),
                        rep(0, x1_size), rep(1, x0_size - x1_size),
                        rep(0, x1_size), rep(1, a0_size - a1_size - x0_size - x1_size),
                        rep(0, x1_size), rep(1, x0_size - x1_size),
                        rep(0, x1_size), rep(1, a1_size - x0_size - x1_size),
                        rep(0, x1_size), rep(1, x0_size - x1_size),
                        rep(0, x1_size), rep(1, n - a0_size - a1_size - x0_size - x1_size))) %>%
      
      add_column(t0  = runif(n, min=t0_min, max=t0_max )) %>%
      
      mutate(ti = 0, t0_rsd = t0,
             g_psi = psi_star_vec[[1]] + x_0*psi_star_vec[[2]],
             temp = pmin((t_a - 0)*exp(-g_psi*a_0), t0_rsd),
             ti = ti + temp*exp(g_psi*a_0),
             t0_rsd = t0_rsd - temp,
             
             g_psi = psi_star_vec[[1]] + x_1*psi_star_vec[[2]],
             temp = t0_rsd,
             ti = ti + temp*exp(g_psi*a_1),
             t0_rsd = t0_rsd - temp
      ) %>% select(-c(t0_rsd, g_psi, temp))
    return(df)
  })
}

fit_treatment_models <- function(df, prm){
  with(prm, {
    df <- df %>% mutate(fit_a0 = glm(a_0 ~ x_0, family=binomial)$fitted.values)
    
    mdl_formula <- "a_1 ~ a_0 + x_1"
    model_temp <- with(df %>% filter(.$ti > t_a), 
                       glm(  as.formula(mdl_formula),
                             family=binomial))
    df$fit_a1 <- predict(model_temp, newdata = df, type="response")
    return(df)
  })
}


calculate_tau_k <- function(df, prm, psi_hat_vec,  a_k=0, ...){
  with(prm, {
    df <- df %>% mutate(ti_temp = pmax(ti - t_a, 0),
                        g_psi = psi_hat_vec[[1]] + psi_hat_vec[[2]]*x_1, 
                        tau_k = exp(-g_psi*a_1)*ti_temp)
    if (a_k == 0) {
      df <- df %>% mutate(ti_temp = pmax(pmin(t_a, ti), 0),
                          g_psi = psi_hat_vec[[1]] + psi_hat_vec[[2]]*x_0, 
                          tau_k = tau_k + exp(-g_psi*a_0)*ti_temp)
    } else {
      #a_k == 1
      df <- df %>% filter(.$ti > t_a)
    }
    
    df <- df %>% select(-c(ti_temp))
    return(df)
  })
}

calculate_tau_rsd_m <- function(df, prm, 
                                psi_hat_vec,
                                m=0, ...){
  with(prm, {
    if (m == 0) {
      df <- df %>% mutate(g_psi = psi_hat_vec[[1]] + psi_hat_vec[[2]]*x_0,
                          ti_temp = pmax(0, pmin(t_a, ti) - 0),
                          tau_rsd = ti_temp*exp(-g_psi*a_0)
      )
    } else {
      df <- df %>% mutate(g_psi = psi_hat_vec[[1]] + psi_hat_vec[[2]]*x_1,
                          ti_temp = pmax(0, ti - t_a),
                          tau_rsd = ti_temp*exp(-g_psi*a_1) 
      )
    }
    
    df <- df %>% select(-c(g_psi, ti_temp))
    
    return(df )
  })
  
}

calculate_score <- function(df, prm, psi_hat_vec){
  
  df_1 <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=1)
  S_0 <- with(df_1, sum((a_1 - fit_a1)*tau_k))
  S_1 <- with(df_1, sum((a_1 - fit_a1)*tau_k*x_1))
  
  df_0 <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=0)
  S_0 <- S_0 + with(df_0, sum((a_0 - fit_a0)*tau_k))
  S_1 <- S_1 + with(df_0, sum((a_0 - fit_a0)*tau_k*x_0))
  
  S_vec <- c(S_0, S_1)
  
  return(S_vec)
}

calculate_jacobian <- function(df, prm, psi_hat_vec){
  
  df_0 <- df %>% calculate_tau_rsd_m(prm=prm, psi_hat_vec=psi_hat_vec, m=0)
  
  jacobi_vec <- with(df_0, c(
    sum((a_0 - fit_a0)*tau_rsd*a_0),
    sum((a_0 - fit_a0)*tau_rsd*x_0*a_0),
    sum((a_0 - fit_a0)*tau_rsd*x_0*a_0),
    sum((a_0 - fit_a0)*tau_rsd*x_0*x_0*a_0)
  ))
  
  df_1 <- df %>% calculate_tau_rsd_m(prm=prm, psi_hat_vec=psi_hat_vec, m=1)
  
  jacobi_vec <- jacobi_vec + with(df_1, c(
    sum((a_0 - fit_a0)*tau_rsd*a_1 + 
          (a_1 - fit_a1)*tau_rsd*a_1),
    sum((a_0 - fit_a0)*tau_rsd*a_1*x_1 + 
          (a_1 - fit_a1)*tau_rsd*a_1*x_1),    
    sum((a_0 - fit_a0)*tau_rsd*a_1*x_0 + 
          (a_1 - fit_a1)*tau_rsd*a_1*x_1),
    sum((a_0 - fit_a0)*tau_rsd*a_1*x_0*x_1 + 
          (a_1 - fit_a1)*tau_rsd*a_1*x_1*x_1)
  ))
  
  return(matrix(jacobi_vec, byrow=T, nrow=2))
  #return(matrix(jacobi_vec, byrow=T, nrow= (length(psi_hat_CT) + length(psi_hat_TOT) - 1)))
}

calculate_variance <- function(df, prm, psi_hat_vec){
  
  df_1 <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k = 1)
  df_0 <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k = 0)
  
  n_0 <- prm$n
  n_1 <- dim(df_1)[[1]]
  
  D_vec <- c(rep(1, n_0), rep(0, n_1),
             df_0$x_0, rep(0, n_1),
             rep(0, n_0), rep(1, n_1),
             rep(0, n_0), df_1$a_0,
             rep(0, n_0), df_1$x_1)
  D <- matrix(D_vec, ncol=5, byrow=F)
  
  
  D_theta_vec <- c(rep(1, n_0 + n_1),
                   df_0$x_0, df_1$x_1)
  D_theta <- matrix(D_theta_vec, byrow=F, ncol=2)
  
  fit <- c(df_0$fit_a0, df_1$fit_a1)
  tau <- c(df_0$tau_k, df_1$tau_k)
  
  Jbb <- t(sapply(1:dim(D)[[2]], function(i){
    sapply(1:dim(D)[[2]], function(j){
      sum(D[,i]*D[,j]*fit*(1-fit))
    })
  }))
  
  Jtt <- t(sapply(1:dim(D_theta)[[2]], function(i){
    sapply(1:dim(D_theta)[[2]], function(j){
      sum(D_theta[,i]*D_theta[,j]*fit*(1-fit)*tau*tau)
    })
  }))
  
  Jtb <- t(sapply(1:dim(D_theta)[[2]], function(i){
    sapply(1:dim(D)[[2]], function(j){
      sum(D_theta[,i]*D[,j]*fit*(1-fit)*tau)
    })
  }))
  
  
  JTT <- Jtt - Jtb%*%solve(Jbb)%*%t(Jtb)
  
  Jacobian <- df_0 %>% calculate_jacobian(prm=prm, psi_hat_vec=psi_hat_vec)
  J_inv <- solve(Jacobian)
  
  (VCV <- J_inv %*% JTT %*% t(J_inv))
  
  return(J_inv %*% JTT %*% t(J_inv))
}


newton_raphson_grad <- function(df, prm,
                                psi_start_vec=c(0,0), 
                                tol=0.001, max_iter=50, psi_max = 10){
  
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
  
  if ((max(abs(psi_hat_vec)) > psi_max) | (steps == max_iter)) {
    fail_flag <- T
    psi_hat_vec <- c(NA, NA)
  }
  return(list(psi_hat_vec, steps, fail_flag))
}




nr_run <- function(psi_star_0=log(2), 
                   psi_star_1=log(1.5),
                   t_a=30,
                   t0_min = 10,
                   t0_max = 100,
                   n = 1000,
                   sims = 10){
  
  
  prm <- list(psi_star_0=psi_star_0,
              psi_star_1=psi_star_1,
              psi_star_vec=c(psi_star_0, psi_star_1),
              t_a=t_a,
              t0_min=t0_min,
              t0_max=t0_max,
              n=n,
              sims=sims)
  
  #prm$t_df <- get_time_df(prm=prm)
  
  cl <- makePSOCKcluster(28)
  registerDoParallel(cl)
  #registerDoParallel()
  
  nr_out <- foreach(icount(sims), .combine=rbind,
                    .export=c("create_sample", "fit_treatment_models",
                              "calculate_tau_k", "calculate_tau_rsd_m",
                              "calculate_score", "calculate_jacobian", 
                              "newton_raphson_grad",
                              "calculate_variance"
                    ), .packages="tidyverse") %dopar%
    {
      df <- create_sample(prm=prm)
      df <- df %>% fit_treatment_models(prm=prm)
      
      nri_out <- df %>% newton_raphson_grad(prm=prm)
      psi_hat_vec <- nri_out[[1]]
      
      (var_hat <- df %>% calculate_variance(psi_hat_vec=psi_hat_vec, prm=prm))
      
      c(psi_hat_vec, diag(var_hat))
      #c(unlist(psi_hat_list), diag(var_hat))
    }
  
  stopCluster(cl)
  
  for (j in 1:2){
    if (j == 1){ psi_lab <- "psi_0"}
    else {psi_lab <- "psi_1"}
    
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





nr_run(psi_star_0=log(2), 
       psi_star_1=log(1.5),
       t_a=30,
       t0_min = 10,
       t0_max = 100,
       n = 200,
       sims = 1000)
