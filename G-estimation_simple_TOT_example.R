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

  
  a_0 <- c(rep(0, floor(prm$n/2)), rep(1, prm$n - floor(prm$n/2)))
  a_1 <- sample(a_0, size=prm$n, replace=F)
  
  df <- tibble(a_0, a_1) %>% 
        add_column(t0  = runif(prm$n, min=prm$t0_min, max=prm$t0_max )) %>%
        mutate(ti = 0, t0_rsd = t0,
              g_psi = prm$psi_star_vec[[1]],
              temp = pmin((prm$t_a - 0)*exp(-g_psi*a_0), t0_rsd),
              ti = ti + temp*exp(g_psi*a_0),
              t0_rsd = t0_rsd - temp,
              
              g_psi = prm$psi_star_vec[[1]] + a_0*prm$psi_star_vec[[2]],
              temp = t0_rsd,
              ti = ti + temp*exp(g_psi*a_1),
              t0_rsd = t0_rsd - temp
        ) %>% select(-c(t0_rsd, g_psi, temp))
  
  return(df)
  #return(select(df, -c(t0_rsd, temp, w_TOT, g_psi)))
}

fit_treatment_models <- function(df, prm){
  
  df <- df %>% mutate(fit_a0 = glm(a_0 ~ 1, family=binomial)$fitted.values)
  
  mdl_formula <- "a_1 ~ a_0"
  model_temp <- with(df %>% filter(.$ti > prm$t_a), 
                     glm(  as.formula(mdl_formula),
                           family=binomial))
  df$fit_a1 <- predict(model_temp, newdata = df, type="response")
  
  return(df)
}





calculate_tau_k <- function(df, prm, psi_hat_vec,  a_k=0, ...){
  
  df <- df %>% mutate(ti_temp = pmax(ti - prm$t_a, 0),
                      tau_k = exp(-(psi_hat_vec[[1]] + psi_hat_vec[[2]]*a_0)*a_1)*ti_temp)
  if (a_k == 0) {
    df <- df %>% mutate(ti_temp = pmax(pmin(prm$t_a, ti), 0),
                        tau_k = tau_k + exp(-(psi_hat_vec[[1]])*a_0)*ti_temp)
  } else {
    #a_k == 1
    df <- df %>% filter(.$ti > prm$t_a)
  }
  
  df <- df %>% select(-c(ti_temp))
  
  return(df)
}

calculate_tau_rsd_m <- function(df, prm, 
                                psi_hat_vec,
                                m=0, ...){
  
  if (m == 0) {
    df <- df %>% mutate(g_psi = psi_hat_vec[[1]],
                        ti_temp = pmax(0, pmin(prm$t_a, ti) - 0),
                        tau_rsd = ti_temp*exp(-g_psi*a_0)
    )
  } else {
    df <- df %>% mutate(g_psi = psi_hat_vec[[1]] + psi_hat_vec[[2]]*a_0,
                        ti_temp = pmax(0, ti - prm$t_a),
                        tau_rsd = ti_temp*exp(-g_psi*a_1) 
    )
  }
  
  df <- df %>% select(-c(g_psi, ti_temp))
  return(df )
}

calculate_score <- function(df, prm, psi_hat_vec){
  
  df_temp <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=1)
  S_0 <- with(df_temp, sum((a_1 - fit_a1)*tau_k))
  S_1 <- with(df_temp, sum((a_1 - fit_a1)*tau_k*a_0))
  
  df_temp <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=0)
  S_0 <- S_0 + with(df_temp, sum((a_0 - fit_a0)*tau_k))
  
  S_vec <- c(S_0, S_1)
  
  return(S_vec)
}

calculate_jacobian <- function(df, prm, psi_hat_vec){
  
  df_temp <- df %>% calculate_tau_rsd_m(prm=prm, psi_hat_vec=psi_hat_vec, m=0)
  
  jacobi_vec <- c(
    with(df_temp, sum((a_0 - fit_a0)*tau_rsd*a_0)),
    0,
    0,
    0
  )
  
  df_temp <- df %>% calculate_tau_rsd_m(prm=prm, psi_hat_vec=psi_hat_vec, m=1)
  
  jacobi_vec <- jacobi_vec + c(
    with(df_temp, sum((a_0 - fit_a0)*tau_rsd*a_1 + 
                      (a_1 - fit_a1)*tau_rsd*a_1)),
    with(df_temp, sum((a_1 - fit_a1)*tau_rsd*a_0*a_1)),
    with(df_temp, sum((a_1 - fit_a1)*tau_rsd*a_0*a_1)),
    with(df_temp, sum((a_1 - fit_a1)*tau_rsd*a_0*a_1))
  )
  
  
  
  return(matrix(jacobi_vec, byrow=T, nrow=2))
  #return(matrix(jacobi_vec, byrow=T, nrow= (length(psi_hat_CT) + length(psi_hat_TOT) - 1)))
}

calculate_variance <- function(df, prm, psi_hat_vec){
  
  df_1 <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k = 1)
  df_0 <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k = 0)
  
  
  
  n_0 <- prm$n
  n_1 <- dim(df_1)[[1]]
  
  D_vec <- c(rep(1, n_0), rep(0, n_1),
             rep(0, n_0), rep(1, n_1),
             rep(0, n_0), df_1$a_0)
  D <- matrix(D_vec, ncol=3, byrow=F)
  
  
  D_theta_vec <- c(rep(1, n_0 + n_1),
                   rep(0, n_0), df_1$a_0)
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
                   psi_star_a0=log(1.5),
                   t_a=30,
                   t0_min = 10,
                   t0_max = 100,
                   
                   n = 1000,
                   sims = 10){
  
  
  prm <- list(psi_star_0=psi_star_0,
              psi_star_a0=psi_star_a0,
              psi_star_vec=c(psi_star_0, psi_star_a0),
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
       psi_star_a0=log(1.5),
       t_a=30,
       t0_min = 10,
       t0_max = 100,
       n = 100,
       sims = 1000)











draw_heatmap(psi_star_0=log(2), 
             psi_star_a0=log(1.5),
             t_a=30,
             t0_min = 10,
             t0_max = 100,
             n = 200,
             psi_0_vals = seq(-1, 2, 0.1),
             psi_1_vals = seq(-1, 2, 0.1)){
  
  prm <- list(psi_star_0=psi_star_0,
              psi_star_a0=psi_star_a0,
              psi_star_vec=c(psi_star_0, psi_star_a0),
              t_a=t_a,
              t0_min=t0_min,
              t0_max=t0_max,
              n=n)
  
  
  df <- create_sample(prm=prm)
  df <- df %>% fit_treatment_models(prm=prm)
  
  plt.df <- tibble(psi_0=numeric(), psi_1=numeric(),
                   theta_0=numeric(), theta_1=numeric())
  
  for (psi_0 in psi_0_vals) {
    for (psi_1 in psi_1_vals) {
      score <- df %>% calculate_score(prm=prm, psi_hat_vec=c(psi_0, psi_1))
      
      plt.df <- plt.df %>% add_row(psi_0=psi_0, psi_1=psi_1,
                                   theta_0=score[[1]], theta_1=score[[2]])
    }
  }
  
  plt.df <- plt.df %>% mutate(theta_0 = pmin(theta_0, 100),
                              theta_1 = pmin(theta_1, 100))
  
}

ggplot(plt.df, aes(psi_0, psi_1, fill=theta_0)) + geom_tile() +
  #scale_fill_distiller(limits=c(-0.1,0.1), palette="OrRd") 
  scale_fill_gradientn(limits=c(-100,100),
                       colours=c("navyblue", "white", "darkorange1"))

