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
  
  # a0_size <- max(floor(n/2),1)
  # a1_size <- max(floor(n/4),1)
  # x_size <-  max(floor(n/8),1)
  # df <- tibble(a_0 = c(rep(0, a0_size), rep(1, n - a0_size))) %>%
  #   add_column(a_1 = c(rep(0, a1_size), rep(1, a0_size - a1_size),
  #                      rep(0, a1_size), rep(1, n - (a0_size + a1_size)))) %>%
  #   add_column(x = c(rep(0, x_size), rep(1, a1_size - x_size),
  #                    rep(0, x_size), rep(1, (a0_size - a1_size) - x_size),
  #                    rep(0, x_size), rep(1, a1_size - x_size),
  #                    rep(0, x_size), rep(1, n - (a0_size + a1_size + x_size)))) %>%
  #   add_column(t0  = runif(n, min=t0_min, max=t0_max )) %>%
  #   mutate(ti = 0, t0_rsd = t0)
  
  a_vars <- length(prm$t_a_vec)
  k_parts <- 2^(a_vars + 1)
  
  n_s <- ceiling(prm$n_trgt/k_parts)
  prm$n <- n_s*k_parts
  
  df <- tibble(x = c(rep(0, n_s*(2^a_vars)), rep(1, n_s*(2^a_vars))))
  
  if(prm$censor){
    df <- df %>% add_column(C_i = prm$censor_date)
  }
  
  for (k in 0:(a_vars - 1)){
    df <- df %>% add_column(a_curr = rep(c(
      rep(0, n_s*(2^(a_vars-(k+1)))),
      rep(1, n_s*(2^(a_vars-(k+1))))
    ), 2^(k+1)))
    names(df)[names(df) == "a_curr"] <- paste0("a_", k)
  }
  #uniform
  #df <- df %>% add_column(t0  = runif(prm$n, min=prm$t0_min, max=prm$t0_max )) %>%
  #exponential
  df <- df %>% add_column(t0  = rexp(prm$n, rate=(1/prm$expmean))) %>%
    mutate(ti = 0, t0_rsd = t0)
  
  for(i in 1:length(prm$t_a_vec)){
    
    t_curr <- prm$t_a_vec[[i]]
    t_next <- ifelse(i < length(prm$t_a_vec), prm$t_a_vec[[i+1]] , Inf)
    beta_1_now <- ifelse(prm$beta_1_track[[i]]==0, 0, prm$psi_1_star[[prm$beta_1_track[[i]]]])
    beta_x_now <- ifelse(prm$beta_x_track[[i]]==0, 0, prm$psi_x_star[[prm$beta_x_track[[i]]]])
    
    names(df)[names(df)==paste0("a_", i-1)] <- "a_curr"
    
    df <- df %>% mutate(
      g_psi =  beta_1_now + x*beta_x_now,
      temp = pmin((t_next - t_curr)*exp(-g_psi*a_curr), t0_rsd),
      ti = ti + temp*exp(g_psi*a_curr),
      t0_rsd = t0_rsd - temp
    )
    
    names(df)[names(df)=="a_curr"] <- paste0("a_", i-1) 
  }
  
  if(prm$censor){
    df <- df %>% mutate(ti = pmin(ti, C_i))
  }
  
  df <- df %>% select(-c(t0_rsd, g_psi, temp))
  
  return(df)
}

fit_treatment_models <- function(df, prm){
  with(prm, {
    df <- df %>% mutate(fit_0 = glm(a_0 ~ 1 + x, family=binomial)$fitted.values)
    
    if(length(t_a_vec) > 1){
      for(i in 1:(length(t_a_vec)-1)){
        mdl_formula <- paste0("a_", i, "~ 1 + x +  ",
                              paste0("a_", 0:(i-1), collapse="+"))
        model_temp <- with(df %>% filter(.$ti > t_a_vec[[i+1]]), 
                           glm(  as.formula(mdl_formula),
                                 family=binomial))
        df$fit_temp <- predict(model_temp, newdata = df, type="response")
        names(df)[names(df) == "fit_temp"] <- paste0("fit_", i)
      }
    }
    
    return(df)
  })
}


calculate_tau_k <- function(df, psi_hat_vec, prm,  a_k=0, ...){
  #by default calculates tau(0) aka T0 from trial start
  
    t_a <- prm$t_a_vec[[a_k+1]]
    
    df <- df %>% filter(.$ti > t_a)%>% 
      mutate(ti_rsd = ti - t_a,
             tau_k = 0)
    
    if (prm$censor) {
      df <- df %>% mutate(
        C_rsd = C_i - t_a,
        C_k = 0
      )
    }
    
    for (i in (a_k+1):length(prm$t_a_vec)){
      
      t_curr <- prm$t_a_vec[[i]]
      t_next <- ifelse(i < length(prm$t_a_vec), prm$t_a_vec[[i+1]] , Inf)
      beta_1_now <- ifelse(prm$beta_1_track[[i]]==0, 0, 
                           psi_hat_vec[[prm$beta_1_track[[i]]]])
      beta_x_now <- ifelse(prm$beta_x_track[[i]]==0, 0, 
                           psi_hat_vec[[prm$beta_x_track[[i]]+length(prm$psi_1_star)]])
      
      names(df)[names(df)==paste0("a_", i-1)] <- "a_curr"
      
      df <- df %>% mutate(
        g_psi =  beta_1_now + x*beta_x_now,
        ti_temp = pmin(t_next - t_curr, ti_rsd),
        tau_k = tau_k + ti_temp*exp(-g_psi*a_curr),
        ti_rsd = ti_rsd - ti_temp
      )
      
      if (prm$censor) {
        t_next <- ifelse(i < length(prm$t_a_vec), prm$t_a_vec[[i+1]] , prm$censor_date)
        df <- df %>% mutate(
          C_psi = ifelse(g_psi < 0, 1, exp(-g_psi)),
          C_k = C_k + pmin(t_next - t_curr, C_rsd)*C_psi,
          C_rsd = C_rsd - pmin(t_next - t_curr, C_rsd)
        )
      }
      
      names(df)[names(df)=="a_curr"] <- paste0("a_", i-1) 
      
    }
    if (prm$censor) {
      df <- df %>% mutate(tau_k_un = tau_k,
                          delta_k = ifelse(tau_k < C_k, 0, 1),
                          tau_k = pmin(tau_k, C_k))
      return(select(df, -c(ti_temp, ti_rsd, g_psi, C_rsd, C_k, C_psi, tau_k_un)))
    } else {
      return(select(df, -c(ti_temp, ti_rsd, g_psi)))
    }
  }

calculate_tau_rsd_m <- function(df, prm, 
                                psi_hat_vec,
                                m=0, ...){
  with(prm, {
    
    t_curr <- t_a_vec[[m+1]]
    t_next <- ifelse(m+1 < length(t_a_vec), t_a_vec[[m+2]] , Inf)
    beta_1_now <- ifelse(beta_1_track[[m+1]]==0, 0,
                         psi_hat_vec[[beta_1_track[[m+1]]]])
    beta_x_now <- ifelse(beta_x_track[[m+1]]==0, 0,
                         psi_hat_vec[[beta_x_track[[m+1]]+length(psi_1_star)]])
    
    names(df)[names(df)==paste0("a_", m)] <- "a_curr"
    
    df <- df %>% mutate(
      g_psi =  beta_1_now + x*beta_x_now,
      ti_temp = pmax(0, pmin(t_next, ti) - t_curr),
      tau_rsd = ti_temp*exp(-g_psi*a_curr)
    )
    
    names(df)[names(df)=="a_curr"] <- paste0("a_", m)
    
    return(df %>% select(-c(g_psi, ti_temp)))
  })
  
}


calculate_C_m <- function(psi_hat_vec, prm, m=0){
  #,a_k=0, j=1){
  C_m <- 0
  
  t_now <- prm$t_a_vec[[m+1]]
  t_next <- ifelse(m+1 < length(prm$t_a_vec), prm$t_a_vec[[m+2]] , prm$censor_date)
  beta_1_now <- ifelse(prm$beta_1_track[[m+1]]==0, 0,
                       psi_hat_vec[[prm$beta_1_track[[m+1]]]])
  #beta_x_now <- ifelse(beta_x_track[[m+1]]==0, 0,
  #                     psi_hat_vec[[beta_x_track[[m+1]]+length(psi_1_star)]])
  if (beta_1_now >= 0){
    C_m <- C_m + exp(-beta_1_now)*(t_next - t_now)
  }
  
  return(C_m)
}

calculate_C_k <- function(psi_hat_vec, prm, a_k=0, j){
  #with(prm, {
  #C_k <- 0
  #for (m in (a_k:(length(t_a_vec)-1))) {
  #  if(psi_hat_vec[[j]] >= 0){
  #    t_now <- t_a_vec[[m+1]]
  #    t_next <- ifelse(m == (length(t_a_vec)-1), #then
  #                     prm$censor_date, #else
  #                     t_a_vec[[m+2]])
  #    C_k <- C_k + exp(-psi_hat_vec[[j]])*(t_next - t_now)
  #  }
  #}
  
  #return(C_k)
  #})
  C_k <- 0
  for (m in (a_k:(length(prm$t_a_vec)-1))) {
    C_k <- C_k + calculate_C_m(psi_hat_vec=psi_hat_vec, prm=prm, m=m)
  }
  return(C_k)
  
}



calculate_score <- function(df, prm, psi_hat_vec){
  
  S_vec <- c()
  for (i in 1:length(prm$psi_1_star)){
    S_curr <- 0
    
    for (a_curr in (0:(length(prm$t_a_vec)-1))){
      if (prm$beta_1_track[[a_curr+1]] == i) {
        
        df_temp <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=a_curr)
        
        names(df_temp)[names(df_temp)==paste0("a_", a_curr)] <- "a_k"
        names(df_temp)[names(df_temp)==paste0("fit_", a_curr)] <- "fit_k"
        
        S_curr <- S_curr + with(df_temp, sum((a_k - fit_k)*tau_k))
      }
    }
    
    S_vec <- c(S_vec, S_curr)
    
  }
  
  if(length(prm$psi_x_star) > 0) {
    for (i in 1:length(prm$psi_x_star)){
      S_curr <- 0
      
      for (a_curr in (0:(length(prm$t_a_vec)-1))){
        if (prm$beta_x_track[[a_curr+1]] == i) {
          
          df_temp <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=a_curr)
          
          names(df_temp)[names(df_temp)==paste0("a_", a_curr)] <- "a_k"
          names(df_temp)[names(df_temp)==paste0("fit_", a_curr)] <- "fit_k"
          
          S_curr <- S_curr + with(df_temp, sum((a_k - fit_k)*x*tau_k))
        }
      }
      
      S_vec <- c(S_vec, S_curr)
      
    }
  }
  return(S_vec)
  
}

calculate_jacobian <- function(df, prm, psi_hat_vec){
  jacobi_vec <- c()
  jacobi_dim <- length(psi_hat_vec)
  for (row_counter in 1:jacobi_dim){
    for(col_counter in 1:jacobi_dim) {
      
      if (row_counter <= length(prm$psi_1_star)) {
        d_num <- "theta_1"
        i <- row_counter
        beta_track_row <- prm$beta_1_track
      } else {
        d_num <- "theta_x"
        i <- row_counter - length(prm$psi_1_star)
        beta_track_row <- prm$beta_x_track
      }
      
      if (col_counter <= length(prm$psi_1_star)) {
        d_wrt <- "theta_1"
        j <- col_counter
        beta_track_col <- prm$beta_1_track
      } else {
        d_wrt <- "theta_x"
        j <- col_counter - length(prm$psi_1_star)
        beta_track_col <- prm$beta_x_track
      }
      
      denom <- 0
      
      for (a_curr in (0:(length(prm$t_a_vec)-1))){
        if (beta_track_row[[a_curr+1]] == i) {
          
          df_k <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=a_curr) %>%
            mutate(Si_temp = 0)
          
          dC_dpsi <- 0
          for (m in (a_curr:(length(prm$t_a_vec)-1))) {
            if ((beta_track_row[[m+1]] == i) &&
                (beta_track_col[[m+1]] == j)) {
              
              df_k <- df_k %>% calculate_tau_rsd_m(prm=prm, psi_hat_vec=psi_hat_vec, m=m)
              
              names(df_k)[names(df_k)==paste0("a_", m)] <- "a_m"
              
              df_k <- df_k %>% mutate(
                x_wrt = (1 - as.integer(d_wrt == "theta_x")) + x * as.integer(d_wrt == "theta_x"),
                Si_temp = Si_temp + a_m * x_wrt * tau_rsd)
              
              names(df_k)[names(df_k)=="a_m"] <- paste0("a_", m)
              
              if(prm$censor){
                C_m <- calculate_C_m(psi_hat_vec=psi_hat_vec, prm=prm, m=m)
                dC_dpsi <- dC_dpsi + (psi_hat_vec[[j]] >= 0)*C_m
              }
            }
          }
          
          df_k <- df_k %>% mutate(
            x_num = (1 - as.integer(d_num == "theta_x")) + x * as.integer(d_num == "theta_x"),

            #C_k = calculate_C_k(psi_hat_vec=psi_hat_vec, prm=prm, a_k=a_curr, j=j)
          )
          if(prm$censor) { 
            df_k <- df_k %>% mutate(
                      Si_censor = delta_k*dC_dpsi + (1-delta_k)*Si_temp,
                      C_k = calculate_C_k(psi_hat_vec=psi_hat_vec, prm=prm, a_k=a_curr)
            )
            #df_k <- df_k %>% mutate(Si_censor = delta_k*(psi_hat_vec[[j]] >= 0)*C_k + 
            #                          (1-delta_k)*Si_temp)
          } else {
            df_k <- df_k %>% mutate(
              Si_censor = Si_temp)
          }
          
          names(df_k)[names(df_k)==paste0("a_", a_curr)] <- "a_k"
          names(df_k)[names(df_k)==paste0("fit_", a_curr)] <- "fit_k"
          
          denom <- denom + sum(with(df_k, (a_k - fit_k)*x_num*Si_censor))
        }
      }
      
      jacobi_vec <- c(jacobi_vec, denom)
    }
  }
  
  return(matrix(jacobi_vec, byrow=T, nrow=jacobi_dim))
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
    psi_hat_vec <- rep(NA, length(prm$psi_star_vec))
  }
  return(list(psi_hat_vec, steps, fail_flag))
}


calculate_variance <- function(df, prm, psi_hat_vec){
  
  n_vec <- c()
  for (t_a in prm$t_a_vec){
    n_vec <- c(n_vec, nrow(df %>% filter(.$ti > t_a)))
  }
  
  D_vec <- c()
  col_count <- 0
  for (k in 0:(length(prm$t_a_vec)-1)) {
    df_temp <- df %>% filter(.$ti > prm$t_a_vec[[k+1]])
    
    for (sub_col in -2:(k-1)) {
      
      for (section in 0:(length(prm$t_a_vec)-1)) {
        if (section == k) {
          if (sub_col == -2) {
            D_vec <- c(D_vec, rep(1, n_vec[[section+1]]))
          } else if (sub_col == -1) {
            D_vec <- c(D_vec, df_temp$x)
          } else {
            D_vec <- c(D_vec, df_temp[[paste0('a_', sub_col)]]) 
          }
        } else {
          D_vec <- c(D_vec, rep(0, n_vec[[section+1]]))
        }
      }
      
      col_count <- col_count + 1
    }
    
  }
  D <- matrix(D_vec, ncol=col_count, byrow=F)
  
  
  D_theta_vec <- c()
  for (beta_1_val in 1:max(prm$beta_1_track)){
    for (k in 0:(length(prm$t_a_vec)-1)){
      if (prm$beta_1_track[[k+1]] == beta_1_val) {
        D_theta_vec <- c(D_theta_vec, rep(1, n_vec[[k+1]]))
      } else {
        D_theta_vec <- c(D_theta_vec, rep(0, n_vec[[k+1]]))
      }
    }  
  }
  if (max(prm$beta_x_track) > 0) {
    for (beta_x_val in 1:max(prm$beta_x_track)){
      for (k in 0:(length(prm$t_a_vec)-1)){
        if (prm$beta_x_track[[k+1]] == beta_x_val) {
          df_temp <- df %>% filter(.$ti > prm$t_a_vec[[k+1]])
          D_theta_vec <- c(D_theta_vec, df_temp$x)
        } else {
          D_theta_vec <- c(D_theta_vec, rep(0, n_vec[[k+1]]))
        }
      }  
    }
  }
  D_theta <- matrix(D_theta_vec, ncol=max(prm$beta_1_track)+max(prm$beta_x_track),
                    byrow=F)
  
  fit <- c()
  tau <- c()
  for (k in 0:(length(prm$t_a_vec)-1)){
    df_temp <- df %>% filter( .$ti > prm$t_a_vec[[k+1]]) %>% 
      calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=k)
    fit <- c(fit, df_temp[[paste0('fit_', k)]])
    tau <- c(tau, df_temp$tau_k)
  }
  
  
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
  
  Jacobian <- df %>% calculate_jacobian(prm=prm, psi_hat_vec=psi_hat_vec)
  J_inv <- solve(Jacobian)
  
  (VCV <- J_inv %*% JTT %*% t(J_inv))
  
  return(VCV)
}


nr_run <- function(prm){
  
  cl <- makePSOCKcluster(28)
  registerDoParallel(cl)
  
  nr_out <- foreach(seed_curr=1:prm$sims, .combine=rbind, .errorhandling='remove',
                    .export=c(
                      "create_sample", "fit_treatment_models",
                      "calculate_tau_k", "calculate_tau_rsd_m",
                      "calculate_score", "calculate_jacobian", 
                      "calculate_C_k", "calculate_C_m",
                      "newton_raphson_grad",
                      "calculate_variance"
                    ), .packages="tidyverse") %dopar%
    {
      
      nr_out_single <- rep(NA, 2*length(prm$psi_star_vec))
      
      set.seed(seed_curr)
      df <- create_sample(prm=prm)
      df <- df %>% fit_treatment_models(prm=prm)
      
      nri_out <- df %>% newton_raphson_grad(prm=prm, 
                                            psi_start_vec=rep(0, length(prm$psi_star_vec)))
      (psi_hat_vec <- nri_out[[1]])
      
      (var_hat <- df %>% calculate_variance(psi_hat_vec=psi_hat_vec, prm=prm))
      
      if(max(diag(var_hat)) < 1) {
        nr_out_single <- c(psi_hat_vec, diag(var_hat))
      }
      
      nr_out_single
    }
  
  stopCluster(cl)
  results_vec <- calculate_results(nr_out, prm)
  
  
  return(list(nr_out=nr_out, results_vec=results_vec))
}

calculate_results <- function(nr_out, prm){
  results_vec <- c()
  for (j in 1:length(prm$psi_star_vec)){
    results_vec <- c(results_vec, 
                     mean(nr_out[,j] , na.rm=T),
                     var(nr_out[,j], na.rm=T),
                     mean(nr_out[,j+length(prm$psi_star_vec)], na.rm=T),
                     var(nr_out[,j+length(prm$psi_star_vec)], na.rm=T)
    ) 
  }
  
  return(results_vec)
}

print_results <- function(results_vec, prm){
  cat(#"---------------------------", "\n",
    prm$sim_label, "\n",
    "(n = ", prm$n, ", censoring=", prm$censor, ") \n", 
    #  "---------------------------", "\n",
    sep="")
  for (j in 1:length(prm$psi_star_vec)){
    cat(
      paste(prm$psi_lab[[j]], "\t", "MEAN", "\t\t", "VAR", collapse="\t"),
      "\n ",
      
      paste(
        c("mu    ", 
          sprintf("%.6f", c(prm$psi_star_vec[[j]]))
        ),
        collapse = "\t"),
      "\n ",
      
      paste(
        c("xbar  ",
          sprintf("%.6f",
                  c(results_vec[4*(j-1)+1],
                    results_vec[4*(j-1)+2]
                  ))
        ),
        collapse = "\t"),
      "\n ",
      
      paste(
        c("aVar  ",
          sprintf("%.6f",
                  c(results_vec[4*(j-1)+3],
                    results_vec[4*(j-1)+4]
                  ))
        ),
        collapse = "\t"),
      "\n", sep="")
  }
  cat("\n\n")
  
}





two_factor_variance_analysis <- function(sims=3000, 
                                         n_trgt_vec = c(100, 400, 800, 2000)){
  
  prm <- list() 
  
  prm$sims <- sims
  
  prm$sim_label <- "(const_1) -> (const_2)"
  prm$t_a_vec <- c(0, 35)
  prm$expmean <- 50
  prm$beta_1_track <- c(1, 2)
  prm$beta_x_track <- c(0, 0)
  prm$psi_lab <- c("psi_1", "psi_2")
  prm$psi_1_star <- c(log(2), log(2))
  prm$psi_x_star <- c()
  prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
  prm$censor <- F
  
  var_df <- tibble(n_trgt=c(), sample_var=c(), asymp_var=c(), var_ratio=c())
  
  variance_vec <- c()
  sample_var <- c()
  asymp_var <- c()
  var_ratio <- c()
  for(n_trgt in n_trgt_vec){
    print(n_trgt)
    prm$n_trgt <- n_trgt
    
    nr_out_list <- nr_run(prm=prm)
    results_vec <- nr_out_list$results_vec
    
    var_df <- var_df %>% add_row(n_trgt=prm$n_trgt,
                                 sample_var=results_vec[2],
                                 asymp_var=results_vec[3],
                                 var_ratio=results_vec[3]/results_vec[2]
    )
    
    sample_var <- c(sample_var, results_vec[2])
    asymp_var <- c(asymp_var, results_vec[3])
    var_ratio <- c(var_ratio, results_vec[3]/results_vec[2])
  }
  
  #mutate(var_type = replace(var_type, var_type=='asymp_var', 'Mean asymptotic variance'),
  #       var_type = replace(var_type, var_type=='asymp_var', sample_var='Sample variance of the sample means')) %>%
  
  var_df %>% select(-c(var_ratio)) %>%
    rename('Sample variance of the sample means'=sample_var,
           'Mean asymptotic variance'=asymp_var) %>%
    pivot_longer(cols = !n_trgt, names_to = "var_type", values_to = "variance") %>% 
    ggplot(aes(x=n_trgt, y=variance, colour=var_type, line_style=var_type)) + 
    geom_line(aes(linetype=var_type)) + labs(x='n') + 
    theme(legend.title = element_blank())
  
  #var_df %>% pivot_longer(cols = !n_trgt, names_to = "var_type", values_to = "variance") %>%
  #  ggplot(aes(x=n_trgt)) + 
  #  geom_line(aes(y=variance), group_by = var_type)
    #geom_line(aes(y=asymp_var), color = "steelblue", linetype="dashed") +
    #geom_line(aes(y=sample_var), color = "darkred")

  
  return(var_df)
}




prm <- list()
#prm <- list(t_a_vec = c(0, 30),
#            
#            t0_min = 10,
#            t0_max = 100,
#            
#            n_trgt = 400
#)

##################################################################
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           (1) -> (1) block                            #   ##
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
##################################################################


prm$sim_label <- "(const_1) -> (const_1)"
prm$t_a_vec <- c(0, 35)
prm$expmean <- 50
prm$n_trgt <- 400
prm$beta_1_track <- c(1, 1)
prm$beta_x_track <- c(0, 0)
prm$psi_lab <- c("psi_1")
prm$psi_1_star <- c(log(2))
prm$psi_x_star <- c()
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$sims <- 3000
prm$censor <- F
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)
prm$censor <- T
prm$censor_date <- 80
nr_out <- nr_run(prm=prm)

prm$sim_label <- "(const_1) -> (const_1)"
prm$t_a_vec <- c(0, 35)
prm$expmean <- 50
prm$n_trgt <- 400
prm$beta_1_track <- c(1, 1)
prm$beta_x_track <- c(0, 0)
prm$psi_lab <- c("psi_1")
prm$psi_1_star <- c(log(1))
prm$psi_x_star <- c()
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$sims <- 3000
prm$censor <- F
nr_out <- nr_run(prm=prm)
prm$censor <- T
prm$censor_date <- 80
nr_out <- nr_run(prm=prm)

prm$sim_label <- "(const_1) -> (const_1)"
prm$t_a_vec <- c(0, 35)
prm$expmean <- 50
prm$n_trgt <- 400
prm$beta_1_track <- c(1, 1)
prm$beta_x_track <- c(0, 0)
prm$psi_lab <- c("psi_1")
prm$psi_1_star <- c(log(1/2))
prm$psi_x_star <- c()
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$sims <- 3000
prm$censor <- F
nr_out <- nr_run(prm=prm)
prm$censor <- T
prm$censor_date <- 80
nr_out <- nr_run(prm=prm)


##################################################################
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           (1) -> (2) block                            #   ##
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
##################################################################


prm$sim_label <- "(const_1) -> (const_2)"
prm$t_a_vec <- c(0, 35)
prm$expmean <- 50
prm$n_trgt <- 800
prm$beta_1_track <- c(1, 2)
prm$beta_x_track <- c(0, 0)
prm$psi_lab <- c("psi_1", "psi_2")
prm$psi_1_star <- c(log(2), log(2))
prm$psi_x_star <- c()
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$sims <- 5000
prm$censor <- F
nr_out <- nr_run(prm=prm)
prm$censor <- T
prm$censor_date <- 80
nr_out <- nr_run(prm=prm)




prm$sim_label <- "(const_1) -> (const_1)"
prm$t_a_vec <- c(0, 35)
prm$expmean <- 50
prm$n_trgt <- 400
prm$beta_1_track <- c(1, 1)
prm$beta_x_track <- c(0, 0)
prm$psi_lab <- c("psi_1")
prm$psi_1_star <- c(log(1))
prm$psi_x_star <- c()
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$sims <- 3000
prm$censor <- F
nr_out <- nr_run(prm=prm)
prm$censor <- T
prm$censor_date <- 80
nr_out <- nr_run(prm=prm)

prm$sim_label <- "(const_1) -> (const_1)"
prm$t_a_vec <- c(0, 35)
prm$expmean <- 50
prm$n_trgt <- 400
prm$beta_1_track <- c(1, 1)
prm$beta_x_track <- c(0, 0)
prm$psi_lab <- c("psi_1")
prm$psi_1_star <- c(log(1/2))
prm$psi_x_star <- c()
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$sims <- 3000
prm$censor <- F
nr_out <- nr_run(prm=prm)
prm$censor <- T
prm$censor_date <- 80
nr_out <- nr_run(prm=prm)














