library(matlib)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(foreach)
library(doParallel)
library(tidyverse)
library(gt)
library(labelled)
library(DescTools)




fit_treatment_models <- function(df, prm){
  
  #trt_models <- list()
  trt_models <- lapply(0:(length(prm$trt_mod_list)-1),
      function(k){      
        mdl_formula <- paste0("a_", k, " ~ ", 
                              paste0(prm$trt_mod_list[[k+1]], collapse = " + "))
        
        #tryCatch({
        df_k <- df %>% filter(.$ti > prm$t_a_vec[[k+1]])
        model_temp <- with(df_k, 
                           glm(  as.formula(mdl_formula),
                                 family=binomial))
        
        #df_k$fit_temp <- predict(model_temp, newdata = df_k, type="response")
        #df
        #names(df)[names(df) == "fit_temp"] <- paste0("fit_", k)
          
          #return(model_temp)
        #}, warning = function(w) {cat("warning at k = ", k)})
      })
  
  
  for(k in 0:(length(trt_models)-1)){
    df_k <- df %>% filter(.$ti > prm$t_a_vec[[k+1]])
    model_temp <- trt_models[[k+1]]
    df_k$fit_temp <- predict(model_temp, newdata = df_k, type="response")
    df_k <- df_k %>% select(c(id, fit_temp))
    df <- df %>% left_join(df_k, by="id")
    names(df)[names(df) == "fit_temp"] <- paste0("fit_", k)
  }
  
  return(list(df=df, trt_models=trt_models))
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
                         psi_hat_vec[[prm$beta_x_track[[i]]+max(prm$beta_1_track)]])
    
    names(df)[names(df)==paste0("a_", i-1)] <- "a_curr"
    
    df <- df %>% mutate(
      g_psi =  beta_1_now + x*beta_x_now,
      ti_temp = pmin(t_next - t_curr, ti_rsd),
      tau_k = ifelse(is.na(a_curr), tau_k, tau_k + ti_temp*exp(-g_psi*a_curr)),
      ti_rsd = ti_rsd - ti_temp,
      
    )
    
    if (prm$censor) {
      t_next <- ifelse(i < length(prm$t_a_vec), 
                       prm$t_a_vec[[i+1]], 
                       prm$censor_date)
      df <- df %>% mutate(
        C_psi = ifelse(g_psi < 0, 1, exp(-g_psi)),
        C_k = C_k + pmin(t_next - t_curr, C_rsd)*C_psi,
        C_rsd = C_rsd - pmin(t_next - t_curr, C_rsd)
      )
    }
    
    names(df)[names(df)=="a_curr"] <- paste0("a_", i-1)
    
    
    #cat("i=",i,"\n na=", sum(is.na(df$tau_k)))
    
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

calculate_tau_rsd_m <- function(df, psi_hat_vec, prm,
                                m=0, ...){
  
    t_curr <- prm$t_a_vec[[m+1]]
    t_next <- ifelse(m+1 < length(prm$t_a_vec), prm$t_a_vec[[m+2]] , Inf)
    beta_1_now <- ifelse(prm$beta_1_track[[m+1]]==0, 0,
                         psi_hat_vec[[prm$beta_1_track[[m+1]]]])
    beta_x_now <- ifelse(prm$beta_x_track[[m+1]]==0, 0,
                         psi_hat_vec[[prm$beta_x_track[[m+1]]+max(prm$beta_1_track)]])
    
    names(df)[names(df)==paste0("a_", m)] <- "a_curr"
    
    df <- df %>% mutate(
      g_psi =  beta_1_now + x*beta_x_now,
      ti_temp = pmax(0, pmin(t_next, ti) - t_curr),
      tau_rsd =  ifelse(is.na(a_curr), 0, ti_temp*exp(-g_psi*a_curr))
    )
    
    names(df)[names(df)=="a_curr"] <- paste0("a_", m)
    
    return(df %>% select(-c(g_psi, ti_temp)))
}

calculate_C_m <- function(psi_hat_vec, prm, m=0){
  #,a_k=0, j=1){
  
  
  t_now <- prm$t_a_vec[[m+1]]
  t_next <- ifelse(m+1 < length(prm$t_a_vec), 
                   prm$t_a_vec[[m+2]] , 
                   prm$censor_date)
  
  beta_1_now <- ifelse(prm$beta_1_track[[m+1]]==0, 0,
                       psi_hat_vec[[prm$beta_1_track[[m+1]]]])
  beta_x_now <- ifelse(prm$beta_x_track[[m+1]]==0, 0,
                       psi_hat_vec[[prm$beta_x_track[[m+1]]+max(prm$beta_1_track)]])
  
  C_m <- exp(-(beta_1_now*(beta_1_now>=0) + beta_x_now*(beta_x_now>=0)))*(t_next - t_now)
  
  #beta_x_now <- ifelse(beta_x_track[[m+1]]==0, 0,
  #                     psi_hat_vec[[beta_x_track[[m+1]]+length(psi_1_star)]])
  #if (beta_1_now >= 0){
  
  #}
  
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

calculate_score <- function(df, psi_hat_vec, prm){
  
  S_vec <- c()
  for (i in 1:max(prm$beta_1_track)){
    
    S_curr <- 0
    
    for (a_curr in (0:(length(prm$t_a_vec)-1))){
      if (prm$beta_1_track[[a_curr+1]] == i) {
        
        df_temp <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=a_curr)
        
        names(df_temp)[names(df_temp)==paste0("a_", a_curr)] <- "a_k"
        names(df_temp)[names(df_temp)==paste0("fit_", a_curr)] <- "fit_k"
        
        df_temp <- df_temp %>% mutate(
                    S_inc = ifelse((is.na(a_k) | is.na(fit_k) | is.na(tau_k)), 0,
                                                     (a_k - fit_k)*tau_k))
        
        #df_temp <- df %>% calculate_tau_rsd_m(., psi_hat_vec, prm)
        #names(df_temp)[names(df_temp)==paste0("a_", a_curr)] <- "a_k"
        #names(df_temp)[names(df_temp)==paste0("fit_", a_curr)] <- "fit_k"
        
        #df_temp <- df_temp %>% mutate(S_inc = ifelse(is.na(a_k), 0,
        #                                             (a_k - fit_k)*tau_rsd))
        S_curr <- S_curr + with(df_temp, sum(S_inc))
      }
    }
    
    S_vec <- c(S_vec, S_curr)
    
  }
  
  if(max(prm$beta_x_track) > 0) {
    for (i in 1:max(prm$beta_x_track)){
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
  
  df_k_list <- lapply(0:(length(prm$t_a_vec)-1), 
                      FUN=function(k){
                        return(df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=k) %>%
                          mutate(Si_temp = 0))
                      })
  
  #for (row_counter in 1:jacobi_dim){
  jacobian <- sapply(1:jacobi_dim, 
  function(row_counter){
    
    #for(col_counter in 1:jacobi_dim) {
    sapply(1:jacobi_dim,
    function(col_counter){         

      
      if (row_counter <= max(prm$beta_1_track)) {
        d_num <- "theta_1"
        i <- row_counter
        beta_track_row <- prm$beta_1_track
      } else {
        d_num <- "theta_x"
        i <- row_counter - max(prm$beta_1_track)
        beta_track_row <- prm$beta_x_track
      }
      
      if (col_counter <= max(prm$beta_1_track)) {
        d_wrt <- "theta_1"
        j <- col_counter
        beta_track_col <- prm$beta_1_track
      } else {
        d_wrt <- "theta_x"
        j <- col_counter - max(prm$beta_1_track)
        beta_track_col <- prm$beta_x_track
      }
      
      denom <- 0
      
      for (a_curr in (0:(length(prm$t_a_vec)-1))){
        if (beta_track_row[[a_curr+1]] == i) {
          
          df_k <- df_k_list[[a_curr+1]]
          #df_k <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=a_curr) %>%
          #  mutate(Si_temp = 0)
          
          dC_dpsi <- 0
          for (m in (a_curr:(length(prm$t_a_vec)-1))) {
            if ((beta_track_row[[m+1]] == i) &&
                (beta_track_col[[m+1]] == j)) {
              
              df_k <- df_k %>% calculate_tau_rsd_m(prm=prm, psi_hat_vec=psi_hat_vec, m=m)
              
              names(df_k)[names(df_k)==paste0("a_", m)] <- "a_m"
              
              df_k <- df_k %>% mutate(
                x_wrt = (1 - as.integer(d_wrt == "theta_x")) + x * as.integer(d_wrt == "theta_x"),
                Si_inc = ifelse(is.na(a_m), 0, a_m * x_wrt * tau_rsd),
                Si_temp = Si_temp + Si_inc)
              
              names(df_k)[names(df_k)=="a_m"] <- paste0("a_", m)
              
              if(prm$censor){
                C_m <- calculate_C_m(psi_hat_vec=psi_hat_vec, prm=prm, m=m)
                dC_dpsi <- dC_dpsi + as.integer(psi_hat_vec[[j]] > 0)*C_m
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
              #C_k = calculate_C_k(psi_hat_vec=psi_hat_vec, prm=prm, a_k=a_curr)
            )
            #df_k <- df_k %>% mutate(Si_censor = delta_k*(psi_hat_vec[[j]] >= 0)*C_k + 
            #                          (1-delta_k)*Si_temp)
          } else {
            df_k <- df_k %>% mutate(
              Si_censor = Si_temp)
          }
          
          names(df_k)[names(df_k)==paste0("a_", a_curr)] <- "a_k"
          names(df_k)[names(df_k)==paste0("fit_", a_curr)] <- "fit_k"
          
          df_k <- df_k %>% mutate(denom_inc = ifelse((is.na(a_k) | is.na(fit_k)), 0, 
                                                     (a_k - fit_k)*x_num*Si_censor))
          
          denom <- denom + sum(with(df_k, denom_inc))
        }
      }
      
      return(denom)
    })
  })
  
  return(jacobian)
}

newton_raphson_grad <- function(df, prm,
                                psi_start_vec=c(0,0), 
                                tol=0.001, max_iter=50, psi_max = 10,
                                print_results=F){
  
  psi_hat_vec <- psi_start_vec
  psi_old_vec <- psi_hat_vec
  steps <- 1L
  fail_flag <- F
  
  while(((sum(abs(psi_hat_vec - psi_old_vec)) > tol) | steps == 1) &
        (steps <= max_iter) & (max(abs(psi_hat_vec)) < psi_max)){
    
    
    print(steps)
    
    psi_old_vec <- psi_hat_vec
    
    (S_vec <- df %>% calculate_score(psi_hat_vec=psi_hat_vec, prm=prm))
    if(print_results){
      cat("\nS\n")
      print(S_vec)
    }
    (J <- df %>% calculate_jacobian(prm=prm, psi_hat_vec=psi_hat_vec))
    if(print_results){
      cat("\nJ\n")
      print(J)
    }
    (psi_hat_vec <- solve(J, S_vec) + psi_hat_vec)
    if(print_results){
      cat("\npsi_hat\n")
      print(psi_hat_vec)
    }
    if(print_results){cat("\n\n")}
    
    steps <- steps + 1
  }
  
  if ((max(abs(psi_hat_vec)) > psi_max) | (steps == max_iter)) {
    fail_flag <- T
    psi_hat_vec <- rep(NA, length(prm$psi_star_vec))
  }
  return(list(psi_hat_vec, steps, fail_flag))
}



calculate_score_piece <- function(df, psi_hat_vec, prm, i=1){
  
    S_curr <- 0
    
    for (a_curr in (0:(length(prm$t_a_vec)-1))){
      if (prm$beta_1_track[[a_curr+1]] == i) {
        
        df_temp <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=a_curr)
        
        names(df_temp)[names(df_temp)==paste0("a_", a_curr)] <- "a_k"
        names(df_temp)[names(df_temp)==paste0("fit_", a_curr)] <- "fit_k"
        
        df_temp <- df_temp %>% mutate(
          S_inc = ifelse((is.na(a_k) | is.na(fit_k) | is.na(tau_k)), 0,
                         (a_k - fit_k)*tau_k))
        
        #df_temp <- df %>% calculate_tau_rsd_m(., psi_hat_vec, prm)
        #names(df_temp)[names(df_temp)==paste0("a_", a_curr)] <- "a_k"
        #names(df_temp)[names(df_temp)==paste0("fit_", a_curr)] <- "fit_k"
        
        #df_temp <- df_temp %>% mutate(S_inc = ifelse(is.na(a_k), 0,
        #                                             (a_k - fit_k)*tau_rsd))
        S_curr <- S_curr + with(df_temp, sum(S_inc))
      }
    }

  return(S_curr)
  
}

calculate_jacobi_piece <- function(df, prm, psi_hat_vec, part_n){

  df_k_list <- lapply(0:(length(prm$t_a_vec)-1), 
                      FUN=function(k){
                        return(df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=k) %>%
                                 mutate(Si_temp = 0))
                      })
  
  row_counter <- part_n
  col_counter <- part_n
  
    
  if (row_counter <= max(prm$beta_1_track)) {
    d_num <- "theta_1"
    i <- row_counter
    beta_track_row <- prm$beta_1_track
  } else {
    d_num <- "theta_x"
    i <- row_counter - max(prm$beta_1_track)
    beta_track_row <- prm$beta_x_track
  }
  
  if (col_counter <= max(prm$beta_1_track)) {
    d_wrt <- "theta_1"
    j <- col_counter
    beta_track_col <- prm$beta_1_track
  } else {
    d_wrt <- "theta_x"
    j <- col_counter - max(prm$beta_1_track)
    beta_track_col <- prm$beta_x_track
  }
  
  denom <- 0
  
  for (a_curr in (0:(length(prm$t_a_vec)-1))){
    if (beta_track_row[[a_curr+1]] == i) {
      
      df_k <- df_k_list[[a_curr+1]]
      #df_k <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=a_curr) %>%
      #  mutate(Si_temp = 0)
      
      dC_dpsi <- 0
      for (m in (a_curr:(length(prm$t_a_vec)-1))) {
        if ((beta_track_row[[m+1]] == i) &&
            (beta_track_col[[m+1]] == j)) {
          
          df_k <- df_k %>% calculate_tau_rsd_m(prm=prm, psi_hat_vec=psi_hat_vec, m=m)
          
          names(df_k)[names(df_k)==paste0("a_", m)] <- "a_m"
          
          df_k <- df_k %>% mutate(
            x_wrt = (1 - as.integer(d_wrt == "theta_x")) + x * as.integer(d_wrt == "theta_x"),
            Si_inc = ifelse(is.na(a_m), 0, a_m * x_wrt * tau_rsd),
            Si_temp = Si_temp + Si_inc)
          
          names(df_k)[names(df_k)=="a_m"] <- paste0("a_", m)
          
          if(prm$censor){
            C_m <- calculate_C_m(psi_hat_vec=psi_hat_vec, prm=prm, m=m)
            dC_dpsi <- dC_dpsi + as.integer(psi_hat_vec[[j]] > 0)*C_m
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
          #C_k = calculate_C_k(psi_hat_vec=psi_hat_vec, prm=prm, a_k=a_curr)
        )
        #df_k <- df_k %>% mutate(Si_censor = delta_k*(psi_hat_vec[[j]] >= 0)*C_k + 
        #                          (1-delta_k)*Si_temp)
      } else {
        df_k <- df_k %>% mutate(
          Si_censor = Si_temp)
      }
      
      names(df_k)[names(df_k)==paste0("a_", a_curr)] <- "a_k"
      names(df_k)[names(df_k)==paste0("fit_", a_curr)] <- "fit_k"
      
      df_k <- df_k %>% mutate(denom_inc = ifelse((is.na(a_k) | is.na(fit_k)), 0, 
                                                 (a_k - fit_k)*x_num*Si_censor))
      
      denom <- denom + sum(with(df_k, denom_inc))
    }
  }
                                
  return(denom)

}

newton_raphson_stepwise <- function(df, prm,
                                psi_start_vec=c(0,0), 
                                tol=0.001, max_iter=50,
                                max_sub_iter=10,
                                psi_max = 10,
                                print_results=F){
  
  psi_hat_vec <- psi_start_vec
  psi_old_vec <- psi_hat_vec
  steps <- 1L
  fail_flag <- F
  
  sub_steps_vec <- numeric()
  
  #while(((sum(abs(psi_hat_vec - psi_old_vec)) > tol) | steps == 1) &
  #      (steps <= max_iter) & (max(abs(psi_hat_vec)) < psi_max)){
    
    for(n in max(prm$beta_1_track):1){
      cat("Solving for psi n = ", n, "\n\n")
      sub_steps <- 1L
      while (((abs(psi_hat_vec[n] - psi_old_vec[n]) > tol) | sub_steps == 1) &
             (sub_steps <= max_iter)) {
        
        cat("Iteration ", sub_steps, "\n\n")
        psi_old_vec <- psi_hat_vec
        
        (S_vec <- df %>% calculate_score(psi_hat_vec=psi_hat_vec, prm=prm))
        cat("\nS\n")
        print(S_vec)
        
        
        (J <- df %>% calculate_jacobian(prm=prm, psi_hat_vec=psi_hat_vec))
        cat("\nJ\n")
        print(J)
        
        (psi_hat_vec[[n]] <- solve(J[n,n], S_vec[n]) + psi_hat_vec[n])
        cat("\npsi_hat_old\n")
        print(psi_old_vec)
        cat("\npsi_hat_new\n")
        print(psi_hat_vec)
        cat("\n\n")

        
        sub_steps <- sub_steps + 1
      }
      
      sub_steps_vec <- c(sub_steps_vec, sub_steps)
    }
  
    
  return(list(psi_hat_vec, sub_steps_vec))
}

newton_raphson_piece <- function(df, prm,
                                    psi_start_vec=c(0,0), 
                                    tol=0.001, max_iter=50,
                                    max_sub_iter=10,
                                    psi_max = 10,
                                    print_results=F){
  
  psi_hat_vec <- psi_start_vec
  psi_old_vec <- psi_hat_vec
  steps <- 1L
  fail_flag <- F
  
  sub_steps_vec <- numeric()
  
  #while(((sum(abs(psi_hat_vec - psi_old_vec)) > tol) | steps == 1) &
  #      (steps <= max_iter) & (max(abs(psi_hat_vec)) < psi_max)){
  
  for(n in max(prm$beta_1_track):1){
    cat("Solving for psi n = ", n, "\n\n")
    sub_steps <- 1L
    while (((abs(psi_hat_vec[n] - psi_old_vec[n]) > tol) | sub_steps == 1) &
           (sub_steps <= max_iter)) {

      psi_old_vec <- psi_hat_vec
      
      (S_piece <- df %>% calculate_score_piece(psi_hat_vec=psi_hat_vec, 
                                               prm=prm, i=n))
      #cat("\nS")
      #print(S_piece)
      
      
      (J_piece <- df %>% calculate_jacobi_piece(prm=prm, psi_hat_vec=psi_hat_vec,
                                          part_n=n))
      #cat("\nJ")
      #print(J_piece)
      
      psi_hat_vec[[n]] <- solve(J_piece, S_piece) + psi_hat_vec[n]
      
      
      cat(paste0("psi_", n), " - Iteration", sub_steps, "\n")
      cat("\npsi_hat_old\n")
      print(psi_old_vec)
      cat("\npsi_hat_new\n")
      print(psi_hat_vec)
      cat("\n\n")
      
      
      sub_steps <- sub_steps + 1
    }
    
    sub_steps_vec <- c(sub_steps_vec, sub_steps)
  }
  
  S_start <-  df %>% calculate_score(psi_hat_vec=psi_start_vec, prm=prm)
  S_final <- df %>% calculate_score(psi_hat_vec=psi_hat_vec, prm=prm)
  
  cat("S - start vs final")
  print(S_start)
  print(S_final)
  
  
  return(list(psi_hat_vec, sub_steps_vec))
}

calculate_variance_SLOW <- function(df, prm, psi_hat_vec){
  
  n_vec <- c()
  for (t_a in prm$t_a_vec){
    n_vec <- c(n_vec, nrow(df %>% filter(.$ti > t_a)))
  }
  
  progressbar <- txtProgressBar(min=0, max=length(prm$t_a_vec)-1,
                                style=3, width=50, char="=")
  D_vec <- c()
  col_count <- 0
  for (k in 0:(length(prm$t_a_vec)-1)) {
    setTxtProgressBar(progressbar, k)
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
  close(progressbar)
  
  
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


calculate_variance <- function(df, prm, psi_hat_vec, trt_models,
                               suppress_progress=F){
  
  #n_vec <- c()
  n_vec <- sapply(prm$t_a_vec, 
                  function(t_a) nrow(df %>% filter(.$ti > t_a)))
  
  design_matrices <- lapply(trt_models,
                    function(model_temp){
                      mdl_mtrx <- model.matrix(model_temp)
                      return(mdl_mtrx %>% StripAttr() %>% matrix(., nrow=dim(mdl_mtrx)[[1]])) 
                    })
  
  D_cols <- sapply(design_matrices, 
                   function(D_matrix){
                     return(D_matrix %>% ncol())
                   })
  D_rows <- sapply(design_matrices, 
                   function(D_matrix){
                     return(D_matrix %>% nrow())
                   })
  
  if(!suppress_progress){cat("Calculating D\n")}
  D <- sapply(0:(length(prm$t_a_vec)-1),
       function(j) {
         sapply(0:(length(prm$t_a_vec)-1),
          function(i) {
            if (i == j) {
              return(design_matrices[[j+1]])
            } else {
              return(matrix(0, nrow=D_rows[[i+1]], ncol=D_cols[[j+1]]))
            }
          }) %>% Reduce(rbind, .)
       }) %>% Reduce(cbind, .)
      
  if(!suppress_progress){cat("Calculating D_theta\n")}
  #for (beta_1_val in 1:max(prm$beta_1_track))
  D_theta <- sapply(1:max(prm$beta_1_track),
                    function(beta_1_val){
                      #for (k in 0:(length(prm$t_a_vec)-1))
                      sapply(0:(length(prm$t_a_vec)-1), beta_1_val = beta_1_val,
                             FUN = function(k, beta_1_val){
                               if (prm$beta_1_track[[k+1]] == beta_1_val) {
                                 return(rep(1, n_vec[[k+1]]))
                               } else {
                                 return(rep(0, n_vec[[k+1]]))
                               }
                             }) %>% Reduce(c, .) %>% matrix(ncol=1)
                    })
  if (max(prm$beta_x_track) > 0) {
    D_theta <- D_theta %>% cbind(
      #for (beta_x_val in 1:max(prm$beta_x_track))
      sapply(1:max(prm$beta_x_track), 
             FUN = function(beta_x_val){
               #for (k in 0:(length(prm$t_a_vec)-1))
               sapply(0:(length(prm$t_a_vec)-1), beta_x_val = beta_x_val,
                      FUN = function(k, beta_x_val) {
                        if (prm$beta_x_track[[k+1]] == beta_x_val) {
                          df_temp <- df %>% filter(.$ti > prm$t_a_vec[[k+1]])
                          return(df_temp$x)
                        } else {
                          return(rep(0, n_vec[[k+1]]))
                        }          
                      }) %>% Reduce(c, .) %>% matrix(ncol=1)
               
             })
    )
  }
  
  
  if(!suppress_progress){cat("Calculating fit and tau\n")}
  fit <- c()
  tau <- c()
  for (k in 0:(length(prm$t_a_vec)-1)){
    df_temp <- df %>% filter( .$ti > prm$t_a_vec[[k+1]]) %>% 
      calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=k)
    fit <- c(fit, df_temp[[paste0('fit_', k)]])
    tau <- c(tau, df_temp$tau_k)
  }
  
  
  
  #Jbb <- apply(D, 2,
  #        function(col_i){
  #          apply(D, 2,
  #          function(col_j){
  #            sum(col_i*col_j*fit*(1-fit))    
  #          }) %>% matrix(nrow=1)
  #               
  #        })
  if(!suppress_progress){cat("Calculating Jbb\n")}
  start.time <- Sys.time()
  Jbb <- t(sapply(1:dim(D)[[2]], function(i){
    sapply(1:dim(D)[[2]], function(j){
      sum(D[,i]*D[,j]*fit*(1-fit), na.rm=T)
    })
  }))
  end.time <- Sys.time()
  if(!suppress_progress){cat(end.time - start.time, "\n")}
  
  if(!suppress_progress){cat("Calculating Jtt\n")}
  Jtt <- t(sapply(1:dim(D_theta)[[2]], function(i){
    sapply(1:dim(D_theta)[[2]], function(j){
      sum(D_theta[,i]*D_theta[,j]*fit*(1-fit)*tau*tau, na.rm=T)
    })
  }))
  
  if(!suppress_progress){cat("Calculating Jtb\n")}
  Jtb <- t(sapply(1:dim(D_theta)[[2]], function(i){
    sapply(1:dim(D)[[2]], function(j){
      sum(D_theta[,i]*D[,j]*fit*(1-fit)*tau, na.rm=T)
    })
  }))
  
  if(!suppress_progress){cat("Calculating JTT\n")}
  JTT <- Jtt - Jtb%*%solve(Jbb)%*%t(Jtb)
  
  
  if(!suppress_progress){cat("Calculating VCV\n")}
  Jacobian <- df %>% calculate_jacobian(prm=prm, psi_hat_vec=psi_hat_vec)
  J_inv <- solve(Jacobian)
  
  (VCV <- J_inv %*% JTT %*% t(J_inv))
  
  return(VCV)
}
