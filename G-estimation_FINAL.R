library(MASS)
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
  
  trt_models <- list()
  trt_models <- lapply(0:(length(prm$trt_mod_list)-1),
      function(k){

   # for(k in 0:(length(prm$trt_mod_list)-1)){
        mdl_formula <- paste0("a_", k, " ~ ", 
                              paste0(prm$trt_mod_list[[k+1]], collapse = " + "))
        
        #tryCatch({
        df_k <- df %>% filter(.$ti > prm$t_a_vec[[k+1]])
        model_temp <- with(df_k, 
                           glm(  as.formula(mdl_formula),
                                 family=binomial))
        
      })
  
  
  
  for(k in 0:(length(trt_models)-1)){
    df_k <- df %>% filter(.$ti > prm$t_a_vec[[k+1]])
    model_temp <- trt_models[[k+1]]
    df_k$fit_temp <- predict(model_temp, newdata = df_k, type="response")
    df_k <- df_k %>% dplyr::select(c(id, fit_temp))
    df <- df %>% left_join(df_k, by="id")
    names(df)[names(df) == "fit_temp"] <- paste0("fit_", k)
  }
  
  return(list(df=df, trt_models=trt_models))
}

fit_treatment_diagnostic <- function(df, prm){
  
  trt_models <- list()
  # trt_models <- lapply(0:(length(prm$trt_mod_list)-1),
  #                      function(k){
                         
  for(k in 0:(length(prm$trt_mod_list)-1)){
    mdl_formula <- paste0("a_", k, " ~ ", 
                          paste0(prm$trt_mod_list[[k+1]], collapse = " + "))
    
    #tryCatch({
    df_k <- df %>% filter(.$ti > prm$t_a_vec[[k+1]])
    model_temp <- with(df_k, 
                       glm(  as.formula(mdl_formula),
                             family=binomial))
  }
  # })
  
  
  
  for(k in 0:(length(trt_models)-1)){
    df_k <- df %>% filter(.$ti > prm$t_a_vec[[k+1]])
    model_temp <- trt_models[[k+1]]
    df_k$fit_temp <- predict(model_temp, newdata = df_k, type="response")
    df_k <- df_k %>% dplyr::select(c(id, fit_temp))
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
      C_rsd = pmax(C_i - t_a, 0),
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
    
    a_curr <- as.name(paste0("a_", i-1))
    
    df <- df %>% mutate(
      g_psi =  beta_1_now + x*beta_x_now,
      ti_temp = pmin(t_next - t_curr, ti_rsd),
      tau_k = ifelse(is.na({{a_curr}}), tau_k, tau_k + ti_temp*exp(-g_psi*{{a_curr}})),
      ti_rsd = ti_rsd - ti_temp,
      
    )
    
    if (prm$censor) {
      t_next <- ifelse(i < length(prm$t_a_vec), 
                       prm$t_a_vec[[i+1]], 
                       prm$censor_max)
      df <- df %>% mutate(
        #C_psi = exp(-abs(g_psi)),
        C_psi = ifelse(g_psi < 0, 1, exp(-g_psi)),
        C_k = C_k + pmin(t_next - t_curr, C_rsd)*C_psi,
        C_rsd = C_rsd - pmin(t_next - t_curr, C_rsd)
      )
    }

    
  }
  if (prm$censor) {
    df <- df %>% mutate(#tau_k_un = tau_k,
                        delta_k = ifelse(tau_k < C_k, 0, 1),
                        tau_k = pmin(tau_k, C_k))
    return(dplyr::select(df, -c(ti_temp, ti_rsd, g_psi, C_rsd, C_k, C_psi)))
  } else {
    return(dplyr::select(df, -c(ti_temp, ti_rsd, g_psi)))
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
    
    return(df %>% dplyr::select(-c(g_psi, ti_temp)))
}

calculate_C_m <- function(df, psi_hat_vec, prm, m=0){
  #,a_k=0, j=1){
  
  
  
  
  t_now <- prm$t_a_vec[[m+1]]
  t_next <- ifelse(m+1 < length(prm$t_a_vec), 
                   prm$t_a_vec[[m+2]], 
                   prm$censor_date)
  
  beta_1_now <- ifelse(prm$beta_1_track[[m+1]]==0, 0,
                       psi_hat_vec[[prm$beta_1_track[[m+1]]]])
  beta_x_now <- ifelse(prm$beta_x_track[[m+1]]==0, 0,
                       psi_hat_vec[[prm$beta_x_track[[m+1]]+max(prm$beta_1_track)]])
  
  df <- df %>% mutate(C_check = ifelse(C_i < t_next, pmax(C_i - t_now, 0), t_next - t_now),
                      C_m = exp(-(beta_1_now*(beta_1_now>=0) + beta_x_now*(beta_x_now>=0)))*C_check)
  
  #C_m <- exp(-(beta_1_now*(beta_1_now>=0) + beta_x_now*(beta_x_now>=0)))*(t_next - t_now)
  
  #beta_x_now <- ifelse(beta_x_track[[m+1]]==0, 0,
  #                     psi_hat_vec[[beta_x_track[[m+1]]+length(psi_1_star)]])
  #if (beta_1_now >= 0){
  
  #}
  
  return(dplyr::select(df, -c(C_check)))
}

calculate_C_k <- function(df, psi_hat_vec, prm, a_k=0){
  
  df <- df %>% mutate(C_k = 0)
  for (m in (a_k:(length(prm$t_a_vec)-1))) {
    df <- df %>% calculate_C_m(psi_hat_vec=psi_hat_vec, prm=prm, m=m) %>%
                mutate(C_k = C_k + C_m)
  }
  return(df)
  
}



calculate_score <- function(df, psi_hat_vec, prm){
  
  S_vec <- c()
  for (i in 1:max(prm$beta_1_track)){
    
    S_curr <- 0
    
    for (a_curr in (0:(length(prm$t_a_vec)-1))){
      if (prm$beta_1_track[[a_curr+1]] == i) {
        
        df_temp <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=a_curr)
        
        a_k <- as.name(paste0("a_", a_curr))
        fit_k <- as.name(paste0("fit_", a_curr))
        #names(df_temp)[names(df_temp)==paste0("a_", a_curr)] <- "a_k"
        #names(df_temp)[names(df_temp)==paste0("fit_", a_curr)] <- "fit_k"
        
        df_temp <- df_temp %>% mutate(
                    S_inc = ifelse((is.na({{a_k}}) | is.na({{fit_k}}) | is.na(tau_k)), 0,
                                                     ({{a_k}} - {{fit_k}})*tau_k))
        
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
          
          
          a_k <- as.name(paste0("a_", a_curr))
          fit_k <- as.name(paste0("fit_", a_curr))
          #names(df_temp)[names(df_temp)==paste0("a_", a_curr)] <- "a_k"
          #names(df_temp)[names(df_temp)==paste0("fit_", a_curr)] <- "fit_k"
          
          df_temp <- df_temp %>% mutate(
                    S_inc = ifelse((is.na({{a_k}}) | is.na({{fit_k}}) | is.na(tau_k)), 0,
                                                     ({{a_k}} - {{fit_k}})*x*tau_k))
        
          S_curr <- S_curr + with(df_temp, sum(S_inc))          
          
          #S_curr <- S_curr + with(df_temp, sum((a_k - fit_k)*x*tau_k))
        }
      }
      
      S_vec <- c(S_vec, S_curr)
      
    }
  }
  return(S_vec)
  
}

calculate_score_piece <- function(df, psi_hat_vec, prm, d=1){
  
  S_curr <- 0
  
  for (a_curr in (0:(length(prm$t_a_vec)-1))){
    if (prm$beta_1_track[[a_curr+1]] == d) {
      
      df_temp <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=a_curr)
      
      a_k <- as.name(paste0("a_", a_curr))
      fit_k <- as.name(paste0("fit_", a_curr))
      
      df_temp <- df_temp %>% mutate(
        #S_i := ({{a_k}} - {{fit_k}})*tau_k
        S_i := ifelse((is.na({{a_k}}) | is.na({{fit_k}}) | is.na(tau_k)), 0,
                        ({{a_k}} - {{fit_k}})*tau_k)
        )
      
      S_inc <- df_temp %>% dplyr::select(S_i) %>% sum()
      S_curr <- S_curr + S_inc
    }
  }
  
  return(S_curr)
  
}

calculate_jacobian <- function(df, prm, psi_hat_vec){
  jacobi_vec <- c()
  jacobi_dim <- length(psi_hat_vec)

  # df_k_list <- lapply(0:(length(prm$t_a_vec)-1), 
  #                     FUN=function(k){
  #                       return(df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=k))
  #                     })
  
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
          
          df_k <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=a_curr) %>%
                        mutate(Si_temp = 0, dC_dpsi = 0)

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
                df_k <- df_k %>% calculate_C_m(psi_hat_vec=psi_hat_vec, prm=prm, m=m)
                df_k <- df_k %>% mutate(
                          dC_dpsi = dC_dpsi + as.integer(psi_hat_vec[[j]] > 0)*C_m
                        )
              }
            }
          }
          
          df_k <- df_k %>% mutate(
            x_num = (1 - as.integer(d_num == "theta_x")) + x * as.integer(d_num == "theta_x"),
          )
          if(prm$censor) { 
            df_k <- df_k %>% mutate(
                      Si_censor = delta_k*dC_dpsi + (1-delta_k)*Si_temp
                      )
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




calculate_jacobi_piece <- function(df, prm, psi_hat_vec, d){

  row_counter <- d
  col_counter <- d
  
    
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
      
      df_k <- df %>% 
                  calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=a_curr) %>%
                  mutate(Si_temp = 0, 
                         dC_dpsi = 0)
      fit_k <- as.name(paste0("fit_", a_curr))
      for (m in (a_curr:(length(prm$t_a_vec)-1))) {
        if ((beta_track_row[[m+1]] == i) &&
            (beta_track_col[[m+1]] == j)) {
          
          df_k <- df_k %>% 
                    calculate_tau_rsd_m(prm=prm, psi_hat_vec=psi_hat_vec, m=m)
          
          a_m <- as.name(paste0("a_", m))
          a_mprev <- as.name(paste0("a_", m-1))
          #names(df_k)[names(df_m)==] <- "a_m"
          
          #df_k <- df_k %>%
          #  dplyr::select(id,{{a_mprev}}, {{a_m}}, x, tau_rsd, Si_temp, dC_dpsi, {{fit_k}})
          
          df_k <- df_k %>% mutate(
            x_wrt = (1 - as.integer(d_wrt == "theta_x")) + x*as.integer(d_wrt == "theta_x"),
            Si_inc = ifelse(is.na({{a_m}}), 0, {{a_m}}*x_wrt*tau_rsd),
            Si_temp = Si_temp + Si_inc)
          
          if(prm$censor){
            df_k <- df_k %>% calculate_C_m(psi_hat_vec=psi_hat_vec, prm=prm, m=m)
            df_k <- df_k %>% mutate(
                      dC_dpsi = dC_dpsi + as.integer(psi_hat_vec[[j]] > 0)*C_m
                    )
          }
        }
      }
      
      df_k <- df_k %>% mutate(
        x_num = (1 - as.integer(d_num == "theta_x")) + x * as.integer(d_num == "theta_x"))
      
      if(prm$censor) { 
        df_k <- df_k %>% mutate(
                    Si_censor = delta_k*dC_dpsi + (1-delta_k)*Si_temp)

      } else {
        df_k <- df_k %>% mutate(
          Si_censor = Si_temp)
      }
      a_k <- as.name(paste0("a_", a_curr))
      fit_k <- as.name(paste0("fit_", a_curr))
      # names(df_k)[names(df_k)==paste0("a_", a_curr)] <- "a_k"
      # names(df_k)[names(df_k)==paste0("fit_", a_curr)] <- "fit_k"
      
      df_k <- df_k %>% mutate(denom_inc = ifelse((is.na({{a_k}}) | is.na({{fit_k}})), 0, 
                                                 ({{a_k}} - {{fit_k}})*x_num*Si_censor))
      
      denom <- denom + sum(with(df_k, denom_inc))
    }
  }
                                
  return(denom)

}



newton_raphson_grad <- function(df, prm,
                                psi_start_vec=NA, 
                                tol=0.001, min_damp=0.1,
                                max_iter=50, 
                                psi_max = 10,
                                psi_max_vec = NA,
                                print_results=F){
  
  if(all(is.na(psi_start_vec))){
    psi_start_vec <- rep(0, max(prm$beta_1_track)+max(prm$beta_x_track))
  }
  
  if(all(is.na(psi_max_vec))){
    psi_max_vec = rep(psi_max, length(psi_start_vec))
  }
  
  psi_hat_vec <- psi_start_vec
  psi_old_vec <- psi_hat_vec
  steps <- 1L
  fail_flag <- F
  
  while(((sum(abs(psi_hat_vec - psi_old_vec)) > tol) |
                steps == 1 |
                (sum(abs(psi_hat_vec)) < tol)) &
        (steps <= max_iter) &
        (max(abs(psi_hat_vec)) < psi_max)){
    
    if(print_results==T){
      cat("(", steps, ")", " --- ", psi_old_vec, "\n")
    } else {
      print(steps)
    }
    
    psi_old_vec <- psi_hat_vec
    
    (S_vec <- df %>% calculate_score(psi_hat_vec=psi_hat_vec, prm=prm))
    # if(print_results){
    #   cat("\nS\n")
    #   print(S_vec)
    # }
    (J <- df %>% calculate_jacobian_alt(prm=prm, psi_hat_vec=psi_hat_vec))
    # if(print_results){
    #   cat("\nJ\n")
    #   print(J)
    # }
    damping_factor <- max(2^(-(steps - 1)/20), min_damp)
    (psi_hat_vec <- solve(J, S_vec)*damping_factor + psi_hat_vec)
    # if(print_results){
    #   cat("\npsi_hat\n")
    #   print(psi_hat_vec)
    # }
    # if(print_results){cat("\n\n")}
    
    if(!all(abs(psi_hat_vec) < psi_max_vec)){
      psi_hat_vec <- psi_start_vec
    }
    
    
    # psi_hat_vec <- sapply(1:length(psi_hat_vec), function(d){
    #                   psi_bit <- psi_hat_vec[d]
    #                   if(psi_bit > -psi_max_vec[d] & psi_bit < psi_max_vec[d]){
    #                     return(psi_bit)
    #                   } else {
    #                     return(0)
    #                   }
    #                 })
    
    steps <- steps + 1
  }
  
  if ((max(abs(psi_hat_vec)) > psi_max) | (steps >= max_iter)) {
    fail_flag <- T
    psi_hat_vec <- rep(NA, length(prm$psi_star_vec))
  }
  
  return(list(psi_hat_vec, steps, fail_flag))
}

newton_raphson_stepwise <- function(df, prm,
                                psi_start_vec=c(0,0), 
                                tol=0.001, max_iter=50,
                                max_sub_iter=10,
                                psi_max = 10,
                                print_results=T){
  
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
                                    psi_start_vec=NA, 
                                    tol=0.001, 
                                    max_iter=100,
                                    min_damp=0.1,
                                    psi_max = 10,
                                    psi_max_vec = NA,
                                    print_results=T){
  
  if(all(is.na(psi_start_vec))){
    psi_start <- rep(0, max(prm$beta_1_track)+max(prm$beta_x_track))
  } else {
    psi_start <- psi_start_vec
  }
  
  psi_hat_vec <- psi_start
  psi_old_vec <- psi_hat_vec
  steps <- 1L
  fail_flag <- F
  damping_factor <- 1
  
  sub_steps_vec <- numeric()
  
  #while(((sum(abs(psi_hat_vec - psi_old_vec)) > tol) | steps == 1) &
  #      (steps <= max_iter) & (max(abs(psi_hat_vec)) < psi_max)){
  
  for(d in max(prm$beta_1_track):1){
    #cat("Solving for psi n = ", n, "\n\n")
    sub_steps <- 1L
    if(all(is.na(psi_max_vec))){
      psi_max_temp <- psi_max
    } else {
      psi_max_temp <- psi_max_vec[d]
    }
    while (((abs(psi_hat_vec[d] - psi_old_vec[d]) > tol*damping_factor) | (sub_steps == 1) |
            (psi_hat_vec[d] < -psi_max_temp+2*tol) | (psi_hat_vec[d] > psi_max_temp - 2*tol)) &
           (sub_steps <= max_iter)) {

      if(print_results){cat(paste0("psi_", d), " - ", sub_steps, "\n")}
      psi_old_vec <- psi_hat_vec
      
      (S_piece <- df %>% calculate_score_piece(psi_hat_vec=psi_hat_vec, 
                                               prm=prm, d=d))
      #cat("\nS")
      #print(S_piece)
      
      
      (J_piece <- df %>% calculate_jacobi_piece(prm=prm, psi_hat_vec=psi_hat_vec,
                                          d=d))
      #cat("\nJ")
      #print(J_piece)
      
      damping_factor <- max(2^(-(sub_steps - 1)/20), min_damp)
      
      (psi_hat_vec[[d]] <- solve(J_piece, S_piece)*damping_factor + psi_hat_vec[d])
      
      psi_hat_vec[d] <- max(psi_hat_vec[d], -psi_max_temp)
      psi_hat_vec[d] <- min(psi_hat_vec[d], psi_max_temp)
      
      if((psi_hat_vec[d] == psi_old_vec[d]) & (abs(psi_hat_vec[d]) == psi_max_temp)){
        psi_hat_vec[d] <- 0
      }
      
      
      psi_piece <- psi_hat_vec[[d]]      
      if(print_results){cat(psi_piece, "\n")}
      #cat("\npsi_hat_old\n")
      #print(psi_old_vec)
      #cat("\npsi_hat_new\n")
      #print(psi_hat_vec)
      #cat("\n\n")
      
      
      sub_steps <- sub_steps + 1
    }
    
    if(print_results){cat("\n", psi_hat_vec, "\n\n")}
    
    sub_steps_vec <- c(sub_steps_vec, sub_steps)
  }
  
  #S_start <-  df %>% calculate_score(psi_hat_vec=psi_start_vec, prm=prm)
  #S_final <- df %>% calculate_score(psi_hat_vec=psi_hat_vec, prm=prm)
  
  #cat("S - start vs final\n")
  #print(S_start)
  #print(S_final)
  
  
  return(list(psi_hat_vec, sub_steps_vec))
}


calculate_variance <- function(df, prm, psi_hat_vec,# trt_models,
                               print_results=T){
  
  #n_vec <- c()
  if(print_results){cat("\nCalculating variance\n")}
  n_vec <- sapply(prm$t_a_vec, 
                  function(t_a) nrow(df %>% filter(.$ti > t_a)))
  
  design_matrices <- lapply(prm$trt_models,
                            function(model_temp){
                              mdl_mtrx <- model.matrix(model_temp)
                              return(mdl_mtrx %>% StripAttr() %>% matrix(., nrow=dim(mdl_mtrx)[[1]])) 
                            })
  # design_dims <- lapply(1:length(design_matrices),
  #                function(k){
  #                  return(tibble(nrows=dim(design_matrices[[k]])[1],
  #                                ncols=dim(design_matrices[[k]])[2]))
  #                }) %>% Reduce(rbind, .)
  
  
  D_cols <- sapply(design_matrices, 
                   function(D_matrix){
                     return(D_matrix %>% ncol())
                   })
  
  D_rows <- sapply(design_matrices, 
                   function(D_matrix){
                     return(D_matrix %>% nrow())
                   })
  
  if(print_results){cat("Calculating D\n")}
  start.time <- Sys.time()
  D <- lapply(0:(length(prm$t_a_vec)-1),
              function(j) {
                lapply(0:(length(prm$t_a_vec)-1),
                       function(i) {
                         if (i == j) {
                           return(design_matrices[[j+1]])
                         } else {
                           return(matrix(0, nrow=D_rows[[i+1]], ncol=D_cols[[j+1]]))
                         }
                       }) %>% Reduce(rbind, .)
              }) %>% Reduce(cbind, .)
  if(print_results){cat(Sys.time() - start.time,"\n")}
  
  if(print_results){cat("Calculating D_theta\n")}
  start.time <- Sys.time()
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
  if(print_results){cat(Sys.time() - start.time,"\n")}
  
  
  if(print_results){cat("Calculating fit and tau\n")}
 start.time <- Sys.time()
  fit <- c()
  tau <- c()
  for (k in 0:(length(prm$t_a_vec)-1)){
    df_temp <- df %>% filter( .$ti > prm$t_a_vec[[k+1]]) %>% 
      calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=k)
    fit <- c(fit, df_temp[[paste0('fit_', k)]])
    tau <- c(tau, df_temp$tau_k)
  }
  fit_list <- lapply(0:(length(prm$t_a_vec)-1),
                function(k){
                  df_temp <- df %>% filter( .$ti > prm$t_a_vec[[k+1]]) %>%
                    calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=k)
                  return(df_temp[[paste0('fit_', k)]])
                })
  if(print_results){cat(Sys.time() - start.time,"\n")}
  
  
  
  #Jbb <- apply(D, 2,
  #        function(col_i){
  #          apply(D, 2,
  #          function(col_j){
  #            sum(col_i*col_j*fit*(1-fit))    
  #          }) %>% matrix(nrow=1)
  #               
  #        })
  #if(!suppress_progress){cat("Calculating Jbb\n")}
  
  # start.time <- Sys.time()
  # Jbb <- t(sapply(1:dim(D)[[2]], function(i){
  #   sapply(1:dim(D)[[2]], function(j){
  #     sum(D[,i]*D[,j]*fit*(1-fit), na.rm=T)
  #   })
  # }))
  # end.time <- Sys.time()
  
  
  #start.time <- Sys.time()
  if(print_results){cat("Calculating Jbb\n")}
  start.time <- Sys.time()
  Jbb_list <- lapply(1:length(design_matrices), function(k){
    sapply(1:dim(design_matrices[[k]])[[2]], function(i){
      sapply(1:dim(design_matrices[[k]])[[2]], function(j){
        fit_k <- fit_list[[k]]
        design_k <- design_matrices[[k]]
        return(sum(design_k[,i]*design_k[,j]*fit_k*(1-fit_k)))
      })
    })
  })
  
  Jbb_dims <- sapply(1:length(Jbb_list),
                function(k){
                  Jbb_part <- Jbb_list[[k]] %>% as.matrix()
                  return(dim(Jbb_part)[1])
                })
  
  Jbb <- lapply(1:length(Jbb_list),
          function(i){
            lapply(1:length(Jbb_list),
            function(j){
              if(j==i){
                return(Jbb_list[[i]])
              } else{
                return(matrix(rep(0, Jbb_dims[i]*Jbb_dims[j]),
                              nrow=Jbb_dims[i]))
              }
            }) %>% Reduce(cbind, .)
          }) %>% Reduce(rbind, .)
  if(print_results){cat(Sys.time() - start.time,"\n")}
  #end.time <- Sys.time()
  #cat(end.time - start.time, "\n")
  #if(!suppress_progress){cat(end.time - start.time, "\n")}
  
  if(print_results){cat("Calculating Jtt\n")}
  start.time <- Sys.time()
  Jtt <- t(sapply(1:dim(D_theta)[[2]], function(i){
    sapply(1:dim(D_theta)[[2]], function(j){
      sum(D_theta[,i]*D_theta[,j]*fit*(1-fit)*tau*tau, na.rm=T)
    })
  }))
  if(print_results){cat(Sys.time() - start.time,"\n")}
  
  if(print_results){cat("Calculating Jtb\n")}
  start.time <- Sys.time()
  Jtb <- lapply(1:dim(D_theta)[[2]], function(i){
          lapply(1:dim(D)[[2]], function(j){
      sum(D_theta[,i]*D[,j]*fit*(1-fit)*tau, na.rm=T)
    }) %>% Reduce(cbind, .)
  }) %>% Reduce(rbind, .)
  if(print_results){cat(Sys.time() - start.time,"\n")}
  
  if(print_results){cat("Calculating JTT\n")}
  start.time <- Sys.time()
  Jbb_inv <- solve(Jbb)
  if(dim(Jbb_inv)[1] == 1 && dim(Jbb_inv)[2] == 1){
    JTT <- Jtt - Jtb %*% t(Jtb) * as.numeric(Jbb_inv)
  } else {
    JTT <- Jtt - Jtb%*%solve(Jbb)%*%t(Jtb)
  }
  if(print_results){cat(Sys.time() - start.time,"\n")}
  
  if(print_results){cat("Calculating VCV\n")}
  start.time <- Sys.time()
  Jacobian <- df %>% calculate_jacobian(prm=prm, psi_hat_vec=psi_hat_vec)
  J_inv <- solve(Jacobian)
  (VCV <- J_inv %*% JTT %*% t(J_inv))
  if(print_results){cat(Sys.time() - start.time,"\n")}
  
  return(VCV)
}






#### LARGE VARIANT


calculate_score_large <- function(df, psi_hat_vec, prm){
  
  S_vec <-  lapply((0:(length(prm$t_a_vec)-1)),
                   function(a_curr){
                     df_temp <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=a_curr)
                     
                     a_k <- as.name(paste0("a_", a_curr))
                     fit_k <- as.name(paste0("fit_", a_curr))
                     
                     df_temp <- df_temp %>% mutate(
                       S_inc = ifelse((is.na({{a_k}}) | is.na({{fit_k}}) | is.na(tau_k)), 0,
                                      ({{a_k}} - {{fit_k}})*tau_k))
                     
                     return(with(df_temp, sum(S_inc)))
                     
                   }) %>% Reduce(rbind, .)

  return(S_vec)
  
}

calculate_jacobian_large <- function(df, prm, psi_hat_vec){
  jacobi_vec <- c()
  jacobi_dim <- length(psi_hat_vec)

  jacobian <- lapply(0:(length(prm$t_a_vec) -1),
              function(a_curr){
                lapply(1:max(prm$beta_1_track),
                function(psi_val){
                  df_k <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=a_curr) %>%
                        mutate(Si_temp = 0, dC_dpsi = 0)
                  
                  for(m in (a_curr:(length(prm$t_a_vec)-1))){
                    if (prm$beta_1_track[[m+1]] == psi_val){
                      df_k <- df_k %>% calculate_tau_rsd_m(prm=prm, psi_hat_vec=psi_hat_vec, m=m)

                      a_m <- as.name(paste0("a_", m))
                      df_k <- df_k %>% mutate(
                        x_wrt = 1,
                        Si_inc = ifelse(is.na({{a_m}}), 0, {{a_m}} * x_wrt * tau_rsd),
                        Si_temp = Si_temp + Si_inc)
                      
                      
                      if(prm$censor){
                        df_k <- df_k %>% calculate_C_m(psi_hat_vec=psi_hat_vec, prm=prm, m=m)
                        df_k <- df_k %>% mutate(
                          dC_dpsi = dC_dpsi + as.integer(psi_hat_vec[[psi_val]] > 0)*C_m
                        )
                      }
                  
                      
                    }
                  }
                  
                  if(prm$censor) { 
                    df_k <- df_k %>% mutate(
                      Si_censor = delta_k*dC_dpsi + (1-delta_k)*Si_temp
                    )
                  } else {
                    df_k <- df_k %>% mutate(
                      Si_censor = Si_temp)
                  }
                  
                  
                  a_k <- as.name(paste0("a_", a_curr))
                  fit_k <- as.name(paste0("fit_", a_curr))
                  df_k <- df_k %>% mutate(denom_inc = ifelse((is.na({{a_k}}) | is.na({{fit_k}})), 0, 
                                                             ({{a_k}} - {{fit_k}})*Si_censor))
                  
                  return(sum(with(df_k, denom_inc)))
                  }) %>% Reduce(cbind, .)
                }) %>% Reduce(rbind, .)
  
  return(jacobian)
}

newton_raphson_grad_large <- function(df, prm,
                                psi_start_vec=NA, 
                                tol=0.001, min_damp=0.1,
                                max_iter=50, 
                                psi_max = 10,
                                psi_max_vec = NA,
                                print_results=F){
  
  if(all(is.na(psi_start_vec))){
    psi_start_vec <- rep(0, max(prm$beta_1_track)+max(prm$beta_x_track))
  }
  
  if(all(is.na(psi_max_vec))){
    psi_max_vec = rep(psi_max, length(psi_start_vec))
  }
  
  psi_hat_vec <- psi_start_vec
  #S_vec <- df %>% calculate_score_alt(psi_hat_vec=psi_hat_vec, prm=prm)
  psi_old_vec <- psi_hat_vec
  steps <- 1L
  fail_flag <- F
  
  while(((sum(abs(psi_hat_vec - psi_old_vec))/(length(psi_hat_vec)) > tol) |
                steps == 1 |
                (sum(abs(psi_hat_vec)) < tol)) &
        (steps <= max_iter) &
        (max(abs(psi_hat_vec)) < psi_max)){
    
    if(print_results==T){
      cat("(", steps, ")", " --- ", psi_old_vec, "\n")
    } else {
      print(steps)
    }
    
    psi_old_vec <- psi_hat_vec
    
    (S_vec <- df %>% calculate_score_large(psi_hat_vec=psi_hat_vec, prm=prm))

    (J <- df %>% calculate_jacobian_large(prm=prm, psi_hat_vec=psi_hat_vec))

    damping_factor <- max(2^(-(steps - 1)/20), min_damp)
    (psi_hat_vec <- solve(t(J) %*% J) %*% t(J)  %*% S_vec *damping_factor + psi_hat_vec)
    # (psi_hat_vec <- solve(t(J) %*% J) %*% t(J)  %*% S_vec *0.2 + psi_hat_vec)

    
    if(!all(abs(psi_hat_vec) < psi_max_vec)){
      psi_hat_vec <- psi_start_vec
    }
    
    steps <- steps + 1
  }
  
  if ((max(abs(psi_hat_vec)) > psi_max) | (steps >= max_iter)) {
    fail_flag <- T
    psi_hat_vec <- rep(NA, length(prm$psi_star_vec))
  }
  
  return(list(psi_hat_vec, steps, fail_flag))
}

calculate_variance_large <- function(df, prm, psi_hat_vec,# trt_models,
                               print_results=T){
  

  n_vec <- sapply(prm$t_a_vec, 
                  function(t_a) nrow(df %>% filter(.$ti > t_a)))
  
  design_matrices <- lapply(prm$trt_models,
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
  
  D <- lapply(0:(length(prm$t_a_vec)-1),
              function(j) {
                lapply(0:(length(prm$t_a_vec)-1),
                       function(i) {
                         if (i == j) {
                           return(design_matrices[[j+1]])
                         } else {
                           return(matrix(0, nrow=D_rows[[i+1]], ncol=D_cols[[j+1]]))
                         }
                       }) %>% Reduce(rbind, .)
              }) %>% Reduce(cbind, .)

  D_theta <- lapply(0:(length(prm$t_a_vec)-1),
              function(k){
                lapply(0:(length(prm$t_a_vec)-1),
                function(j){
                  if(k==j){
                    return(rep(1, n_vec[[j+1]]))
                  } else {
                    return(rep(0, n_vec[[j+1]]))
                  }
                }) %>% Reduce(c, .)
              }) %>% Reduce(cbind, .)
  
  # D_theta <- sapply(1:max(prm$beta_1_track),
  #                   function(beta_1_val){
  #                     #for (k in 0:(length(prm$t_a_vec)-1))
  #                     sapply(0:(length(prm$t_a_vec)-1), beta_1_val = beta_1_val,
  #                            FUN = function(k, beta_1_val){
  #                              if (prm$beta_1_track[[k+1]] == beta_1_val) {
  #                                return(rep(1, n_vec[[k+1]]))
  #                              } else {
  #                                return(rep(0, n_vec[[k+1]]))
  #                              }
  #                            }) %>% Reduce(c, .) %>% matrix(ncol=1)
  #                   })
  
  
  fit <- c()
  tau <- c()
  for (k in 0:(length(prm$t_a_vec)-1)){
    df_temp <- df %>% filter( .$ti > prm$t_a_vec[[k+1]]) %>% 
      calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=k)
    fit <- c(fit, df_temp[[paste0('fit_', k)]])
    tau <- c(tau, df_temp$tau_k)
  }
  fit_list <- lapply(0:(length(prm$t_a_vec)-1),
                function(k){
                  df_temp <- df %>% filter( .$ti > prm$t_a_vec[[k+1]]) %>%
                    calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=k)
                  return(df_temp[[paste0('fit_', k)]])
                })
  
  Jbb_list <- lapply(1:length(design_matrices), function(k){
    sapply(1:dim(design_matrices[[k]])[[2]], function(i){
      sapply(1:dim(design_matrices[[k]])[[2]], function(j){
        fit_k <- fit_list[[k]]
        design_k <- design_matrices[[k]]
        return(sum(design_k[,i]*design_k[,j]*fit_k*(1-fit_k)))
      })
    })
  })
  
  Jbb_dims <- sapply(1:length(Jbb_list),
                function(k){
                  Jbb_part <- Jbb_list[[k]] %>% as.matrix()
                  return(dim(Jbb_part)[1])
                })
  
  Jbb <- lapply(1:length(Jbb_list),
          function(i){
            lapply(1:length(Jbb_list),
            function(j){
              if(j==i){
                return(Jbb_list[[i]])
              } else{
                return(matrix(rep(0, Jbb_dims[i]*Jbb_dims[j]),
                              nrow=Jbb_dims[i]))
              }
            }) %>% Reduce(cbind, .)
          }) %>% Reduce(rbind, .)
  
  
  Jtt <- lapply(0:(length(prm$t_a_vec) - 1),
          function(k){
            df_k <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=k)
            
            a_k <- as.name(paste0("a_", k))
            fit_k <- as.name(paste0("fit_", k))
            
            df_k <- df_k %>% mutate(theta_bit = {{fit_k}}*(1-{{fit_k}})*tau_k*tau_k)
            
            return(sum(df_k$theta_bit))
          }) %>% diag()
  
  Jtb <- lapply(1:dim(D_theta)[[2]], function(i){
          lapply(1:dim(D)[[2]], function(j){
      sum(D_theta[,i]*D[,j]*fit*(1-fit)*tau, na.rm=T)
    }) %>% Reduce(cbind, .)
  }) %>% Reduce(rbind, .)

  Jbb_inv <- solve(Jbb)
  if(dim(Jbb_inv)[1] == 1 && dim(Jbb_inv)[2] == 1){
    JTT <- Jtt - Jtb %*% t(Jtb) * as.numeric(Jbb_inv)
  } else {
    JTT <- Jtt - Jtb%*%solve(Jbb)%*%t(Jtb)
  }

  #Jacobian <- df %>% calculate_jacobian(prm=prm, psi_hat_vec=psi_hat_vec)
  Jacobian <- df %>% calculate_jacobian_large(prm=prm, psi_hat_vec=psi_hat_vec)
  J_inv <- MASS::ginv(Jacobian)
  (VCV <- J_inv %*% JTT %*% t(J_inv))
  #(VCV_alt <- J_inv_alt %*% JTT %*% t(J_inv_alt))

  return(VCV)
}





# ALT VARIANT






calculate_score_alt <- function(df, psi_hat_vec, prm){
  
  S_vec <- sapply((0:(length(prm$t_a_vec)-1)),
                  function(a_curr){
    df_k <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=a_curr)  
  
    a_k <- as.name(paste0("a_", a_curr))
    fit_k <- as.name(paste0("fit_", a_curr))
    
    df_k <- df_k %>% mutate(
                    S_inc = ifelse((is.na({{a_k}}) | is.na({{fit_k}}) | is.na(tau_k)), 0,
                                                     ({{a_k}} - {{fit_k}})*tau_k))
    
    return(df_k %>% dplyr::select(S_inc) %>% sum())
  })
  
  return(S_vec)
  
}


calculate_jacobian_alt <- function(df, prm, psi_hat_vec){
  jacobi_vec <- c()
  jacobi_dim <- length(psi_hat_vec)
  
  #for (row_counter in 1:jacobi_dim){
  jacobian <- lapply(1:jacobi_dim, 
  function(row_counter){
    
    #for(col_counter in 1:jacobi_dim) {
    lapply(1:jacobi_dim,
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
          
          df_k <- df %>% 
                    calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=a_curr) %>%
                    mutate(Si_temp = 0, dC_dpsi = 0)

          for (m in (a_curr:(length(prm$t_a_vec)-1))) {
            if (#(beta_track_row[[m+1]] == i) &&
                (beta_track_col[[m+1]] == j)) {
              
              df_k <- df_k %>% calculate_tau_rsd_m(prm=prm, psi_hat_vec=psi_hat_vec, m=m)
              
              a_m <- as.name(paste0("a_", m))
              #names(df_k)[names(df_k)==paste0("a_", m)] <- "a_m"
              
              df_k <- df_k %>% mutate(
                x_wrt = (1 - as.integer(d_wrt == "theta_x")) + x * as.integer(d_wrt == "theta_x"),
                Si_inc = ifelse(is.na({{a_m}}), 0, {{a_m}} * x_wrt * tau_rsd),
                Si_temp = Si_temp + Si_inc)
              
              #names(df_k)[names(df_k)=="a_m"] <- paste0("a_", m)
              
              if(prm$censor){
                df_k <- df_k %>% calculate_C_m(psi_hat_vec=psi_hat_vec, prm=prm, m=m)
                df_k <- df_k %>% mutate(
                          dC_dpsi = dC_dpsi + as.integer(psi_hat_vec[[j]] > 0)*C_m
                        )
              }
            }
          }
          
          df_k <- df_k %>% mutate(
            x_num = (1 - as.integer(d_num == "theta_x")) + x * as.integer(d_num == "theta_x"),
          )
          if(prm$censor) { 
            df_k <- df_k %>% mutate(
                      Si_censor = delta_k*dC_dpsi + (1-delta_k)*Si_temp
                      )
          } else {
            df_k <- df_k %>% mutate(
              Si_censor = Si_temp)
          }
          
          a_k <- as.name(paste0("a_", a_curr))
          fit_k <- as.name(paste0("fit_", a_curr))
          #names(df_k)[names(df_k)==paste0("a_", a_curr)] <- "a_k"
          #names(df_k)[names(df_k)==paste0("fit_", a_curr)] <- "fit_k"
          
          df_k <- df_k %>% mutate(denom_inc = ifelse((is.na({{a_k}}) | is.na({{fit_k}})), 0, 
                                                     ({{a_k}} - {{fit_k}})*x_num*Si_censor))
          
          denom <- denom + sum(with(df_k, denom_inc))
        }
      }
      
      return(denom)
    }) %>% Reduce(cbind, .)
  }) %>% Reduce(rbind, .)
  
  return(jacobian)
}




newton_raphson_grad_alt <- function(df, prm,
                                psi_start_vec=NA, 
                                tol=0.001, min_damp=0.1,
                                max_iter=50, 
                                psi_max = 10,
                                psi_max_vec = NA,
                                print_results=F){
  
  if(all(is.na(psi_start_vec))){
    psi_start_vec <- rep(0, max(prm$beta_1_track)+max(prm$beta_x_track))
  }
  
  if(all(is.na(psi_max_vec))){
    psi_max_vec = rep(psi_max, length(psi_start_vec))
  }
  
  psi_hat_vec <- psi_start_vec
  S_vec <- df %>% calculate_score_alt(psi_hat_vec=psi_hat_vec, prm=prm)
  psi_old_vec <- psi_hat_vec
  steps <- 1L
  fail_flag <- F
  
  while(((sum(abs(psi_hat_vec - psi_old_vec)) > tol) |
                steps == 1 |
                (sum(abs(psi_hat_vec)) < tol)) &
        (steps <= max_iter) &
        (max(abs(psi_hat_vec)) < psi_max)){
    
    if(print_results==T){
      cat("(", steps, ")", " --- ", psi_old_vec, "\n")
    } else {
      print(steps)
    }
    
    psi_old_vec <- psi_hat_vec
    
    (S_vec <- df %>% calculate_score_alt(psi_hat_vec=psi_hat_vec, prm=prm))
    # if(print_results){
    #   cat("\nS\n")
    #   print(S_vec)
    # }
    (J <- df %>% calculate_jacobian_alt(prm=prm, psi_hat_vec=psi_hat_vec))
    # if(print_results){
    #   cat("\nJ\n")
    #   print(J)
    # }
    damping_factor <- max(2^(-(steps - 1)/20), min_damp)
    (psi_hat_vec <- solve(t(J) %*% J) %*% t(J)  %*% S_vec *damping_factor + psi_hat_vec)
    # if(print_results){
    #   cat("\npsi_hat\n")
    #   print(psi_hat_vec)
    # }
    # if(print_results){cat("\n\n")}
    
    if(!all(abs(psi_hat_vec) < psi_max_vec)){
      psi_hat_vec <- psi_start_vec
    }
    
    
    # psi_hat_vec <- sapply(1:length(psi_hat_vec), function(d){
    #                   psi_bit <- psi_hat_vec[d]
    #                   if(psi_bit > -psi_max_vec[d] & psi_bit < psi_max_vec[d]){
    #                     return(psi_bit)
    #                   } else {
    #                     return(0)
    #                   }
    #                 })
    
    steps <- steps + 1
  }
  
  if ((max(abs(psi_hat_vec)) > psi_max) | (steps >= max_iter)) {
    fail_flag <- T
    psi_hat_vec <- rep(NA, length(prm$psi_star_vec))
  }
  
  return(list(psi_hat_vec, steps, fail_flag))
}


calculate_jacobi_piece <- function(df, prm, psi_hat_vec, d){

  row_counter <- d
  col_counter <- d
  
    
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
      
      df_k <- df %>% 
                  calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=a_curr) %>%
                  mutate(Si_temp = 0, 
                         dC_dpsi = 0)
      fit_k <- as.name(paste0("fit_", a_curr))
      for (m in (a_curr:(length(prm$t_a_vec)-1))) {
        if ((beta_track_row[[m+1]] == i) &&
            (beta_track_col[[m+1]] == j)) {
          
          df_k <- df_k %>% 
                    calculate_tau_rsd_m(prm=prm, psi_hat_vec=psi_hat_vec, m=m)
          
          a_m <- as.name(paste0("a_", m))
          a_mprev <- as.name(paste0("a_", m-1))
          #names(df_k)[names(df_m)==] <- "a_m"
          
          #df_k <- df_k %>%
          #  dplyr::select(id,{{a_mprev}}, {{a_m}}, x, tau_rsd, Si_temp, dC_dpsi, {{fit_k}})
          
          df_k <- df_k %>% mutate(
            x_wrt = (1 - as.integer(d_wrt == "theta_x")) + x*as.integer(d_wrt == "theta_x"),
            Si_inc = ifelse(is.na({{a_m}}), 0, {{a_m}}*x_wrt*tau_rsd),
            Si_temp = Si_temp + Si_inc)
          
          if(prm$censor){
            df_k <- df_k %>% calculate_C_m(psi_hat_vec=psi_hat_vec, prm=prm, m=m)
            df_k <- df_k %>% mutate(
                      dC_dpsi = dC_dpsi + as.integer(psi_hat_vec[[j]] > 0)*C_m
                    )
          }
        }
      }
      
      df_k <- df_k %>% mutate(
        x_num = (1 - as.integer(d_num == "theta_x")) + x * as.integer(d_num == "theta_x"))
      
      if(prm$censor) { 
        df_k <- df_k %>% mutate(
                    Si_censor = delta_k*dC_dpsi + (1-delta_k)*Si_temp)

      } else {
        df_k <- df_k %>% mutate(
          Si_censor = Si_temp)
      }
      a_k <- as.name(paste0("a_", a_curr))
      fit_k <- as.name(paste0("fit_", a_curr))
      # names(df_k)[names(df_k)==paste0("a_", a_curr)] <- "a_k"
      # names(df_k)[names(df_k)==paste0("fit_", a_curr)] <- "fit_k"
      
      df_k <- df_k %>% mutate(denom_inc = ifelse((is.na({{a_k}}) | is.na({{fit_k}})), 0, 
                                                 ({{a_k}} - {{fit_k}})*x_num*Si_censor))
      
      denom <- denom + sum(with(df_k, denom_inc))
    }
  }
                                
  return(denom)

}




calculate_variance_alt <- function(df, prm, psi_hat_vec,# trt_models,
                               print_results=T){
  

  n_vec <- sapply(prm$t_a_vec, 
                  function(t_a) nrow(df %>% filter(.$ti > t_a)))
  
  design_matrices <- lapply(prm$trt_models,
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
  
  D <- lapply(0:(length(prm$t_a_vec)-1),
              function(j) {
                lapply(0:(length(prm$t_a_vec)-1),
                       function(i) {
                         if (i == j) {
                           return(design_matrices[[j+1]])
                         } else {
                           return(matrix(0, nrow=D_rows[[i+1]], ncol=D_cols[[j+1]]))
                         }
                       }) %>% Reduce(rbind, .)
              }) %>% Reduce(cbind, .)


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
  
  
  fit <- c()
  tau <- c()
  for (k in 0:(length(prm$t_a_vec)-1)){
    df_temp <- df %>% filter( .$ti > prm$t_a_vec[[k+1]]) %>% 
      calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=k)
    fit <- c(fit, df_temp[[paste0('fit_', k)]])
    tau <- c(tau, df_temp$tau_k)
  }
  fit_list <- lapply(0:(length(prm$t_a_vec)-1),
                function(k){
                  df_temp <- df %>% filter( .$ti > prm$t_a_vec[[k+1]]) %>%
                    calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=k)
                  return(df_temp[[paste0('fit_', k)]])
                })
  
  Jbb_list <- lapply(1:length(design_matrices), function(k){
    sapply(1:dim(design_matrices[[k]])[[2]], function(i){
      sapply(1:dim(design_matrices[[k]])[[2]], function(j){
        fit_k <- fit_list[[k]]
        design_k <- design_matrices[[k]]
        return(sum(design_k[,i]*design_k[,j]*fit_k*(1-fit_k)))
      })
    })
  })
  
  Jbb_dims <- sapply(1:length(Jbb_list),
                function(k){
                  Jbb_part <- Jbb_list[[k]] %>% as.matrix()
                  return(dim(Jbb_part)[1])
                })
  
  Jbb <- lapply(1:length(Jbb_list),
          function(i){
            lapply(1:length(Jbb_list),
            function(j){
              if(j==i){
                return(Jbb_list[[i]])
              } else{
                return(matrix(rep(0, Jbb_dims[i]*Jbb_dims[j]),
                              nrow=Jbb_dims[i]))
              }
            }) %>% Reduce(cbind, .)
          }) %>% Reduce(rbind, .)
  
  Jtt <- t(sapply(1:dim(D_theta)[[2]], function(i){
    sapply(1:dim(D_theta)[[2]], function(j){
      sum(D_theta[,i]*D_theta[,j]*fit*(1-fit)*tau*tau, na.rm=T)
    })
  }))
  
  Jtb <- lapply(1:dim(D_theta)[[2]], function(i){
          lapply(1:dim(D)[[2]], function(j){
      sum(D_theta[,i]*D[,j]*fit*(1-fit)*tau, na.rm=T)
    }) %>% Reduce(cbind, .)
  }) %>% Reduce(rbind, .)

  Jbb_inv <- solve(Jbb)
  if(dim(Jbb_inv)[1] == 1 && dim(Jbb_inv)[2] == 1){
    JTT <- Jtt - Jtb %*% t(Jtb) * as.numeric(Jbb_inv)
  } else {
    JTT <- Jtt - Jtb%*%solve(Jbb)%*%t(Jtb)
  }

  #Jacobian <- df %>% calculate_jacobian(prm=prm, psi_hat_vec=psi_hat_vec)
  Jacobian <- df %>% calculate_jacobian_alt(prm=prm, psi_hat_vec=psi_hat_vec)
  J_inv <- solve(Jacobian)
  (VCV <- J_inv %*% JTT %*% t(J_inv))
  #(VCV_alt <- J_inv_alt %*% JTT %*% t(J_inv_alt))

  return(VCV)
}
