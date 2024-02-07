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

        mdl_formula <- paste0("a_", k, " ~ ", 
                              paste0(prm$trt_mod_list[[k+1]], collapse = " + "))
        
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

calculate_tau_k <- function(df, psi_hat_vec, prm,  a_k=0, ...){

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
      tau_k = ifelse(is.na({{a_curr}}), tau_k, tau_k + ti_temp*exp(g_psi*{{a_curr}})),
      ti_rsd = ti_rsd - ti_temp,
      
    )
    
    if (prm$censor) {
      t_next <- ifelse(i < length(prm$t_a_vec), 
                       prm$t_a_vec[[i+1]], 
                       prm$censor_max)
      df <- df %>% mutate(
        C_psi = ifelse(g_psi >= 0, 1, exp(g_psi)),
        C_k = C_k + pmin(t_next - t_curr, C_rsd)*C_psi,
        C_rsd = C_rsd - pmin(t_next - t_curr, C_rsd)
      )
    }

    
  }
  if (prm$censor) {
    df <- df %>% mutate(delta_k = ifelse(tau_k < C_k, 0, 1),
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
    
    a_curr <- as.name(paste0("a_", m))
    
    df <- df %>% mutate(
      g_psi =  beta_1_now + x*beta_x_now,
      ti_temp = pmax(0, pmin(t_next, ti) - t_curr),
      tau_rsd =  ifelse(is.na({{a_curr}}), 0, ti_temp*exp(g_psi*{{a_curr}}))
    )

    return(df %>% dplyr::select(-c(g_psi, ti_temp)))
}

calculate_C_m <- function(df, psi_hat_vec, prm, m=0){

  t_now <- prm$t_a_vec[[m+1]]
  t_next <- ifelse(m+1 < length(prm$t_a_vec), 
                   prm$t_a_vec[[m+2]], 
                   prm$censor_date)
  
  beta_1_now <- ifelse(prm$beta_1_track[[m+1]]==0, 0,
                       psi_hat_vec[[prm$beta_1_track[[m+1]]]])
  beta_x_now <- ifelse(prm$beta_x_track[[m+1]]==0, 0,
                       psi_hat_vec[[prm$beta_x_track[[m+1]]+max(prm$beta_1_track)]])
  
  df <- df %>% mutate(C_check = ifelse(C_i < t_next, pmax(C_i - t_now, 0), t_next - t_now),
                      C_m = exp((beta_1_now*(beta_1_now<0) + beta_x_now*(beta_x_now<0)))*C_check)
  
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

          df_temp <- df_temp %>% mutate(
                    S_inc = ifelse((is.na({{a_k}}) | is.na({{fit_k}}) | is.na(tau_k)), 0,
                                                     ({{a_k}} - {{fit_k}})*x*tau_k))
        
          S_curr <- S_curr + with(df_temp, sum(S_inc))          
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
  
  jacobian <- lapply(1:jacobi_dim, 
  function(row_counter){
    
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
            if ((beta_track_col[[m+1]] == j)) {
              
              df_k <- df_k %>% calculate_tau_rsd_m(prm=prm, psi_hat_vec=psi_hat_vec, m=m)
              
              a_m <- as.name(paste0("a_", m))

              df_k <- df_k %>% mutate(
                x_wrt = (1 - as.integer(d_wrt == "theta_x")) + x * as.integer(d_wrt == "theta_x"),
                Si_inc = ifelse(is.na({{a_m}}), 0, {{a_m}} * x_wrt * tau_rsd),
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

    (J <- df %>% calculate_jacobian(prm=prm, psi_hat_vec=psi_hat_vec))

    damping_factor <- max(2^(-(steps - 1)/20), min_damp)
    (psi_hat_vec <- psi_hat_vec - solve(J, S_vec)*damping_factor)

    
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






calculate_variance <- function(df, prm, psi_hat_vec,
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
      sapply(1:max(prm$beta_x_track), 
             FUN = function(beta_x_val){
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

  Jacobian <- df %>% calculate_jacobian(prm=prm, psi_hat_vec=psi_hat_vec)
  J_inv <- solve(Jacobian)
  (VCV <- J_inv %*% JTT %*% t(J_inv))

  return(VCV)
}
