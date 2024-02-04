library(MASS, include.only=c("ginv"))
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


prm_times <- function(prm){
 
  prm$t_vec <- c(prm$t_a_vec, prm$t_psi_vec) %>% sort() %>% unique()
  
  
  psi_tracker <- NULL
  ak_tracker <- NULL
  for(time in prm$t_vec){
    psi_tracker <- c(psi_tracker,
                   (prm$t_psi_vec <= time) %>% as.integer() %>% sum())
    ak_tracker <- c(ak_tracker,
                    (prm$t_a_vec <= time) %>% as.integer() %>% sum() - 1)
  }
  
  prm$psi_tracker <- psi_tracker
  prm$ak_tracker <- ak_tracker
  
  return(prm)
  
}

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

calculate_tau_k <- function(df, psi_hat_vec, prm,  a_k=0, ...){
  #by default calculates tau(0) aka T0 from trial start
  
  C_psi_vec <- calculate_C_psi_vec(psi_hat_vec=psi_hat_vec, prm=prm)
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
  
  #time_vec <- prm$t_vec[prm$t_vec >= t_a]
  
  for (m in 0:(length(prm$t_vec)-1)){
    t_curr <- prm$t_vec[[m+1]]
    if(t_curr >= t_a){
      t_next <- ifelse(m < (length(prm$t_vec)-1), prm$t_vec[[m+2]], Inf)
      
      psi_now <- psi_hat_vec[[prm$psi_tracker[[m+1]]]]
      a_k <- as.name(paste0("a_", prm$ak_tracker[[m+1]]))
      
      df <- df %>% mutate(
        ti_temp = pmin(t_next - t_curr, ti_rsd),
        tau_k = ifelse(is.na({{a_k}}), tau_k, tau_k + ti_temp*exp(-psi_now*{{a_k}})),
        ti_rsd = ti_rsd - ti_temp
      )
      
      if (prm$censor) {
        t_next <- ifelse(m < (length(prm$t_vec)-1),
                         prm$t_vec[[m+2]],
                         prm$censor_max)
        df <- df %>% mutate(
          C_psi = exp(-C_psi_vec[m+1]),
          C_k = C_k + pmin(t_next - t_curr, C_rsd)*C_psi,
          C_rsd = C_rsd - pmin(t_next - t_curr, C_rsd)
        )
      }
    }
  }

  if (prm$censor) {
    df <- df %>% mutate(#tau_k_un = tau_k,
                        delta_k = ifelse(tau_k < C_k, 0, 1),
                        tau_k = pmin(tau_k, C_k))
    return(dplyr::select(df, -c(ti_temp, ti_rsd, C_rsd, C_k, C_psi)))
  } else {
    return(dplyr::select(df, -c(ti_temp, ti_rsd)))
  }
}

calculate_tau_rsd_m <- function(df, psi_hat_vec, prm,
                                m=0, ...){
  
    t_curr <- prm$t_vec[[m+1]]
    t_next <- ifelse(m+1 < length(prm$t_vec), prm$t_vec[[m+2]] , Inf)
    psi_now <- psi_hat_vec[[prm$psi_tracker[[m+1]]]]
    
    a_k <- as.name(paste0("a_", prm$ak_tracker[[m+1]]))
    
    df <- df %>% mutate(
      ti_temp = pmax(0, pmin(t_next, ti) - t_curr),
      tau_rsd =  ifelse(is.na({{a_k}}), 0, ti_temp*exp(-psi_now*{{a_k}}))
    )
    
    return(df %>% dplyr::select(-c(ti_temp)))
}

calculate_C_psi_vec <- function(psi_hat_vec, prm){
  psi_now_vec <- sapply(0:(length(prm$t_vec)-1),
                  function(m){
                    return(psi_hat_vec[prm$psi_tracker[m+1]])
                  })
  m_time <- sapply(0:(length(prm$t_vec)-1),
              function(m){
                t_now <- prm$t_vec[m+1]
                t_next <- ifelse(m < length(prm$t_vec)-1, 
                                   prm$t_vec[[m+2]], 
                                   prm$censor_date)
                return(t_next - t_now)
              })
  
  C_psi <- lapply(0:(length(prm$t_a_vec)-1),
            function(k){
              is_k <- (prm$ak_tracker == k)
              scaled_C <- (exp(-psi_now_vec[is_k])*m_time[is_k]) %>% sum()
              unscaled_C <- m_time[is_k] %>% sum()
              if(scaled_C < unscaled_C){
                return(psi_now_vec[is_k])
              } else {
                return(rep(0, length(psi_now_vec[is_k])))
              }
          }) %>% Reduce(c, .)

  return(C_psi)
}


calculate_C_m <- function(df, psi_hat_vec, prm, m=0){
  
  C_psi_vec <- calculate_C_psi_vec(psi_hat_vec=psi_hat_vec, prm=prm)
  
  t_now <- prm$t_vec[[m+1]]
  t_next <- ifelse(m+1 < length(prm$t_a_vec), 
                   prm$t_a_vec[[m+2]], 
                   prm$censor_date)
  
  psi_now <- psi_hat_vec[[prm$psi_tracker[[m+1]]]]
  
  df <- df %>% mutate(C_check = ifelse(C_i < t_next, pmax(C_i - t_now, 0), t_next - t_now),
                      C_m = exp(-C_psi_vec[m+1])*C_check)
  
  
  return(dplyr::select(df, -c(C_check)))
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
                lapply(1:max(prm$psi_tracker),
                function(psi_val){
                  df_k <- df %>% calculate_tau_k(prm=prm, psi_hat_vec=psi_hat_vec, a_k=a_curr) %>%
                        mutate(Si_temp = 0, dC_dpsi = 0)
                  
                  for(m in (0:(length(prm$t_vec)-1))){
                    if (prm$psi_tracker[[m+1]] == psi_val){
                      df_k <- df_k %>% calculate_tau_rsd_m(prm=prm, psi_hat_vec=psi_hat_vec, m=m)

                      a_m <- as.name(paste0("a_", prm$ak_tracker[[m+1]]))
                      df_k <- df_k %>% mutate(
                        Si_inc = ifelse(is.na({{a_m}}), 0, {{a_m}} * tau_rsd),
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
    psi_start_vec <- rep(0, max(prm$psi_tracker))
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



