# Testing GIT
# Testing GIT 

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


### DATASET CREATION / PRELIMINARY DATA WRANGLING ###



get_time_df <- function(prm, ...){
  
  t_a_vec <- prm$t_a_vec
  t_psi_vec <- prm$t_psi_vec
  xi_vec <- prm$xi_vec
  
  times_vec <- sort(unique(c(t_psi_vec, 
                              as.vector(sapply(xi_vec, function(xi){t_a_vec + xi})))))
  time_df <- tibble(t=numeric(), t_next=numeric(),
                    a_val=integer(), psi_counter=integer())
  for(i in 1:length(times_vec)){
    time <- times_vec[[i]]
    if(i < length(times_vec)){
      time_next <- times_vec[[i+1]]
    } else {
      time_next <- Inf
    }
    a_val       <- as.integer(length(t_a_vec)-1 - sum(t_a_vec > time))
    psi_counter <- as.integer(length(t_psi_vec)-1 - sum(t_psi_vec > time))
    
    time_df <- time_df %>% add_row(t=time, t_next=time_next,
                                   a_val=a_val,
                                   psi_counter=psi_counter)
  }
  return(time_df)
}

create_m_k_map <- function(t_df){
  a_val_curr <- 0
  m_k_map <- c(0)
  for (i in 1:nrow(t_df)){
    time_row <- t_df[i,]
    if (time_row$a_val > a_val_curr){
      m_k_map <- c(m_k_map, i-1)
      a_val_curr <- time_row$a_val
    }
  }
  return(m_k_map)
}

create_sample <- function(prm, psi_star_list){
  
  t_df <- prm$t_df
  xi_vec <- prm$xi_vec
  psi_star_CT <- prm$psi_star_CT
  psi_star_TOT <- prm$psi_star_TOT
  n <- prm$n
  t0_min <- prm$t0_min
  t0_max <- prm$t0_max
  
  # creates treatment variables
  a_0 <- c(rep(0, floor(n/2)), rep(1, n - floor(n/2)))
  df <- tibble(a_0)
  if (max(t_df$a_val) > 0) {
    for(i in c(1:max(t_df$a_val))){
      a_i <- sample(a_0, size=n, replace=F)
      df <- df %>% add_column(a_i)
      names(df)[names(df)=="a_i"] <- paste("a_",i, sep="")
    }
  }
  
  # generates T0 times
  df <- df %>% add_column(t0  = runif(n, min=t0_min, max=t0_max )) %>%
    mutate(ti = 0, t0_rsd = t0)
  
  # pre-fills TOT variables
  df <- df %>% add_column(w_TOT = rep(0,n))
  for(m in 0:(nrow(t_df)-1)){
    df <- df %>% add_column(w_m = rep(0,n))
    names(df)[names(df)=="w_m"] <- paste0("w_m",m)
  }
  
  
  # Creates a variable "w_m{}" for each sub-period m
  # that records the index of the highest TOT psi activated
  # by the individual
  for(m in 0:(nrow(t_df)-1)){
    time_row <- t_df[m+1,]
    t_curr <- time_row$t
    a_curr <- time_row$a_val
    psi_curr <- time_row$psi_counter
    
    psi_val_CT <- psi_star_CT[[psi_curr+1]]
    
    if(m == (nrow(t_df)-1)){
      t_max <- Inf
    } else {
      t_max <- t_df[m+1,]$t_next
    }
    
    names(df)[names(df)==paste0("w_m",m)] <- "w_m_curr"
    names(df)[names(df)==paste0("a_", a_curr)] <- "a_vec"
    
    df <- df %>% mutate(
      w_m_curr = sapply(w_TOT, function(q){sum(q >= xi_vec) - 1}),
      g_psi = psi_val_CT + rowSums(sapply(0:(length(psi_star_TOT)-1),
                                          function(q){
                                            psi_star_TOT[[q+1]]*(w_m_curr >= q)
                                          })),
      temp = pmin((t_max - t_curr)*exp(-g_psi*a_vec), t0_rsd),
      ti = ti + temp*exp(g_psi*a_vec),
      t0_rsd = t0_rsd - temp,
      w_TOT = w_TOT + a_vec*(t_max - t_curr)
    )
    
    names(df)[names(df)=="w_m_curr"] <- paste0("w_m", m)
    names(df)[names(df)=="a_vec"] <- paste0("a_", a_curr) 
  }
  
  return(select(df, -c(t0_rsd, temp, w_TOT, g_psi)))
}

fit_treatment_models <- function(df, prm){
  
  t_a_vec <- prm$t_a_vec
  #t_psi_vec <- prm$t_psi_vec
  #xi_vec <- prm$xi_vec
  
  df <- df %>% mutate(fit_0 = glm(a_0 ~ 1, family=binomial)$fitted.values)
  
  if (length(t_a_vec)-1 > 0) {
    for(i in 1:(length(t_a_vec)-1)) {
      #     mdl_formula <- paste0("a_", i, "~",
      #                                      paste0("a_", 0:(i-1), collapse="+"))
      #     if (length(xi_vec) > 1) {
      #       mdl_formula <- paste0(mdl_formula, "+", paste0("w", 0:(length(xi_vec)-1), collapse="+"))
      #       
      #       subperiod_index <- sum((t_df$a_val < i))
      #       names(df)[names(df) == paste0("w_m", subperiod_index)] <- "w_compare"
      #       for (j in 0:(length(xi_vec)-1)){
      #         df <- df %>% mutate(w_curr = as.integer((j <= w_compare)))
      #         names(df)[names(df) == "w_curr"] <- paste0("w", j)
      #       }
      #       names(df)[names(df) == "w_compare"] <- paste0("w_m", subperiod_index)
      #     }
      
      mdl_formula <- paste0("a_", i, "~",
                            paste0("a_", 0:(i-1), collapse="+"))
      model_temp <- with(df %>% filter(.$ti > t_a_vec[[i+1]]), 
                         glm(  as.formula(mdl_formula),
                               family=binomial))
      df$fit_temp <- predict(model_temp, newdata = df, type="response")
      names(df)[names(df) == "fit_temp"] <- paste0("fit_", i)
    }
  }
  
  
  return(df)
}



### AFT CALCULATIONS ###
calculate_tau_k <- function(df, psi_hat_list, prm,  a_k=0, ...){
  #by default calculates tau(0) aka T0 from trial start
  
  t_df <- prm$t_df
  t_a_vec <- prm$t_a_vec
  
  psi_hat_CT <- psi_hat_list[[1]]
  psi_hat_TOT <- psi_hat_list[[2]]
  
  t_a <- t_a_vec[[a_k+1]]
  
  t_df_k <- t_df %>%
    filter(.$t >= t_a)
  
  df <- df %>% filter(.$ti > t_a)%>% 
    mutate(ti_rsd = ti - t_a,
           tau_k = 0)
  
  for(i in 1:nrow(t_df_k)){
    time_row <- t_df_k[i,]
    t_curr <- time_row$t
    t_max <- time_row$t_next
    a_curr <- time_row$a_val
    psi_curr <- time_row$psi_counter
    
    psi_val_CT <- psi_hat_CT[[psi_curr+1]]
    
    m_val <- sum((t_df$t <= t_curr))-1
    
    names(df)[names(df)==paste0("w_m", m_val)] <- "w_m_curr"
    names(df)[names(df)==paste0("a_", a_curr)] <- "a_k"
    
    df <- df %>% mutate(
      g_psi = psi_val_CT + 
        rowSums(sapply(0:(length(psi_hat_TOT)-1),
                  function(q){
                    psi_hat_TOT[[q+1]]*(w_m_curr >= q)
                  })),
      temp = pmin(t_max - t_curr, ti_rsd),
      tau_k = tau_k + temp*exp(-g_psi*a_k),
      ti_rsd = ti_rsd - temp
    )
    
    names(df)[names(df)=="w_m_curr"] <- paste0("w_m", m_val)
    names(df)[names(df)=="a_k"] <- paste0("a_", a_curr)
  }
  #return(df)
  return(select(df, -c(temp, ti_rsd, g_psi)))
}

calculate_tau_rsd <- function(df, psi_hat_list,
                              prm, 
                              a_k=0, psi_k=0, ...){
  #by default calculates tau(0) aka T0 from trial start
  
  t_df <- prm$t_df
  
  psi_hat_CT <- psi_hat_list[[1]]
  psi_hat_TOT <- psi_hat_list[[2]]
  
  time_row_df <- t_df %>%
    filter(.$a_val == a_k & .$psi_counter == psi_k)
  
  t_curr <- time_row_df$t[[1]]
  t_max <- time_row_df$t_next[[1]]
  
  psi_val_CT <- psi_hat_CT[[psi_curr+1]]
  m_val <- sum((t_df$t <= t_curr))-1
  
  names(df)[names(df)==paste0("w_m", m_val)] <- "w_m_curr"
  names(df)[names(df)==paste0("a_", a_k)] <- "a_k"
  
  df <- df %>% mutate(
    g_psi = psi_val_CT + 
      rowSums(sapply(0:(length(psi_hat_TOT)-1),
                     function(q){
                       psi_hat_TOT[[q+1]]*(w_m_curr >= q)
                     })),
    temp = pmax(0, pmin(t_max, ti) - t_curr),
    tau_rsd = temp*exp(-g_psi*a_k)
  )
  
  names(df)[names(df)=="w_m_curr"] <- paste0("w_m", m_val)
  names(df)[names(df)=="a_k"] <- paste0("a_", a_curr)
  
  #return(df)
  return(select(df, -c(temp, g_psi)))
}

calculate_tau_rsd_m <- function(df, psi_hat_list,
                                prm,
                                m=0, ...){
  #by default calculates tau(0) aka T0 from trial start
  t_df <- prm$t_df
  
  psi_hat_CT <- psi_hat_list[[1]]
  psi_hat_TOT <- psi_hat_list[[2]]
  
  time_row_df <- t_df[m+1,]
  
  t_curr <- time_row_df$t
  t_max <- time_row_df$t_next
  a_k <- time_row_df$a_val
  psi_curr <- time_row_df$psi_counter
  
  psi_val_CT <- psi_hat_CT[[psi_curr+1]]
  
  names(df)[names(df)==paste0("w_m", m)] <- "w_m_curr"
  names(df)[names(df)==paste0("a_", a_k)] <- "a_k"
  
  df <- df %>% mutate(
    g_psi = psi_val_CT + 
      rowSums(sapply(0:(length(psi_hat_TOT)-1),
                     function(q){
                       psi_hat_TOT[[q+1]]*(w_m_curr >= q)
                     })),
    temp = pmax(0, pmin(t_max, ti) - t_curr),
    tau_rsd = temp*exp(-g_psi*a_k)
  )
  
  names(df)[names(df)=="w_m_curr"] <- paste0("w_m", m)
  names(df)[names(df)=="a_k"] <- paste0("a_", a_k)
  
  #return(df)
  return(select(df, -c(temp, g_psi)))
}

calculate_dtau_dpsi <- function(df, psi_hat, 
                                prm,
                                a_k=0, psi_k=0, ...){
  #by default calculates tau(0) aka T0 from trial start
  
  t_df <- prm$t_df
  
  t_df_tmp <- t_df %>%
    filter(.$a_val >= a_k & .$psi_counter == psi_k)
  psi_val <- psi_hat[[psi_k+1]]
  
  df <- df %>% mutate(dtau_dpsi = 0)
  for (l in 1:nrow(t_df_tmp)) {
    t_min_l <- t_df_tmp$t[[l]]
    t_max_l <- t_df_tmp$t_next[[l]]
    a_curr <- t_df_tmp$a_val[[l]]
    
    names(df)[names(df)==paste0("a_", a_curr)] <- "a_dash"
    df <- df %>% mutate(dtau_dpsi = dtau_dpsi +  
                          pmax(0, pmin(ti, t_max_l) - t_min_l)*
                          exp(-psi_val*a_dash)*a_dash)
    names(df)[names(df)=="a_dash"] <- paste0("a_", a_curr)
  }
  
  return(df)
  #return(select(df, -c(temp)))
}


### LR CALCULATIONS ###
calculate_score <- function(df, prm, psi_hat_list){
  
  t_psi_vec <- prm$t_psi_vec
  t_df <- prm$t_df
  xi_vec <- prm$xi_vec
  
  S_CT <- sapply(0:(length(t_psi_vec)-1), function(psi_choose){
    num <- 0
    a_check <- -1
    for(i in 1:nrow(t_df)){
      if(t_df$psi_counter[[i]] == psi_choose){
        a_curr <- t_df$a_val[[i]]
        if (a_check < a_curr){
          df_temp <- calculate_tau_k(df, psi_hat_list=psi_hat_list, 
                                     prm=prm, a_k=a_curr)
          
          names(df_temp)[names(df_temp)==paste0("a_", a_curr)] <- "a_k"
          names(df_temp)[names(df_temp)==paste0("fit_", a_curr)] <- "fit_k"
          
          num <- num + with(df_temp, sum((a_k-fit_k)*tau_k))
        }
        a_check <- a_curr
      }
    }
    return(num)
  }) 
  
  #TOT varying component
  if (length(xi_vec) > 1) {
    S_TOT <-sapply(1:(length(xi_vec)-1), function(psi_TOT_choose){
      num <- 0
      a_check <- -1
      for(m in 0:(nrow(t_df)-1)){
        a_curr <- t_df$a_val[[m+1]]
        if (a_check < a_curr){
          
          df_temp <- df %>% calculate_tau_k(psi_hat_list=psi_hat_list,
                                            prm=prm, a_k=a_curr)
          
          names(df_temp)[names(df_temp)==paste0("w_m", m)] <- "w_m_k"
          df_temp <- df_temp %>% mutate(
            w_nu = as.integer(w_m_k >= psi_TOT_choose))
          names(df_temp)[names(df_temp)==paste0("a_", a_curr)] <- "a_k"
          names(df_temp)[names(df_temp)==paste0("fit_", a_curr)] <- "fit_k"
          
          num <- num + with(df_temp, sum((a_k-fit_k)*w_nu*tau_k))
        }
        a_check <- a_curr
      }
      return(num)
    })
  } else {
    S_TOT <- c()
  }
  
  S_vec <- c(S_CT, S_TOT)
  
  return(S_vec)
  
  
}

calculate_jacobian <- function(df, prm, psi_hat_list){
  
  t_df <- prm$t_df
  
  psi_hat_CT <- psi_hat_list[[1]]
  psi_hat_TOT <- psi_hat_list[[2]]
  
  jacobi_vec <- c()
  
  jacobi_dim <- length(psi_hat_CT) + (length(psi_hat_TOT)-1)
  
  for (row_counter in 1:jacobi_dim){
    for (col_counter in 1:jacobi_dim){
      
      denom <- 0
      
      if (col_counter <= length(psi_hat_CT)){
        d_wrt <- "CT"
        j <- col_counter - 1
      } else {
        d_wrt <- "TOT"
        j <- col_counter - length(psi_hat_CT)
      }
      
      #dCT/dCT
      if (row_counter <= length(psi_hat_CT)){
        # SCORE IS CT COMPONENT (a_k - pi_k)*T_0
        i <- row_counter - 1

        psi_choose <- i
        t_df_psi <- t_df %>% filter( .$psi_counter == psi_choose)
        t_df_a <- t_df %>% filter(.$a_val %in% t_df_psi$a_val)
        
        for(a_curr in t_df_a$a_val){
          for(m in 0:(nrow(t_df)-1)){
            a_dash <- t_df$a_val[[m+1]]
            psi_dash <- t_df$psi_counter[[m+1]]
            if ((a_dash %in% t_df_psi$a_val) 
                & (a_dash >= a_curr) 
                & ((psi_dash == j) | (d_wrt == "TOT"))
                ){
              df_temp <- df %>% calculate_tau_rsd_m(prm=prm,
                                             psi_hat_list=psi_hat_list,
                                             m=m)
              
              names(df_temp)[names(df_temp)==paste0("w_m", m)] <- "w_m_curr"
              df_temp <- df_temp %>% mutate(w_nu_j = as.integer(w_m_curr >= j))
              
              names(df_temp)[names(df_temp)==paste0("a_", a_curr)] <- "a_k"
              names(df_temp)[names(df_temp)==paste0("a_", a_dash)] <- "a_dish"
              names(df_temp)[names(df_temp)==paste0("fit_", a_curr)] <- "fit_k"
              
              if (d_wrt == "CT") {
                if(a_dash == a_curr){
                  denom <- denom + with(df_temp, sum((a_k-fit_k)*a_k*tau_rsd))
                } else {
                  denom <- denom + with(df_temp, sum((a_k-fit_k)*a_dish*tau_rsd))
                }
              } else {
                #d_wrt == "TOT"
                if(a_dash == a_curr){
                  denom <- denom + with(df_temp, sum((a_k-fit_k)*a_k*w_nu_j*tau_rsd))
                } else {
                  denom <- denom + with(df_temp, sum((a_k-fit_k)*a_dish*w_nu_j*tau_rsd))
                }
              }
            }
          }
        }
        
      }
      
      else {
        
        i <- row_counter - length(psi_hat_CT)
        for(a_curr in t_df$a_val){
          for(m in 0:(nrow(t_df)-1)){
            psi_dash <- t_df$psi_counter[[m+1]]
            a_dash <- t_df$a_val[[m+1]]
            
            if((a_dash >= a_curr)
               & ((psi_dash == j) | (d_wrt == "TOT"))
            ){
              
              df_temp <- df %>% calculate_tau_rsd_m(prm=prm,
                                             psi_hat_list=psi_hat_list,
                                             m=m)
              
              m_k <- sum(t_df$a_val < a_curr)
              names(df_temp)[names(df_temp)==paste0("w_m", m_k)] <- "w_m_k"
              df_temp <- df_temp %>% mutate(w_nu_i = as.integer(w_m_k >= i))
              names(df_temp)[names(df_temp)=="w_m_k"] <- paste0("w_m", m_k)
              
              names(df_temp)[names(df_temp)==paste0("w_m", m)] <- "w_m_curr"
              df_temp <- df_temp %>% mutate(w_nu_j = as.integer(w_m_curr >= j))
              
              names(df_temp)[names(df_temp)==paste0("a_", a_curr)] <- "a_k"
              names(df_temp)[names(df_temp)==paste0("a_", a_dash)] <- "a_dish"
              names(df_temp)[names(df_temp)==paste0("fit_", a_curr)] <- "fit_k"
              
              if (d_wrt == "CT") {
                if(a_dash == a_curr){
                  denom <- denom + with(df_temp, sum((a_k-fit_k)*a_k*w_nu_i*tau_rsd))
                } else {
                  denom <- denom + with(df_temp, sum((a_k-fit_k)*a_dish*w_nu_i*tau_rsd))
                }
              } else {
                #d_wrt == "TOT"
                if(a_dash == a_curr){
                  denom <- denom + with(df_temp, 
                                        sum((a_k-fit_k)*a_k*w_nu_i*w_nu_j*tau_rsd))
                } else {
                  denom <- denom + with(df_temp, sum((a_k-fit_k)*a_dish*w_nu_i*w_nu_j*tau_rsd))
                }
              }
              
              

            }
          }
        }
      }
      
      jacobi_vec <- c(jacobi_vec, denom)
    }
  }
  
  
  return(matrix(jacobi_vec, byrow=T, nrow= (length(psi_hat_CT) + length(psi_hat_TOT) - 1)))
  
}


calculate_hessian <- function(df, prm, psi_hat_list){
  
  t_df <- prm$t_df
  
  psi_hat_CT <- psi_hat_list[[1]]
  psi_hat_TOT <- psi_hat_list[[2]]
  
  hessian_vec <- c()
  
  hessian_dim <- length(psi_hat_CT) + (length(psi_hat_TOT)-1)
  
  for (row_counter in 1:hessian_dim){
    for (col_counter in 1:hessian_dim){
      
      denom <- 0
      
      if (row_counter <= length(psi_hat_CT)){
        d_fun <- "CT"
        i <- row_counter - 1
      } else {
        d_fun <- "TOT"
        i <- row_counter - length(psi_hat_CT)
      }
      
      if (col_counter <= length(psi_hat_CT)){
        d_wrt <- "CT"
        j <- col_counter - 1
      } else {
        d_wrt <- "TOT"
        j <- col_counter - length(psi_hat_CT)
      }
      
      for(a_curr in 0:(length(prm$t_a_vec)-1)){
        df_temp <- df %>% calculate_tau_k(psi_hat_list=psi_hat_list,
                                          prm=prm, a_k=a_curr)
        names(df_temp)[names(df_temp)==paste0("a_", a_curr)] <- "a_k"
        names(df_temp)[names(df_temp)==paste0("fit_", a_curr)] <- "fit_k"
        
        m_k <- sum(t_df$a_val < a_curr)
        
        if (d_fun == "TOT"){
          names(df_temp)[names(df_temp)==paste0("w_m", m_k)] <- "w_m_k"
          df_temp <- df_temp %>% mutate(w_nu_i = as.integer(w_m_k >= i))
          names(df_temp)[names(df_temp)=="w_m_k"] <- paste0("w_m", m_k)
        }
        if (d_wrt == "TOT") {
          names(df_temp)[names(df_temp)==paste0("w_m", m_k)] <- "w_m_k"
          df_temp <- df_temp %>% mutate(w_nu_j = as.integer(w_m_k >= j))
          names(df_temp)[names(df_temp)=="w_m_k"] <- paste0("w_m", m_k)
        }
        
        t_df_psi_i <- t_df %>% filter( .$psi_counter == i)
        t_df_psi_j <- t_df %>% filter( .$psi_counter == j)
        
        if ((d_fun == "CT") & (a_curr %in% t_df_psi_i$a_val)){
          if ((d_wrt == "CT") & (a_curr %in% t_df_psi_j$a_val)){
            denom <- denom + with(df_temp, sum(fit_k*(1-fit_k)*tau_k*tau_k))
          }
          if (d_wrt == "TOT") {
            denom <- denom + with(df_temp, sum(fit_k*(1-fit_k)*tau_k*tau_k*w_nu_j))
          }
        }
        if (d_fun == "TOT") {
          if ((d_wrt == "CT") & (a_curr %in% t_df_psi_j$a_val)){
            denom <- denom + with(df_temp, sum(fit_k*(1-fit_k)*tau_k*tau_k*w_nu_i))
          }
          if (d_wrt == "TOT") {
            denom <- denom + with(df_temp, sum(fit_k*(1-fit_k)*tau_k*tau_k*w_nu_i*w_nu_j))
          }
        }
      }
        
      hessian_vec <- c(hessian_vec, denom)
    }
  }
  
  return(matrix(hessian_vec, byrow=T, nrow=(length(psi_hat_CT) + length(psi_hat_TOT) - 1)))
}





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
### ESTIMATION ###
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
newton_raphson_psi <- function(df, t_df,
                               psi_start, psi_choose=0, 
                               tol=0.001, max_iter=20, psi_max = 10){
  psi_hat <- psi_start
  psi_next <- psi_hat[[psi_choose+1]]
  
  psi_old <- psi_next
  steps <- 1L
  fail_flag <- F
  
  while(((sum(abs(psi_next - psi_old)) > tol) | steps == 1) &
        (steps <= max_iter) & (max(abs(psi_next)) < psi_max)){
    
    psi_old <- psi_next
    num <- 0
    denom <- 0
    a_check <- -1
    
    for(i in 1:nrow(t_df)){
      if(t_df$psi_counter[[i]] == psi_choose){
        a_curr <- t_df$a_val[[i]]
        if (a_check < a_curr){
          df_temp <- calculate_tau_k(df, psi_hat=psi_hat,
                                     prm=prm, a_k=a_curr)
          
          names(df_temp)[names(df_temp)==paste0("a_", a_curr)] <- "a_k"
          names(df_temp)[names(df_temp)==paste0("fit_", a_curr)] <- "fit_k"
          
          num <- num + with(df_temp, sum((a_k-fit_k)*tau_k))
        }
        a_check <- a_curr
      }
    }
    
    t_df_psi <- t_df %>% filter( .$psi_counter == psi_choose)
    for(a_curr in t_df_psi$a_val){
      for(i in 1:nrow(t_df_psi)){
        a_dash <- t_df_psi$a_val[[i]]
        psi_dash <- t_df_psi$psi_counter[[i]]
        
        if((a_dash >= a_curr) & (psi_dash == psi_choose)){
          df_temp <- calculate_tau_rsd(df, t_df=t_df,
                                       psi_hat=psi_hat,
                                       a_k=a_dash, psi_k=psi_dash)
          
          names(df_temp)[names(df_temp)==paste0("a_", a_curr)] <- "a_k"
          names(df_temp)[names(df_temp)==paste0("a_", a_dash)] <- "a_dish"
          names(df_temp)[names(df_temp)==paste0("fit_", a_curr)] <- "fit_k"
          
          if(a_dash == a_curr){
            denom <- denom + with(df_temp, sum((a_k-fit_k)*a_k*tau_rsd))
          } else {
            denom <- denom + with(df_temp, sum((a_k-fit_k)*a_dish*tau_rsd))
          }
        }
      }
    }
    
    psi_next <- psi_old + num/denom
    steps <- steps + 1
    psi_hat[[psi_choose+1]] <- psi_next
  }
  
  if (max(abs(psi_next)) >= psi_max){
    fail_flag <- T
  }
  
  
  return(list(psi_hat, steps, fail_flag))
}


newton_raphson_iterative <- function(df, t_df, 
                                     psi_init=rep(0, length(t_psi_vec)),
                                     tol=0.001, max_iter=20){
  
  psi_current <- psi_init
  psi_old <- psi_current
  steps <- 1L
  
  follow_df <- tibble(psi_follow=integer(),
                      psi_before=double(), psi_after=double(),
                      nr_steps=integer())
  fail_flag_iter <- FALSE
  while(((sum(abs(psi_current - psi_old)) > tol) | steps == 1) &
        (steps <= max_iter) & isTRUE(!fail_flag_iter)){
    
    psi_old <- psi_current
    
    for (i in ((length(t_psi_vec):1)-1)){
      nr_out <- newton_raphson_psi(df, t_df=t_df, 
                                   psi_start=psi_current,
                                   psi_choose=i)
      psi_current <- nr_out[[1]]
      
      fail_flag_iter <- nr_out[[3]]
      
      
      follow_df <- follow_df %>% add_row(psi_follow=i,
                                         psi_before=psi_old[[i+1]],
                                         psi_after=psi_current[[i+1]],
                                         nr_steps=nr_out[[2]]) 
    }
    steps <- steps + 1
  }
  
  return(list(psi_current, steps, follow_df, fail_flag_iter))
}






newton_raphson_grad <- function(df, prm,
                                psi_start_list=NA, 
                                tol=0.001, max_iter=20, psi_max = 10){
  t_psi_vec <- prm$t_psi_vec
  xi_vec <- prm$xi_vec
  
  if(is.na(psi_start_list)) {
    psi_start_list <- list(rep(0, length(t_psi_vec)),
                            rep(0, length(xi_vec)))
  }
  
  t_df <- prm$t_df
  
  psi_hat_list <- psi_start_list
  psi_old_list <- psi_hat_list
  steps <- 1L
  fail_flag <- F
  
  while(((sum(abs(unlist(psi_hat_list) - unlist(psi_old_list))) > tol) | steps == 1) &
        (steps <= max_iter) & (max(abs(unlist(psi_hat_list))) < psi_max)){
    
    psi_old_list <- psi_hat_list
    
    S_vec <- df %>% calculate_score(prm=prm, psi_hat_list=psi_hat_list)
    J <- df %>% calculate_jacobian(prm=prm, psi_hat_list=psi_hat_list)

    
    psi_hat_vec <- psi_hat_list[[1]]
    if (length(psi_hat_list[[2]]) > 1){
      psi_hat_vec <- c(psi_hat_vec, psi_hat_list[[2]][2:length(psi_hat_list[[2]])])
    }
    
    psi_hat_vec <- solve(J, S_vec) + psi_hat_vec
    psi_hat_list <- list(psi_hat_vec[1:length(t_psi_vec)],
                          c(0, if(length(psi_hat_vec) > length(t_psi_vec)){
                                  psi_hat_vec[(length(t_psi_vec)+1):length(psi_hat_vec)]}
                            else{c()}))
    psi_hat_vec
    
    steps <- steps + 1
  }
  
  #if (max(abs(psi_hat_list)) >= psi_max){
  #  fail_flag <- T
  #}
  return(list(psi_hat_list, steps, fail_flag))
}




calculate_variance_OLD <- function(df, psi_hat, t_a_vec, t_psi_vec, t_df){
  
  D_vec <- c()
  row_count <- 0
  for (i in 0:(length(t_a_vec)-1)) {
    for (j in -1:(i-1)){
      for (k in 0:(length(t_a_vec)-1)) {
        df_temp <- df %>% filter( .$ti > t_a_vec[[k+1]])
        nk <- nrow(df_temp)
        if (k == i) {
          if (j == -1){
            D_vec <- c(D_vec, rep(1, nk))
          } else {
            D_vec <-  c(D_vec, df_temp[[paste0('a_', j)]])
          }
        } else {
          D_vec <- c(D_vec, rep(0, nk))
        }
      }
      row_count <- row_count + 1
    }
  }
  D <- matrix(D_vec, nrow=row_count, byrow=T)
  D <- t(D)
  
  
  fit <- c()
  for (i in 0:(length(t_a_vec)-1)){
    df_temp <- df %>% filter( .$ti > t_a_vec[[i+1]])
    fit <- c(fit, df_temp[[paste0('fit_', i)]])
  }
  
  Jbb_vec <- c()
  for(i in 1:dim(D)[2]){
    for(j in 1:dim(D)[2]){
      Jbb_vec <- c(Jbb_vec, sum(D[,i]*D[,j]*fit*(1-fit)))
    }
  }
  Jbb <- matrix(Jbb_vec, nrow=sqrt(length(Jbb_vec)), byrow=T)
  
  tau_vec <- c()
  for (i in 0:(length(t_psi_vec)-1)){
    for (k in 0:(length(t_a_vec)-1)){
      df_temp <- df %>% filter( .$ti >= t_a_vec[[k+1]])
      nk <- nrow(df_temp)
      if (ifelse(k == (length(t_a_vec)-1), T, 
                 (t_psi_vec[[i+1]] < t_a_vec[[k+1+1]])) & 
          ifelse(i == (length(t_psi_vec)-1), T, 
                 t_psi_vec[[i+1+1]] > t_a_vec[[k+1]])) {
        df_temp <- calculate_tau_k(df_temp, psi_hat=psi_hat,
                                   prm=prm, a_k=k)
        tau_vec <- c(tau_vec, df_temp$tau_k)
      } else {
        tau_vec <- c(tau_vec, rep(0,nk))
      }
    }
  }
  tau <- matrix(tau_vec, ncol=length(t_psi_vec), byrow=F)
  
  
  
  Jtb <- t(sapply(1:dim(tau)[2],
                  function(i){sapply(1:dim(D)[2], function(j){
                    sum(D[,j]*tau[,i]*fit*(1-fit))})
                  }, simplify=T))
  
  
  Jtt <- sapply(1:dim(tau)[2],
                function(i){ sapply(1:dim(tau)[2],
                                    function(j) { sum(tau[,i]*tau[,j]*fit*(1-fit)) 
                                    }
                )}
  )
  
  
  JTT <- Jtt - Jtb%*%solve(Jbb)%*%t(Jtb)
  
  
  
  Stp_vec <- c()
  for (i in 0:(length(t_psi_vec)-1)){
    for (j in 0:(length(t_psi_vec)-1)){
      t_df_i <- t_df %>% filter( .$psi_counter == i)
      ds_dpsi <- 0
      for (a_i in unique(t_df_i$a_val)){
        t_df_a <- t_df %>% filter( .$a_val == a_i)
        psi_min <- min(t_df_a$psi_counter)
        
        if (j >= psi_min){
          df_rsd <- calculate_dtau_dpsi(df=df, t_df=t_df, psi_hat=psi_hat,
                                        a_k=a_i,
                                        psi_k=j)
          
          names(df_rsd)[names(df_rsd)==paste0("a_", a_i)] <- "a_k"
          names(df_rsd)[names(df_rsd)==paste0("fit_", a_i)] <- "fit_k"
          ds_dpsi <- ds_dpsi - sum(with(df_rsd, (a_k - fit_k)*dtau_dpsi))
        }
      }
      Stp_vec <- c(Stp_vec, ds_dpsi)
    }
  }
  Stp <- matrix(Stp_vec, nrow=length(t_psi_vec), byrow=T)
  Stp_inv <- solve(Stp)
  
  #Stp_1_inv%*%JTT%*%t(Stp_1_inv)
  #return(solve(Stp_1)%*%JTT%*%t(solve(Stp_1)))
  return(solve(Stp)%*%JTT%*%t(solve(Stp)))
}



calculate_variance <- function(df, psi_hat_list, prm){
  
  t_a_vec <- prm$t_a_vec
  t_psi_vec <- prm$t_psi_vec
  xi_vec <- prm$xi_vec
  
  t_df <- prm$t_df
  
  psi_hat_CT <- psi_hat_list[[1]]
  psi_hat_TOT <- psi_hat_list[[2]]
  
  D_vec <- c()
  row_count <- 0
  for (i in 0:(length(t_a_vec)-1)) {
    for (j in -1:(i-1)){
      for (k in 0:(length(t_a_vec)-1)) {
        df_temp <- df %>% filter( .$ti > t_a_vec[[k+1]])
        nk <- nrow(df_temp)
        if (k == i) {
          if (j == -1){
            D_vec <- c(D_vec, rep(1, nk))
          } else {
            D_vec <-  c(D_vec, df_temp[[paste0('a_', j)]])
          }
        } else {
          D_vec <- c(D_vec, rep(0, nk))
        }
      }
      row_count <- row_count + 1
    }
  }
  D <- matrix(D_vec, nrow=row_count, byrow=T)
  D <- t(D)
  

  D_psi_vec <- c()
  for (j in 1:(length(psi_hat_CT) + length(psi_hat_TOT) - 1)){
    for (a_curr in 0:(length(t_a_vec)-1)){
      df_temp <- df %>% filter( .$ti > t_a_vec[[a_curr+1]])
      
      if (j <= length(psi_hat_CT)) {
        psi_choose <- j - 1
        t_df_psi <- t_df %>% filter( .$psi_counter == psi_choose)
        if (a_curr %in% t_df_psi$a_val) {
          D_psi_vec <- c(D_psi_vec, rep(1, dim(df_temp)[[1]]))
        } else {
          D_psi_vec <- c(D_psi_vec, rep(0, dim(df_temp)[[1]]))
        }
      } else {
        nu <- j - length(psi_hat_CT)
        m_k <- sum(t_df$a_val < a_curr)
        
        names(df_temp)[names(df_temp)==paste0("w_m", m_k)] <- "w_m_k"
        df_temp <- df_temp %>% mutate(w_nu = as.integer(w_m_k >= nu))
        
        D_psi_vec <- c(D_psi_vec, df_temp$w_nu)
      }
    }
  }
  D_psi <- matrix(D_psi_vec, byrow=F, ncol=length(psi_hat_CT) + length(psi_hat_TOT) - 1)
  
  
  
  fit <- c()
  tau_vec <- c()
  for (a_curr in 0:(length(t_a_vec)-1)){
    df_temp <- df %>% filter( .$ti > t_a_vec[[a_curr+1]]) %>% 
      calculate_tau_k(psi_hat_list=psi_hat_list, prm=prm, a_k = a_curr)
    fit <- c(fit, df_temp[[paste0('fit_', a_curr)]])
    tau_vec <- c(tau_vec, df_temp$tau_k)
  }
  
  
  Jbb_vec <- c()
  for(i in 1:dim(D)[2]){
    for(j in 1:dim(D)[2]){
      Jbb_vec <- c(Jbb_vec, sum(D[,i]*D[,j]*fit*(1-fit)))
    }
  }
  Jbb <- matrix(Jbb_vec, nrow=sqrt(length(Jbb_vec)), byrow=T)
  
  
  #Jtt <- df %>% calculate_hessian(prm=prm, psi_hat_list=psi_hat_list)
  
  Jtt <- t(sapply(1:dim(D_psi)[[2]], function(i){
           sapply(1:dim(D_psi)[[2]], function(j){
              sum(D_psi[,i]*D_psi[,j]*tau_vec*tau_vec*fit*(1-fit))
            })
          }))
  
  Jtb <- t(sapply(1:dim(D_psi)[[2]], function(i){
            sapply(1:dim(D)[[2]], function(j){
              sum(D_psi[,i]*D[,j]*tau_vec*fit*(1-fit))
            })
          }))
  
  
  JTT <- Jtt - Jtb%*%solve(Jbb)%*%t(Jtb)
  
  Jacobian <- df %>% calculate_jacobian(prm=prm, psi_hat_list=psi_hat_list)
  J_inv <- solve(Jacobian)
  
  return(J_inv %*% JTT %*% t(J_inv))
}




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
### SPOOKY TESTING REALM ###
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

test_plot_1 <- function(x_coords = seq(-2, 2, 0.2), y_coords = seq(-2, 2, 0.2)){
  
  n <- 1000
  t0_min <- 10
  t0_max <- 100
  
  psi_star_CT <- c(log(2))
  t_psi_vec <- c(0)
  
  psi_star_TOT <- c(0, log(1.5))
  xi_vec <- c(0, 30)
  
  psi_star_list <- list(psi_star_CT, psi_star_TOT)
  t_a_vec   <- c(0, 30)
  
  TOT_times_vec <- sort(unique(as.vector(sapply(xi_vec, function(q){t_a_vec + q}))))
  t_df <- get_time_df(t_a_vec, t_psi_vec, TOT_times_vec)
  
  df <- create_sample(t_df, psi_star_list=psi_star_list)
  df <- df %>% fit_treatment_models(t_a_vec=t_a_vec, 
                                    t_psi_vec=t_psi_vec,
                                    xi_vec=xi_vec)
  
  score_CT <- c()
  score_TOT <- c()
  for (psi_CT in x_coords){
    for (psi_TOT in y_coords){
      S <- df %>% calculate_score(t_df=t_df, psi_hat_list=list(c(psi_CT), c(0,psi_TOT)))
      score_CT <- c(score_CT, S[[1]])
      score_TOT <- c(score_TOT, S[[2]])
    }
  }
  
  score_CT <- score_CT / max(score_CT)
  score_TOT <- score_TOT / max(score_TOT)
  
  old_pars <- par(mfrow=c(2,2))
  image(x_coords, y_coords, matrix(score_CT, nrow=length(x_coords), byrow=T))
  contour(x_coords, y_coords, matrix(score_CT, nrow=length(x_coords), byrow=T),
          add = TRUE, levels = seq(-1,1,0.2))
  points(psi_star_list[[1]][[1]], psi_star_list[[2]][[2]])
  
  image(x_coords, y_coords, matrix(score_TOT, nrow=length(x_coords), byrow=T))
  contour(x_coords, y_coords, matrix(score_TOT, nrow=length(x_coords), byrow=T),
           add = TRUE, levels = seq(-1,1,0.2))
  points(psi_star_list[[1]][[1]], psi_star_list[[2]][[2]])
  
  image(x_coords, y_coords, matrix(score_CT, nrow=length(x_coords), byrow=T))
  contour(x_coords, y_coords, matrix(score_CT, nrow=length(x_coords), byrow=T),
          add = TRUE, levels = seq(-1,1,0.2))
  contour(x_coords, y_coords, matrix(score_TOT, nrow=length(x_coords), byrow=T),
          add = TRUE, levels = seq(-1,1,0.2))
  points(psi_star_list[[1]][[1]], psi_star_list[[2]][[2]])
  
  image(x_coords, y_coords, matrix(score_TOT, nrow=length(x_coords), byrow=T))
  contour(x_coords, y_coords, matrix(score_CT, nrow=length(x_coords), byrow=T),
          add = TRUE, levels = seq(-1,1,0.2))
  contour(x_coords, y_coords, matrix(score_TOT, nrow=length(x_coords), byrow=T),
          add = TRUE, levels = seq(-1,1,0.2))
  points(psi_star_list[[1]][[1]], psi_star_list[[2]][[2]])
  par(old_pars)
}

test_zone <- function(){
  
  psi_star_CT <- c(log(2))
  t_psi_vec <- c(0)
  psi_star_TOT <- c(0, log(1.5))
  xi_vec <- c(0, 30)
  t_a_vec  <- c(0, 30)
  
  n <- 10
  t0_min <- 10
  t0_max <- 100
  sims <- 1
  
  prm <- list(psi_star_CT=psi_star_CT, 
              t_psi_vec=t_psi_vec,
              psi_star_TOT=psi_star_TOT, 
              xi_vec=xi_vec,
              psi_star_list=list(psi_star_CT, psi_star_TOT),
              t_a_vec=t_a_vec,
              n=n,
              t0_min=t0_min,
              t0_max=t0_max,
              sims=sims)
  
  prm$t_df <- get_time_df(prm=prm)
  
  df <- create_sample(prm=prm)
  df <- df %>% fit_treatment_models(prm=prm)
  
  nri_out <- df %>% newton_raphson_grad(prm=prm)
  
  psi_hat_list <- nri_out[[1]]
  
  df %>% calculate_variance(psi_hat_list=psi_hat_list, prm=prm)
  
}

nr_run <- function(psi_star_CT, t_psi_vec, 
                   psi_star_TOT, xi_vec,
                   t_a_vec,
                   n = 1000,
                   t0_min = 10,
                   t0_max = 100,
                   sims = 10){
  

  prm <- list(psi_star_CT=psi_star_CT, 
              t_psi_vec=t_psi_vec,
              psi_star_TOT=psi_star_TOT, 
              xi_vec=xi_vec,
              psi_star_list=list(psi_star_CT, psi_star_TOT),
              t_a_vec=t_a_vec,
              n=n,
              t0_min=t0_min,
              t0_max=t0_max,
              sims=sims)
  
  prm$t_df <- get_time_df(prm=prm)
  
  cl <- makePSOCKcluster(28)
  registerDoParallel(cl)
  registerDoParallel()

  nr_out <- foreach(i=c(1:sims), .combine=rbind,
                    .export=c("create_sample", "fit_treatment_models",
                              "calculate_tau_k", "calculate_tau_rsd",
                              "calculate_tau_rsd_m", "calculate_jacobian", 
                              "calculate_hessian",
                              "calculate_score", "newton_raphson_grad",
                              "calculate_variance"),
                    .packages="tidyverse") %dopar%
    
    {
      df <- create_sample(prm=prm)
      df <- df %>% fit_treatment_models(prm=prm)
      
      nri_out <- df %>% newton_raphson_grad(prm=prm)
      psi_hat_list <- nri_out[[1]]
      
      var_hat <- df %>% calculate_variance(psi_hat_list=psi_hat_list, prm=prm)
      
      c(unlist(psi_hat_list), diag(var_hat))
    }
  
  stopCluster(cl)
  
  psi_star_comb <- c(psi_star_CT, psi_star_TOT)
  
  offset_CT <- length(t_psi_vec)
  offset_psi <- length(psi_star_comb)
  
  for (j in 1:length(psi_star_comb)){
    if (j <= offset_CT) {
      psi_lab <- paste0("psi CT (", j-1, ")")
      var_index <- j + offset_psi
    }
    else if (j >= offset_CT + 2) {
      psi_lab <- paste0("psi TOT (", j-1 - offset_CT, ")")
      var_index <- j + offset_psi - 1
    } else {
      next
    }
    
    cat(psi_lab, "\n",
        paste(
          format(c(psi_star_comb[[j]], 
                   mean(nr_out[,j], na.rm=T))
                 , digits = 5),
          collapse = "\t"),
        "\n",
        paste(
          format(c(var(nr_out[,j]),
                   mean(nr_out[,var_index], na.rm=T))
                 , digits = 5),
          collapse = "\t"),
        "\n")
  }
  
  #return(NULL)
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
### MAIN # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

nr_run(psi_star_CT = c(log(2)), 
       t_psi_vec = c(0),
       psi_star_TOT = c(0, log(1.5)),
       xi_vec = c(0, 30),
       t_a_vec  = c(0, 30),
       n = 1000,
       sims=1000
)


nr_run(psi_star_CT = c(log(2)), 
       t_psi_vec = c(0),
       psi_star_TOT = c(0),
       xi_vec = c(0),
       t_a_vec  = c(0, 30),
       n = 1000,
       sims=5000
)


nr_run(psi_star_CT = c(log(2), log(1.5)), 
       t_psi_vec = c(0, 30),
       psi_star_TOT = c(0),
       xi_vec = c(0),
       t_a_vec  = c(0, 30),
       n = 1000,
       sims=5000
)



nr_run(psi_star_CT = c(log(2)), 
       t_psi_vec = c(0),
       psi_star_TOT = c(0, log(1.5), log(0.5)),
       xi_vec = c(0, 30, 60),
       t_a_vec  = c(0, 30, 60),
       t0_max = 120,
       n = 1000,
       sims=1000
)



# prm <- list(psi_star_CT = c(log(2)), 
#        t_psi_vec = c(0),
#        psi_star_TOT = c(0, log(1.5)),
#        xi_vec = c(0, 30),
#        t_a_vec  = c(0, 30),
#        n = 1000,
#        t0_min = 10,
#        t0_max = 100,
#        sims=100
# )
# 
# prm$psi_star_list <- list(prm$psi_star_CT, prm$psi_star_TOT)





function_test(psi_star_CT = c(log(2)), 
              t_psi_vec = c(0),
              psi_star_TOT = c(0, log(1.5)),
              xi_vec = c(0, 30),
              t_a_vec  = c(0, 30))
       

 # PSI's for no TOT dependency
#psi_star <- c(log(2))
#psi_star <- c(log(1.8), log(1.5), log(2))
#psi_star <- c(log(1.8), log(1.5), log(2))

#PSI's when CT and TOT dependency (should also work in general)
psi_star_CT <- c(log(2))
t_psi_vec <- c(0)
#psi_star_CT <- c(log(2), log(1.5))
#t_psi_vec <- c(0, 40)

#PSI's for TOT dependency. 
# NOTE: 1st element must be zero.
# AND the treatment effect is the cumulative sum

#psi_star_TOT <- c(0)
#xi_vec <- c(0)
psi_star_TOT <- c(0, log(1.5))
xi_vec <- c(0, 30)
#psi_star_TOT <- c(0, log(1.2), log(0.3))
#xi_vec <- c(0, 5, 20)



t_a_vec   <- c(0, 30)

psi_star_list <- list(psi_star_CT, psi_star_TOT)
t_df <- get_time_df(t_a_vec, t_psi_vec, TOT_times_vec)

TOT_times_vec <- sort(unique(as.vector(sapply(xi_vec, function(q){t_a_vec + q}))))





#df <- create_sample(t_df, psi_star_list=psi_star_list)
#df <- df %>% fit_treatment_models(t_a_vec=t_a_vec, 
#                                  t_psi_vec=t_psi_vec,
#                                  xi_vec=xi_vec)

#df %>% newton_raphson_grad(t_df=t_df)


nr_out <- tibble(psi_hat_CT_0=numeric(),
                 psi_hat_TOT_0=numeric(),
                 iters=integer()
                 )

nr_out <- foreach(i=c(1:sims), .combine=rbind) %dopar%
  
  {
    
    df <- create_sample(t_df, psi_star_list=psi_star_list)
    df <- df %>% fit_treatment_models(t_a_vec=t_a_vec, 
                                      t_psi_vec=t_psi_vec,
                                      xi_vec=xi_vec)
    
    #nri_out <- newton_raphson_iterative(df, t_df=t_df)
    nri_out <- newton_raphson_grad(df, t_df=t_df)
    
    unlist(nri_out[[1]])
    
  }




nr_out <- foreach(i=c(1:sims), 
                  .combine=rbind,
                  .packages='tidyverse') %dopar% 
  {
    
    df <- create_sample(t_df, psi_star_list=psi_star_list)
    df <- df %>% fit_treatment_models(t_a_vec=t_a_vec, 
                                      t_psi_vec=t_psi_vec,
                                      xi_vec=xi_vec)
    
    #nri_out <- newton_raphson_iterative(df, t_df=t_df)
    nri_out <- newton_raphson_grad(df, t_df=t_df)
    
    tibble(psi_hat_CT_0 = nri_out[[1]][[1]],
           psi_hat_TOT_0 = nri_out[[1]][[2]][[2]],
           inters = nri_out[[2]]
    )
  }

cat("psi(0) mean","\n","true mean --- s.mean of est. means","\n",
    paste(format(
      with(nr_out, c(psi_star_CT[[1]], mean(psi_hat_CT_0, na.rm=T)))
      , digits = 5), collapse = "   "),

    "\n\n","psi(1) mean", "\n", "true mean --- s.mean of est. means", "\n",  
    paste(format(
      with(nr_out, c(psi_star_TOT[[2]], mean(psi_hat_TOT_0, na.rm=T)))
      , digits = 5), collapse = "   "),

    sep="")













nr_out <- tibble(psi_hat_0=numeric(), var_psi_0=numeric()
                 , psi_hat_1=numeric(), var_psi_1=numeric(),
                 fail_flag=logical()
                 #, var_psi_0_alt=numeric(), var_psi_1_alt=numeric()
                 #,psi_hat_2=numeric(), var_psi_2=numeric()
)
















nr_out <- foreach(i=c(1:sims), 
                  .combine=rbind,
                  .packages='tidyverse') %dopar% 
  {
    output_check <- F
    while(!isTRUE(output_check)){
      nr_sim <- tryCatch(
        expr = {
          df <- create_sample(t_df, psi_star_CT=psi_star_CT)
          df <- df %>% fit_treatment_models(t_a_vec=t_a_vec, 
                                            t_psi_vec=t_psi_vec)
          
          #nri_out <- newton_raphson_iterative(df, t_df=t_df)
          nri_out <- newton_raphson_grad(df, t_df=t_df)
          psi_hat <- nri_out[[1]]
          fail_check <- nri_out[[3]]
          #fail_check <- nri_out[[4]]
          
          VCV <- calculate_variance(df=df, psi_hat=psi_hat,
                                    t_a_vec=t_a_vec, t_psi_vec=t_psi_vec,
                                    t_df=t_df)
          if(max(abs(VCV)) >= max(abs(psi_hat))){
            fail_check <- T
          }
          
          tibble(psi_hat_0 = psi_hat[[1]], var_psi_0 = VCV[[1,1]]
                 , psi_hat_1 = psi_hat[[2]], var_psi_1 = VCV[[2,2]]
                 , fail_flag = fail_check
                 #, var_psi_0_alt = (solve(VCV)[[1,1]])^-1
                 #, var_psi_1_alt = (solve(VCV)[[2,2]])^-1
                 #, psi_hat_2 = psi_hat[[3]], var_psi_2 = VCV[[3,3]]
          )
        },
        error = function(e){conditionMessage(e)},
        warning = function(w){conditionMessage(w)}
      )
      
      if(typeof(nr_sim) == "list"){output_check <- T}
    }
    nr_sim
  }

#print("psi(1) mean", "test")
nr_out_pass <- nr_out  %>% filter(.$fail_flag == F)
cat("psi(0) mean","\n","true mean --- s.mean of est. means","\n",
    paste(format(
      with(nr_out_pass, c(psi_star[[1]], mean(psi_hat_0, na.rm=T)))
      , digits = 5), collapse = "   "),
    "\n\n","psi(0) var", "\n", "var of s.mean --- s.mean of est. var", "\n",
    paste(format(
      with(nr_out_pass, c(var(psi_hat_0, na.rm=T), mean(var_psi_0, na.rm=T)))
      , digits = 5), collapse = "   "),
    "\n\n","psi(1) mean", "\n", "true mean --- s.mean of est. means", "\n",  
    paste(format(
      with(nr_out_pass, c(psi_star[[2]], mean(psi_hat_1, na.rm=T)))
      , digits = 5), collapse = "   "),
    "\n\n","psi(1) var", "\n", "var of s.mean --- s.mean of est. var", "\n",
    paste(format(
      with(nr_out_pass, c(var(psi_hat_1, na.rm=T), mean(var_psi_1, na.rm=T)))
      , digits = 5), collapse = "   "),
    paste("\n\n", "FAILED ", with(nr_out, sum(fail_flag)), " OF ", sims,
          sep=""),
    sep="")

#with(nr_out, c(psi_star[[3]], mean(psi_hat_2)))
#with(nr_out, c(var(psi_hat_2), mean(var_psi_2)))

