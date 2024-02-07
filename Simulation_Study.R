rm(list=ls())
gc()

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


source("G-estimation_FINAL.R")


create_sample <- function(prm){
  
  a_vars <- length(prm$t_a_vec)
  k_parts <- 2^(a_vars + 1)
  
  n_s <- ceiling(prm$n_trgt/k_parts)
  prm$n <- n_s*k_parts
  
  df <- tibble(id=1:prm$n, 
               x = c(rep(0, n_s*(2^a_vars)), rep(1, n_s*(2^a_vars))))
  
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
      temp = pmin((t_next - t_curr)*exp(g_psi*a_curr), t0_rsd),
      ti = ti + temp*exp(-g_psi*a_curr),
      t0_rsd = t0_rsd - temp
    )
    
    names(df)[names(df)=="a_curr"] <- paste0("a_", i-1) 
  }
  
  if(prm$censor){
    df <- df %>% mutate(ti = pmin(ti, C_i))
  }
  
  df <- df %>% dplyr::select(-c(t0_rsd, g_psi, temp))
  
  return(df)
}

create_sample_ANZDATA <- function(prm){

  prm$n <- prm$n_trgt
  
  df <- tibble(id=1:prm$n, 
               x = rep(0, prm$n))
  
  if(prm$censor){
    df <- df %>% add_column(C_i = prm$censor_date,
                            a_0 = rbinom(prm$n, size=1, prob=0.25))
  }
  
  for (k in 1:(length(prm$t_a_vec)-1)){
      a_prev <- as.name(paste0("a_", k-1))
      a_k <- as.name(paste0("a_", k))
      df <- df %>% add_column(switches = rbinom(prm$n, size=1, prob=0.05)) %>%
                mutate("{a_k}" := switches*(1-{{a_prev}}) + (1-switches)*{{a_prev}}) %>%
                select(-c(switches))
  }

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
      temp = pmin((t_next - t_curr)*exp(g_psi*a_curr), t0_rsd),
      ti = ti + temp*exp(-g_psi*a_curr),
      t0_rsd = t0_rsd - temp
    )
    
    names(df)[names(df)=="a_curr"] <- paste0("a_", i-1) 
  }
  
  if(prm$censor){
    df <- df %>% mutate(ti = pmin(ti, C_i))
  }
  
  df <- df %>% dplyr::select(-c(t0_rsd, g_psi, temp))
  
  return(df)
}



nr_run <- function(prm, collapse=T, piecewise=F){
  
  cl <- makePSOCKcluster(28)
  registerDoParallel(cl)
  
  nr_out_list <- foreach(seed_curr=1:prm$sims, 
                         .packages=c("tidyverse", "DescTools"),
                          .errorhandling='remove',
                    .export=ls(envir=globalenv())[ls(envir=globalenv()) != "prm"]
  ) %dopar%
    {
    
    nr_out_single <- rep(NA, 2*length(prm$psi_star_vec))
    
    set.seed(seed_curr)
    df <- create_sample(prm=prm)
    
    fit_trt_out <- df %>% fit_treatment_models(prm=prm)
    df <- fit_trt_out[[1]]
    prm$trt_models <- fit_trt_out[[2]]
    
    
    nri_out <- df %>% newton_raphson_grad(prm=prm, print_results=T, 
                                          psi_max_vec=prm$psi_max_vec, max_iter=10)
    (psi_hat_vec <- nri_out[[1]])
    (var_hat <- df %>% calculate_variance(psi_hat_vec=psi_hat_vec, prm=prm))
    
    

    list(psi_hat_vec=psi_hat_vec, var_hat=var_hat, 
         n_steps = nri_out[[2]], fail_flag=nri_out[[3]])
    
    }
  
   stopCluster(cl)

   if(collapse == T){
   return(lapply(1:length(nr_out_list), function(i){
      c(nr_out_list[[i]]$psi_hat_vec, diag(nr_out_list[[i]]$var_hat))}))
    } else {
     return(nr_out_list)
   }

}

calculate_and_display <- function(prm, param_list, calc_var_var=F){
  cat(prm$sim_label, "\n", "n = ", prm$n, "\n\n")
  
  for(item in param_list){
    prm$psi_1_star <- item$psi_1_star
    prm$psi_x_star <- item$psi_x_star
    prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
    prm$censor <- item$censoring
    
    nr_out <- nr_run(prm=prm) %>% Reduce(rbind, .)
  
    cat("censoring=", prm$censor, "\n", sep="")
    for (j in 1:length(prm$psi_star_vec)){
          psi_star <- prm$psi_star_vec[j]
          m_psi <- mean(nr_out[,j] , na.rm=T)
          v_psi <- var(nr_out[,j], na.rm=T)
          m_sigmasq <- mean(nr_out[,j+length(prm$psi_star_vec)], na.rm=T)
          rem <- ifelse(psi_star == 0, NA,  (m_psi - psi_star)/psi_star)
          rev <- (m_sigmasq - v_psi)/v_psi 
      cat(prm$psi_lab[[j]],
          psi_star,
          m_psi,
          v_psi,
          m_sigmasq,
          rem,
          rev,
          sep=",   ")
      if(calc_var_var){cat("  varvar=", var(nr_out[,j+length(prm$psi_star_vec)], na.rm=T))}
      cat("\n")
      }
    cat("\n")
  }
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



#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           (1) -> (1) block                            #   ##
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##

prm <- list()
prm$sim_label <- "(const_1) -> (const_1)"
prm$t_a_vec <- c(0, 35)
prm$expmean <- 50
prm$n_trgt <- 400
prm$beta_1_track <- c(1, 1)
prm$beta_x_track <- c(0, 0)
prm$psi_lab <- c("psi_1")
prm$psi_x_star <- c()
prm$sims <- 3000
prm$censor_date <- 80
prm$censor_max <- prm$censor_date
prm$trt_mod_list <- list(
  c("1"),
  c("1")
)


calculate_and_display(prm=prm,
  list(
    list(psi_1_star = c(log(2)), psi_x_star = NULL, censoring=F),
    list(psi_1_star = c(log(2)), psi_x_star = NULL, censoring=T),
    list(psi_1_star = c(log(1)), psi_x_star = NULL, censoring=F),
    list(psi_1_star = c(log(1)), psi_x_star = NULL, censoring=T),
    list(psi_1_star = c(log(1/2)), psi_x_star = NULL, censoring=F),
    list(psi_1_star = c(log(1/2)), psi_x_star = NULL, censoring=T)
  ))



#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           (1, x1) -> (1, x1)  test block
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##

prm <- list()
prm$sim_label <- "(const_1, x_1) -> (const_1, x_1)"
prm$t_a_vec <- c(0, 35)
prm$expmean <- 50
prm$n_trgt <- 400
prm$beta_1_track <- c(1, 1)
prm$beta_x_track <- c(1, 1)
prm$psi_lab <- c("psi_1", "psi_x1")
prm$sims <- 3000
prm$censor_date <- 80
prm$censor_max <- prm$censor_date
prm$trt_mod_list <- list(
  c("1", "x"),
  c("1", "x")
)
prm$psi_max_vec <- c(3,3)

calculate_and_display(prm=prm,
  list(
    list(psi_1_star = c(log(2)), psi_x_star = c(log(1.5)), censoring=F),
    list(psi_1_star = c(log(2)), psi_x_star = c(log(1.5)), censoring=T),
    list(psi_1_star = c(log(2)), psi_x_star = c(log(1/1.5)), censoring=F),
    list(psi_1_star = c(log(2)), psi_x_star = c(log(1/1.5)), censoring=T),
    list(psi_1_star = c(log(1/2)), psi_x_star = c(log(1.5)), censoring=F),
    list(psi_1_star = c(log(1/2)), psi_x_star = c(log(1.5)), censoring=T),
    list(psi_1_star = c(log(1/2)), psi_x_star = c(log(1/1.5)), censoring=F),
    list(psi_1_star = c(log(1/2)), psi_x_star = c(log(1/1.5)), censoring=T)
  ))



#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           (1) -> (2) block                            #   ##
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##

prm <- list()
prm$sim_label <- "(const_1) -> (const_2)"
prm$t_a_vec <- c(0, 35)
prm$expmean <- 50
prm$n_trgt <- 400
prm$beta_1_track <- c(1, 2)
prm$beta_x_track <- c(0, 0)
prm$psi_lab <- c("psi_1", "psi_2")
prm$sims <- 3000
prm$censor_date <- 80
prm$censor_max <- prm$censor_date
prm$trt_mod_list <- list(
  c("1"),
  c("1")
)
prm$psi_max_vec <- c(2,2)

prm$psi_1_star <- rep(log(2), 2)
prm$psi_x_star <- NULL
prm$psi_star_vec <- prm$psi_1_star


calculate_and_display(prm=prm,
  list(
    list(psi_1_star = c(log(2), log(2)), psi_x_star = NULL, censoring=F),
    list(psi_1_star = c(log(2), log(2)), psi_x_star = NULL, censoring=T),
    list(psi_1_star = c(log(2), log(1/2)), psi_x_star = NULL, censoring=F),
    list(psi_1_star = c(log(2), log(1/2)), psi_x_star = NULL, censoring=T),
    list(psi_1_star = c(log(1/2), log(2)), psi_x_star = NULL, censoring=F),
    list(psi_1_star = c(log(1/2), log(2)), psi_x_star = NULL, censoring=T),
    list(psi_1_star = c(log(1/2), log(1/2)), psi_x_star = NULL, censoring=F),
    list(psi_1_star = c(log(1/2), log(1/2)), psi_x_star = NULL, censoring=T)
  ), calc_var_var=T)

for(n in c(400, 600, 800, 1000, 1200)){
 prm$n_trgt <- n
 calculate_and_display(prm=prm,
  list(
    list(psi_1_star = c(log(2), log(2)), psi_x_star = NULL, censoring=F)
  ))
}







#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           (1) -> ... -> (k) block                            #   ##
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##

prm <- list()
prm$sim_label <- "(1) -> ... -> (k)"
prm$expmean <- 50
prm$n_trgt <- 4000
prm$sims <- 3000
prm$censor_date <- 100
prm$censor_max <- prm$censor_date

for (k in c(3,5,7,9)){
  prm$t_a_vec <- (80 %/% (k))*(0:(k-1))
  prm$beta_1_track <- 1:k
  prm$beta_x_track <- rep(0,k)
  prm$psi_lab <- paste0("psi_", 1:k)
  prm$trt_mod_list <- lapply(1:k, function(x){c("1")})
  prm$psi_max_vec <- rep(2, k)
  prm$censor <- T
  
  calculate_and_display(prm=prm,
      list(
        list(psi_1_star = rep(log(1/2), k),
           psi_x_star = NULL,
           censoring=T)
      ))
  calculate_and_display(prm=prm,
      list(
        list(psi_1_star = rep(log(2), k),
           psi_x_star = NULL,
           censoring=T)        
      ))
}





#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           ANZDATA                           #   ##
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##

prm <- list()
prm$sim_label <- "ANZDATA"
prm$expmean <- 17
prm$n_trgt <- 12000
prm$sims <- 100
prm$censor <- T
prm$t_a_vec <- 0:32
prm$beta_1_track <- 1:33
prm$beta_x_track <- rep(0, 33)
prm$trt_mod_list <- lapply(1:33, function(x){c("1")})
prm$psi_max_vec <- rep(1.5*log(4), 33)
prm$psi_1_star <- seq(from=log(4), to=log(1.5), length.out=33)
prm$psi_x_star <- NULL
prm$psi_star_vec <- prm$psi_1_star
prm$censor_date <- 33
prm$censor_max <- prm$censor_date

df <- create_sample_ANZDATA(prm)


for (k in c(3,5,7,9)){
  prm$t_a_vec <- (80 %/% (k))*(0:(k-1))
  prm$beta_1_track <- 1:k
  prm$beta_x_track <- rep(0,k)
  prm$psi_lab <- paste0("psi_", 1:k)
  prm$trt_mod_list <- lapply(1:k, function(x){c("1")})
  prm$psi_max_vec <- rep(2, k)
  prm$censor <- T
  
  calculate_and_display(prm=prm,
      list(
        list(psi_1_star = rep(log(1/2), k),
           psi_x_star = NULL,
           censoring=T)
      ))
  calculate_and_display(prm=prm,
      list(
        list(psi_1_star = rep(log(2), k),
           psi_x_star = NULL,
           censoring=T)        
      ))
}