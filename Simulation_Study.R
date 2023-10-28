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
library(gt)


source("G-estimation_FINAL.R")


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
  
  df <- tibble(id=1:prm$n_trgt, 
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


nr_run <- function(prm){
  
  cl <- makePSOCKcluster(28)
  registerDoParallel(cl)
  
  nr_out_list <- foreach(seed_curr=1:prm$sims, 
                         .packages=c("tidyverse", "DescTools"), #.combine='comb',
                          .errorhandling='remove',
                    .export=c(
                      "create_sample", "fit_treatment_models",
                      "calculate_tau_k", "calculate_tau_rsd_m",
                      "calculate_score", "calculate_jacobian",
                      "calculate_C_k", "calculate_C_m",
                      "newton_raphson_grad",
                      "calculate_variance"
                    )) %dopar%
    {
  
  # nr_out_list <- lapply(1:prm$sims, function(simnum){
    
    # if(simnum==1){
    #   progressbar <- txtProgressBar(min=1, max=prm$sims,
    #                                 style=3, width=50, char="=")
    # } else {
    #   setTxtProgressBar(progressbar, simnum)
    # }
    
    simnum <- seed_curr
    cat(simnum)
    nr_out_single <- rep(NA, 2*length(prm$psi_star_vec))
    
    set.seed(simnum)
    df <- create_sample(prm=prm)
    
    fit_trt_out <- df %>% fit_treatment_models(prm=prm)
    df <- fit_trt_out[[1]]
    prm$trt_models <- fit_trt_out[[2]]
    
    nri_out <- df %>% newton_raphson_piece(prm=prm, print_results=F) 
    #psi_start_vec=rep(0, length(prm$psi_star_vec)))
    (psi_hat_vec <- nri_out[[1]])
    
    (var_hat <- df %>% calculate_variance(psi_hat_vec=psi_hat_vec, prm=prm, print_results=F))
    # trt_models = trt_models))
    #(var_hat_fast <- df %>% calculate_variance_fast(psi_hat_vec=psi_hat_vec, prm=prm))
    
    if(max(diag(var_hat)) < 1) {
      nr_out_single <- c(psi_hat_vec, diag(var_hat))
    }
    
    # if(simnum==prm$sims){
    #   close(progressbar)
    # }
    
    # })
    }
  
   stopCluster(cl)
  
  nr_out <- nr_out_list %>% Reduce(rbind, .)
  results_vec <- calculate_results(nr_out, prm)
  #results_df <- calculate_results_df(nr_out, prm)
  
  return(list(nr_out=nr_out, results_vec=results_vec, results_df=results_df))
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


calculate_results_df <- function(nr_out, prm){
  
  results_df <- tibble(name=character(), mean=numeric(), var=numeric())
  
  for (j in 1:length(prm$psi_star_vec)){
    results_df <- results_df %>% 
      add_row(name="true",
              mean=prm$psi_star_vec[j]) %>%
      add_row(name="mean", 
              mean=mean(nr_out[,j] , na.rm=T),
              var=var(nr_out[,j], na.rm=T)) %>%
      add_row(name="aVar",
              mean=mean(nr_out[,j+length(prm$psi_star_vec)], na.rm=T),
              var=var(nr_out[,j+length(prm$psi_star_vec)], na.rm=T))
  }
  
  #results_df <- results_df %>% column_to_rownames("name")
  return(results_df)
}


print_results_df <- function(results_df, latex=FALSE){
  
  if(latex){
    results_df %>% gt(rowname_col="name") %>%
      fmt_number(decimals=5) %>%
      as_latex() %>% as.character() %>% cat()
  } else {
    results_df %>% gt(rowname_col="name")
  }
  
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



one_factor_bias_analysis <- function(sims=3000, 
                                     n_trgt_vec = c(100, 400, 800, 2000),
                                     psi_1_star = c(log(2))){
  
  prm <- list() 
  
  prm$sims <- sims
  
  prm$sim_label <- "(const_1) -> (const_2)"
  prm$t_a_vec <- c(0, 35)
  prm$expmean <- 50
  prm$beta_1_track <- c(1, 1)
  prm$beta_x_track <- c(0, 0)
  prm$psi_lab <- c("psi_1", "psi_2")
  prm$psi_1_star <- psi_1_star
  prm$psi_x_star <- c()
  prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
  prm$censor_date <- 80
  
  results_df <- tibble(n_trgt=numeric(), psi_hat=numeric(),
                       censor=logical())
  censor_vec <- c(T,F)
  for(censor_val in censor_vec){
    prm$censor <- censor_val
    for(n_trgt in n_trgt_vec){
      print(n_trgt)
      prm$n_trgt <- n_trgt
      
      nr_out_list <- nr_run(prm=prm)
      results_vec <- nr_out_list$results_vec
      
      results_df <- results_df %>% add_row(n_trgt=prm$n_trgt,
                                           censor=prm$censor,
                                           psi_hat=results_vec[1]
      )
      
    }
  }
  
  #mutate(var_type = replace(var_type, var_type=='asymp_var', 'Mean asymptotic variance'),
  #       var_type = replace(var_type, var_type=='asymp_var', sample_var='Sample variance of the sample means')) %>%
  
  results_df %>% 
    ggplot(aes(x=n_trgt, y=psi_hat, line_style=censor)) + 
    geom_line(aes(linetype=censor)) + labs(x='n', y='psi-hat sample mean')
  #rename('Sample variance of the sample means'=sample_var,
  #       'Mean asymptotic variance'=asymp_var) %>%
  #pivot_longer(cols = !n_trgt, names_to = "var_type", values_to = "variance") %>% 
  #ggplot(aes(x=n_trgt, y=variance, colour=var_type, line_style=var_type)) + 
  #geom_line(aes(linetype=var_type)) + labs(x='n') + 
  #theme(legend.title = element_blank())
  
  #var_df %>% pivot_longer(cols = !n_trgt, names_to = "var_type", values_to = "variance") %>%
  #  ggplot(aes(x=n_trgt)) + 
  #  geom_line(aes(y=variance), group_by = var_type)
  #geom_line(aes(y=asymp_var), color = "steelblue", linetype="dashed") +
  #geom_line(aes(y=sample_var), color = "darkred")
  
  
  return(results_df)
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


#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           (1) test block
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##

prm$sim_label <- "(const_1)"
prm$t_a_vec <- c(0)
prm$expmean <- 50
prm$n_trgt <- 400
prm$beta_1_track <- c(1)
prm$beta_x_track <- c(0)
prm$psi_lab <- c("psi_1")
prm$psi_x_star <- c()
prm$sims <- 300
prm$censor_date <- 80
prm$censor_max <- prm$censor_date
prm$trt_mod_list <- list(
  c("1"))


prm$psi_1_star <- c(log(2))
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$censor <- F
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)
prm$censor <- T
nr_out_list <- nr_run(prm=prm)
# print_results_df(nr_out_list$results_df)
print_results(nr_out_list$results_vec, prm)


prm$psi_1_star <- c(log(1))
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$censor <- F
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)
prm$censor <- T
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)


prm$psi_1_star <- c(log(1/2))
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$censor <- F
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)
prm$censor <- T
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)




#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           (1, x1) test block
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##

prm$sim_label <- "(const_1, x_1)"
prm$t_a_vec <- c(0)
prm$expmean <- 50
prm$n_trgt <- 400
prm$beta_1_track <- c(1)
prm$beta_x_track <- c(1)
prm$psi_lab <- c("psi_1", "psi_x1")
prm$sims <- 3000
prm$censor_date <- 80


prm$psi_1_star <- c(log(2))
prm$psi_x_star <- c(log(1.5))
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$censor <- F
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)
prm$censor <- T
nr_out_list <- nr_run(prm=prm)
#print_results_df(nr_out_list$results_df)
print_results(nr_out_list$results_vec, prm)


prm$psi_1_star <- c(log(1))
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$censor <- F
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)
prm$censor <- T
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)


prm$psi_1_star <- c(log(1/2))
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$censor <- F
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)
prm$censor <- T
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)





#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           (1) -> (1) block                            #   ##
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##



prm$sim_label <- "(const_1) -> (const_1)"
prm$t_a_vec <- c(0, 35)
prm$expmean <- 50
prm$n_trgt <- 400
prm$beta_1_track <- c(1, 1)
prm$beta_x_track <- c(0, 0)
prm$psi_lab <- c("psi_1")
prm$sims <- 3000
prm$censor_date <- 80

prm$psi_1_star <- c(log(2))
prm$psi_x_star <- c()
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$censor <- F
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)
prm$censor <- T
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)

prm$psi_1_star <- c(log(1))
prm$psi_x_star <- c()
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$censor <- F
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)
prm$censor <- T
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)


prm$psi_1_star <- c(log(1/2))
prm$psi_x_star <- c()
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$censor <- F
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)
prm$censor <- T
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)


#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           (1, x1) -> (1, x1) block                            #   ##
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##



prm$sim_label <- "(const_1, x_1) -> (const_1, x_1)"
prm$t_a_vec <- c(0, 35)
prm$expmean <- 50
prm$n_trgt <- 400
prm$beta_1_track <- c(1, 1)
prm$beta_x_track <- c(1, 1)
prm$psi_lab <- c("psi_1", "psi_x1")
prm$sims <- 3000
prm$censor_date <- 80


prm$psi_1_star <- c(log(2))
prm$psi_x_star <- c(log(1.5))
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$censor <- F
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)
prm$censor <- T
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)


prm$psi_1_star <- c(log(2))
prm$psi_x_star <- c(log(1/1.5))
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$censor <- F
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)
prm$censor <- T
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)

prm$psi_1_star <- c(log(1/2))
prm$psi_x_star <- c(log(1.5))
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$censor <- F
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)
prm$censor <- T
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)

prm$psi_1_star <- c(log(1/2))
prm$psi_x_star <- c(log(1/1.5))
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$censor <- F
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)
prm$censor <- T
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)



#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           (1) -> (2) block                            #   ##
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##



# 
# prm$sim_label <- "(const_1) -> (const_2)"
# prm$t_a_vec <- c(0, 5)
# prm$expmean <- 50
# prm$n_trgt <- 10000
# prm$beta_1_track <- c(1, 2)
# prm$beta_x_track <- c(0, 0)
# prm$psi_lab <- c("psi_1", "psi_2")
# prm$sims <- 1000
# prm$censor_date <- 80
# prm$trt_mod_list <- list(
#   c("1"), 
#   c("1", "a_0")
# )


prm$sim_label <- "(const_1) -> (const_2)"
prm$t_a_vec <- c(0, 35)
prm$expmean <- 50
prm$n_trgt <- 400
prm$beta_1_track <- c(1, 2)
prm$beta_x_track <- c(0, 0)
prm$psi_lab <- c("psi_1", "psi_2")
prm$sims <- 100
prm$censor_date <- 80
prm$trt_mod_list <- list(
  c("1"), 
  c("1", "a_0")
)



prm$psi_1_star <- c(log(2), log(2))
prm$psi_x_star <- c()
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$censor <- F
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)
prm$censor <- T
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)

prm$psi_1_star <- c(log(1/2), log(1/2))
prm$psi_x_star <- c()
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$censor <- F
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)
prm$censor <- T
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)

prm$psi_1_star <- c(log(2), log(1/2))
prm$psi_x_star <- c()
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$censor <- F
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)
prm$censor <- T
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)

prm$psi_1_star <- c(log(1/2), log(2))
prm$psi_x_star <- c()
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$censor <- F
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)
prm$censor <- T
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)




#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           (1, x1) -> (2, x2) block                            #   ##
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##



prm$sim_label <- "(const_1, x_1) -> (const_2, x_2)"
prm$t_a_vec <- c(0, 35)
prm$expmean <- 50
prm$n_trgt <- 400
prm$beta_1_track <- c(1, 2)
prm$beta_x_track <- c(1, 2)
prm$psi_lab <- c("psi_1", "psi_2", "psi_x1", "psi_x2")
prm$trt_mod_list <- list(
                      c("1", "x"),
                      c("1", "a_0", "x")
                    )
prm$trt_effect_list <- list(
  c("1", "x"),
  c("1", "x")
)
prm$sims <- 3000
prm$censor_date <- 80


prm$psi_1_star <- c(log(2), log(2))
prm$psi_x_star <- c(log(1.5), log(1.5))
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$censor <- F
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)
prm$censor <- T
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)


prm$psi_1_star <- c(log(2), log(1/2))
prm$psi_x_star <- c(log(1.5), log(1/1.5))
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$censor <- F
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)
prm$censor <- T
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)


prm$psi_1_star <- c(log(2))
prm$psi_x_star <- c(log(1/1.5))
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$censor <- F
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)
prm$censor <- T
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)

prm$psi_1_star <- c(log(1/2))
prm$psi_x_star <- c(log(1.5))
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$censor <- F
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)
prm$censor <- T
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)

prm$psi_1_star <- c(log(1/2))
prm$psi_x_star <- c(log(1/1.5))
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$censor <- F
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)
prm$censor <- T
nr_out_list <- nr_run(prm=prm)
print_results(nr_out_list$results_vec, prm)











