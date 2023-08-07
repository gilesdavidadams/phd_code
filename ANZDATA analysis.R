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
library(haven)
library(survival)
library(stats)

source("G-estimation_FINAL.R")


options(dplyr.summarise.inform = FALSE)




#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           DATA CLEANING
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##


#read in data
df_all <- read_dta("C:/Users/gdadams/OneDrive - The University of Melbourne/ANZDATA Analysis/mergedJK_20131018_FINAL.dta")

#choose variables/columns that we want to keep
df_1 <- df_all %>% select(c(id, sequno, firstdate, finaldate, deathdate, deadind,
                            survtime, rxhomeva, "_t", "_t0")) %>%
  select(-c(sequno, firstdate, finaldate, deathdate, survtime)) %>%
  rename(t0 = "_t0", t="_t")


#optional filter for TESTING PURPOSES
#df <- df %>% filter(id %in% c("P7000201", "P7000212"))

#split the dataset at the period cut-off points
df_split <- survSplit(Surv(t0, t, deadind)~., df_1, cut=seq(90,3650,90))
#create a variable indicating which period each row belongs to
df_split <- df_split %>% mutate(period= t0 %/% 90)

#this chunk calculates which treatment is on for the longest in the 90 days
rxhomeva_cols <- c()
rxhomeva_names <- c()
for(i in 0:max(df_split$rxhomeva, na.rm=TRUE)){
  rxhomeva_cols <- c(rxhomeva_cols, paste0("rxhomeva_", i))
  rxhomeva_names <- c(rxhomeva_names, as.name(paste0("rxhomeva_", i)))
  df_split <- df_split %>% 
    mutate("rxhomeva_{i}":=ifelse(rxhomeva==i, t-t0, 0))
}
rxhomeva_sum_cols <- c()
progressbar <- txtProgressBar(min=0, max=max(df_split$rxhomeva, na.rm=TRUE),
                              style=3, width=50, char="=")
for(i in 0:max(df_split$rxhomeva, na.rm=TRUE)){
  setTxtProgressBar(progressbar, i)
  rxhomeva_sum_cols <- c(rxhomeva_sum_cols, paste0("rxhomeva_sum_", i))
  col_name <- as.name(paste0("rxhomeva_", i))
  df_split_grp <- df_split %>% group_by(id, period) %>%
    summarise("rxhomeva_sum_{i}":= sum({{col_name}}))
  
  df_split <- df_split %>% 
    left_join(df_split_grp,
              by=c("id", "period")) %>%
    select(-c({{col_name}}))
}
close(progressbar)



df_split <- df_split %>% 
  mutate(rxhomeva_max = 
           colnames(df_split[, rxhomeva_sum_cols])[max.col(df_split[, rxhomeva_sum_cols])]) %>%
  select(-all_of(rxhomeva_sum_cols)) %>%
  mutate(rxhomeva_90=as.integer(substring(rxhomeva_max, 14)))

df_split <- df_split %>% group_by(id) %>% mutate(ti=max(t)) %>% ungroup()

df_90 <- df_split %>%
  group_by(id, period) %>% 
  summarise(rxhomeva_90=max(rxhomeva_90),
            deadind=max(deadind),
            ti=max(ti),
            tk=sum(t-t0)) %>%
  ungroup()

df_90 <- df_90 %>% group_by(id) %>% 
  mutate(dies=max(deadind)) %>% 
  ungroup() %>% select(-c(deadind, tk)) %>%
  rename(ai = rxhomeva_90) %>% 
  mutate(ai_NA = as.numeric(is.na(ai)))

df_90_NA <- df_90 %>% group_by(id) %>%
  summarise(has_NA = max(ai_NA))

df_90 <- df_90 %>% left_join(df_90_NA,
                             by="id")

df_90 <- df_90 %>% filter(has_NA != 1)

df_90 <- df_90 %>% rename(ai_finest=ai)


df_home_vs_facility_long <- df_90 %>% mutate(ai_HvF = ifelse(ai_finest %in% c(0,1,2), 1, 0))
df_3_vs_not <- df_90 %>% mutate(ai_HvF = ifelse(ai_finest %in% c(3), 1, 0))

df_HvF_wide <- df_home_vs_facility_long %>% select(-c(has_NA, ai_NA, ai_finest)) %>%
  pivot_wider(names_from = period,
              names_prefix = "a_",
              values_from = ai_HvF)  

df_3_wide <- df_3_vs_not %>% select(-c(has_NA, ai_NA, ai_finest)) %>%
  pivot_wider(names_from = period,
              names_prefix = "a_",
              values_from = ai_HvF) 




#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           DATA ANALYSIS
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##




prm <- list()
prm$sim_label <- "(const_1)"
prm$t_a_vec <- seq(0, 90*33, 90)
#prm$expmean <- 50
#prm$n_trgt <- 400
prm$beta_1_track <- rep(1,34)
prm$beta_x_track <- rep(0, 34)
prm$psi_lab <- c("psi_1")
#prm$psi_x_star <- c()
#prm$sims <- 3000
prm$censor_date <- 33*90
prm$psi_1_star <- c(log(2))
prm$psi_x_star <- c()
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$trt_mod_list <- sapply(0:33, function(k) {
  if(k==0){ 
    return(c("1"))
  } else {
    return( c("1", paste0("a_", k)) )
  }
})

prm$censor <- T


#df<- df_HvF_wide %>% fit_treatment_models(prm=prm)
df <- df_3_wide


df <- df %>% fit_treatment_models(prm=prm)
df <- df %>% mutate(C_i = prm$censor_date, x = 0)



nr_out <- df %>% newton_raphson_grad(prm=prm, 
                                     psi_start_vec=rep(0, length(prm$psi_star_vec)))
(psi_hat_vec <- nr_out[[1]])
(var_hat <- df %>% calculate_variance(psi_hat_vec=psi_hat_vec, prm=prm))

for(i in 1:33){
  col_name <- as.name(paste0("a_", i))
  df <- df %>% mutate(has_as = has_as + 1 - as.integer(is.na({{col_name}})) )
}

df <- df %>% mutate(na_cutoff = 90 %/% has_as)
df <- df %>% mutate(is_good = (ti <= 90 * (has_as + 1)))
