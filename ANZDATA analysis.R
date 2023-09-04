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
library(labelled)
library(janitor)

source("G-estimation_FINAL.R")


options(dplyr.summarise.inform = FALSE)




#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           DATA CLEANING
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##


#read in data
df_all <- read_dta("C:/Users/gdadams/OneDrive - The University of Melbourne/ANZDATA Analysis/mergedJK_20131018_FINAL.dta")

#choose variables/columns that we want to keep
df_reduced <- df_all %>% select(c(id, sequno, firstdate, finaldate, deathdate, deadind,
                            survtime, rxhomeva, "_t", "_t0",
                            age, sex, racereduced, sercreat, cig,
                            referral3, rxhome)
                            ) %>%
  #select(-c(sequno, firstdate, finaldate, deathdate, survtime)) %>%
  rename(t0 = "_t0", t="_t", smoke=cig, latereferral=referral3)



#optional filter for TESTING PURPOSES
#df <- df %>% filter(id %in% c("P7000201", "P7000212"))

#split the dataset at the period cut-off points
df_split <- survSplit(Surv(t0, t, deadind)~., df_reduced, cut=seq(90,3650,90))
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

df_90 <- df_90 %>%
          left_join(df_90_NA, by="id") %>%
          filter(has_NA != 1) %>% 
          rename(ai_finest=ai) %>%
          mutate(is.PD = (ai_finest == 6))

df_90_PD <- df_90 %>% group_by(id) %>%
      summarise(has.PD = max(as.numeric(is.PD)))

df_90 <- df_90 %>%
          left_join(df_90_PD, by="id") %>%
          filter(has.PD != 1)

#df_90_wide <- df_90 %>% select(-c(ai_NA, has_NA, is.PD, has.PD)) %>%
#              pivot_wider(names_from = period,
#                    names_prefix = "a_",
#                    values_from = ai_finest)

df_90_home <- df_90 %>% mutate(HvF = ifelse(ai_finest %in% c(0,1,2), 1, 0)) %>%
                        select(c(id, period, HvF)) %>%
                        pivot_wider(names_from = period,
                        names_prefix = "HvF_",
                        values_from = HvF)

df_CVC_vs_AVFG_long <- df_90 %>% mutate(ai = ifelse(ai_finest %in% c(2, 5), 1, 0))


# df_home_vs_facility_long <- df_90 %>% mutate(ai_HvF = ifelse(ai_finest %in% c(0,1,2), 1, 0))
# df_3_vs_not <- df_90 %>% mutate(ai_3vn = ifelse(ai_finest %in% c(3), 1, 0))
# df_home_vs_not <- df_90 %>% mutate(ai_Hvn = ifelse(ai_finest %in% c(0,1,2), 1, 0))
# df_facility_vs_not <- df_90 %>% mutate(ai_HvF = ifelse(ai_finest %in% c(3,4,5), 1, 0))

#df_HvF_wide <- df_home_vs_facility_long %>% select(-c(has_NA, ai_NA, ai_finest)) %>%
#  pivot_wider(names_from = period,
#              names_prefix = "a_",
#              values_from = ai_HvF)


# df_home_vs_not_wide <- df_home_vs_not %>% select(-c(has_NA, ai_NA, ai_finest)) %>%
#                         pivot_wider(names_from = period,
#                                       names_prefix = "a_",
#                                       values_from = ai_Hvn)
# 
# df_3_wide <- df_3_vs_not %>% select(-c(has_NA, ai_NA, ai_finest)) %>%
#   pivot_wider(names_from = period,
#               names_prefix = "a_",
#               values_from = ai_3vn) 


df_CVC_vs_AVFG_wide <- df_CVC_vs_AVFG_long %>% 
                        select(c(id, period, ti, dies, ai)) %>%
                        pivot_wider(names_from = period,
                                    names_prefix = "a_",
                                    values_from = ai)


# Collect baseline variable
t0_var <- as.name("_t0")
df_baseline <- df_all %>% 
                group_by(id) %>% 
                summarise(age=max(age),
                           sex=max(sex),
                           race=max(racereduced),
                           firstdate=max(firstdate),
                           lastperiod=max({{t0_var}} %/% 90),
                           sercreat=max(sercreat),
                           smoke=max(cig),
                           latereferral=max(referral3),
                           cad=max(BEcoronary),
                           lung=max(BElung),
                           pvd=max(BEpvd),
                           diabetes=max(diabB),
                           bmi=max(bmi),
                           at_home=max(rxhome)
                )
df_baseline$bmi_cat <- cut(df_baseline$bmi, c(0, 20, 25, 30, 100))


 #Append baseline variables to wide datasets



#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           DATA ANALYSIS
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##

df_augment <- function(df, prm){
  
  df <- df %>% left_join(df_baseline, by="id")
  df <- df %>% mutate(sex=as.factor(sex),
                      race=factor(race, labels = c("NZ maori/pacific" = 1,
                                                   "Aus Indigenous" = 2,
                                                   "Asian" = 3,
                                                   "White/other" = 4)),
                      firstyear = as.factor(format(firstdate,"%Y")),
                      smoke = as.factor(smoke)
  )
  
  df <- df %>% filter(!is.na(sercreat)) %>% 
    filter(!latereferral=="") %>%
    filter(smoke!=4) %>%
    filter(!is.na(bmi_cat))
 
  df <- df %>% mutate(x = 0) %>%
    mutate(C_i = difftime(prm$censor_date, firstdate,  "days"))
  
  df <- df %>% left_join(df_90_home, by="id")
  
  
  
  #fit_trt_out <- df %>% fit_treatment_models(prm=prm)
  #df <- fit_trt_out[[1]]
  #trt_models <- fit_trt_out[[2]]
  
  return(df)
  #return(list(df=df, trt_models=trt_models))
}



k_max_chosen <- 33
k_max_actual <- 33

prm <- list()
prm$sim_label <- "(const_1)"
prm$t_a_vec <- c(seq(0, 90*(k_max_chosen), 90))
#prm$expmean <- 50
#prm$n_trgt <- 400
prm$beta_1_track <- rep(1, k_max_chosen+1)
prm$beta_x_track <- rep(0, k_max_chosen+1)
prm$psi_lab <- c("psi_1")
#prm$psi_x_star <- c()
#prm$sims <- 3000
prm$censor_date <- as.Date("2011-12-31")
prm$trt_mod_list <- sapply(0:33, function(k) {
  if(k==0){ 
    return(c("1", "age", "sex", "race", "smoke", "firstyear"#, "sercreat"
    ))
  } else if (k %in% 30:33) {
    return( c("1"))
  } else if (k %in% 26:33) {
    return( c("1", paste0("a_", k-1)))
  } else {
    return( c("1", paste0("a_", k-1), 
              "age", "sex", "race", "smoke"))
  }
})
prm$censor <- T

df_list <- df_CVC_vs_AVFG_wide %>% df_augment(prm=prm)
df <- df_list$df
trt_models <- df_list$trt_models


# All periods together, single psi
nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                     psi_start_vec=rep(0, max(prm$beta_1_track) + max(prm$beta_x_track)))
(psi_hat_vec <- nr_out[[1]])


first_part <- 4
prm$beta_1_track <- c(rep(1, first_part), rep(2, 34-first_part))
nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                     psi_start_vec=rep(0, max(prm$beta_1_track) + max(prm$beta_x_track)))
(psi_hat_vec <- nr_out[[1]])





first_part <- c(rep(1,4), rep(2,4))#, rep(3,2))
prm$beta_1_track <- c(first_part, rep(max(first_part)+1, 34-length(first_part)))
#  c(rep(1, 4), 2, rep(3, 34-(4+1)))

nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                     psi_start_vec=rep(0, max(prm$beta_1_track) + max(prm$beta_x_track)))
(psi_hat_vec <- nr_out[[1]])




first_part <- c(rep(1,4), rep(2,4), rep(3,4))
prm$beta_1_track <- c(first_part, rep(max(first_part)+1, 34-length(first_part)))
#  c(rep(1, 4), 2, rep(3, 34-(4+1)))

nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                     psi_start_vec=rep(0, max(prm$beta_1_track) + max(prm$beta_x_track)))
(psi_hat_vec <- nr_out[[1]])



first_part <- c(rep(1,4), rep(2,4), rep(3,4), rep(4,4))
prm$beta_1_track <- c(first_part, rep(max(first_part)+1, 34-length(first_part)))
#  c(rep(1, 4), 2, rep(3, 34-(4+1)))

nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                     psi_start_vec=rep(0, max(prm$beta_1_track) + max(prm$beta_x_track)))
(psi_hat_vec <- nr_out[[1]])



first_part <- c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4))
prm$beta_1_track <- c(first_part, rep(max(first_part)+1, 34-length(first_part)))
#  c(rep(1, 4), 2, rep(3, 34-(4+1)))

nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                     psi_start_vec=rep(0, max(prm$beta_1_track) + max(prm$beta_x_track)))
(psi_hat_vec <- nr_out[[1]])



first_part <- c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4))
prm$beta_1_track <- c(first_part, rep(max(first_part)+1, 34-length(first_part)))
#  c(rep(1, 4), 2, rep(3, 34-(4+1)))

nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                     psi_start_vec=rep(0, max(prm$beta_1_track) + max(prm$beta_x_track)))
(psi_hat_vec <- nr_out[[1]])



(VCV <- df %>% calculate_variance(prm, psi_hat_vec, trt_models))
nr_out_list <- df_3_wide_run(prm)





#df<- df_HvF_wide %>% fit_treatment_models(prm=prm)
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           HOME VS NOT                           #   ##
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
df <- df_home_vs_not_wide
df <- df %>% left_join(df_baseline, by="id")
df <- df %>% mutate(sex=as.factor(sex),
                    race=factor(race, labels = c("NZ maori/pacific" = 1,
                                                 "Aus Indigenous" = 2,
                                                 "Asian" = 3,
                                                 "White/other" = 4)),
                    firstyear = as.factor(format(firstdate,"%Y")),
                    smoke = as.factor(smoke)
)

df <- df %>% filter(!is.na(sercreat))
prm$trt_mod_list <- sapply(0:33, function(k) {
  if(k==0){ 
    return(c("1"))
  } else if (k %in% c(30, 31, 32, 33)) {
    return( c("1"))
  } else {
    return( c("1", paste0("a_", k-1)))
  }
})

fit_trt_out <- df %>% fit_treatment_models(prm=prm)
df <- fit_trt_out[[1]]
trt_models <- fit_trt_out[[2]]
prm$beta_1_track <- rep(1,34)
df <- df %>% mutate(x = 0) %>%
  mutate(C_i = difftime(prm$censor_date, firstdate,  "days"))
nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                     psi_start_vec=rep(0, max(prm$beta_1_track) + max(prm$beta_x_track)))
(psi_hat_vec <- nr_out[[1]])









prm <- list()
prm$sim_label <- "(const_1)"
prm$t_a_vec <- c(seq(0, 90*33, 90))
#prm$expmean <- 50
#prm$n_trgt <- 400
prm$beta_1_track <- rep(1, (33+1))
prm$beta_x_track <- rep(0, (33+1))
prm$psi_lab <- paste0("psi_", 1:(33+1))
#prm$psi_x_star <- c()
#prm$sims <- 3000
prm$censor_date <- as.Date("2011-12-31")
prm$psi_1_star <- rep(log(2))
prm$psi_x_star <- c()
prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
prm$trt_mod_list <- sapply(0:33, function(k) {
  if(k==0){ 
    return(c("1", "age", "sex", "race", "smoke", "firstyear"#, "sercreat"
             ))
  } else if (k %in% c(30, 31, 32, 33)) {
    return( c("1", paste0("a_", k-1)))
  } else {
    return( c("1", paste0("a_", k-1), "age", "sex", "race"))
  }
})
prm$censor <- T





q <- lapply(1:10, function(k){
  cat("k ", k, "\n")
  prm$beta_1_track <- c(1:k, rep(k+1, max(c(33+1-k,0))))
  prm$psi_1_star <- rep(log(2), k+1)
  prm$psi_star_vec <- c(prm$psi_1_star, prm$psi_x_star)
  
  nr_out <- df %>% newton_raphson_grad(prm=prm, 
                                       psi_start_vec=rep(0, length(prm$psi_star_vec)))
  cat("\n", nr_out[[1]], "\n")
  return(nr_out[[1]])
})

tibble(a = 0:7) %>% cbind(
lapply(1:20, function(k){
  thevar <- as.name(paste0("a_", k))
  df_90_wide %>% tabyl({{thevar}}) %>% select(c(n))
  }) %>% Reduce(cbind, .)) 


(df_all$rxhomeva %>% attributes())$labels %>% str()





#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           TESTING MODELS FOR BASELINE COVARIATES
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##



df <- df_CVC_vs_AVFG_wide
df <- df %>% left_join(df_baseline, by="id")
df <- df %>% mutate(sex=as.factor(sex),
                    race=factor(race, labels = c("NZ maori/pacific" = 1,
                                                 "Aus Indigenous" = 2,
                                                 "Asian" = 3,
                                                 "White/other" = 4)),
                    firstyear = as.factor(format(firstdate,"%Y")),
                    smoke = as.factor(smoke)
)
df <- df %>% filter(!is.na(sercreat)) %>% 
  filter(!latereferral=="") %>%
  filter(smoke!=4) %>%
  filter(!is.na(bmi_cat))

# These models converge up to an including k=25 (k=26 fails)

##### TRY MODEL UP TO K=25


prm <- list()
prm$sim_label <- "(const_1)"
prm$t_a_vec <- c(seq(0, 90*25, 90))
prm$beta_1_track <- rep(1, 26)
prm$beta_x_track <- rep(0, 26)
prm$psi_lab <- c("psi_1")
prm$censor_date <- as.Date("2003-10-01") + 27*90
prm$trt_mod_list <- sapply(0:25, function(k) {
  mdl_k <- c("1", "age", "sex", "race", "smoke", "firstyear",
             "sercreat", "latereferral", "cad", "lung", "pvd",
             "diabetes", "bmi_cat", "at_home")
  if(k > 0){ mdl_k <- c(mdl_k, paste0("a_", k-1))}
  return(mdl_k)
})
prm$censor <- T

df <- df %>% filter(firstdate < prm$censor_date)

fit_trt_out <- df %>% fit_treatment_models(prm=prm)
df <- fit_trt_out[[1]]
trt_models <- fit_trt_out[[2]]
df <- df %>% mutate(x = 0) %>%
  mutate(C_i = difftime(prm$censor_date, firstdate,  "days"))

nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                     psi_start_vec=rep(0, max(prm$beta_1_track) + max(prm$beta_x_track)))
(psi_hat_vec <- nr_out[[1]])



#prm$beta_1_track <- (c(rep(1,4), rep(2, 26-4)))
#nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
#                                     psi_start_vec=rep(0, max(prm$beta_1_track) + max(prm$beta_x_track)))
#(psi_hat_vec <- nr_out[[1]])


prm$beta_1_track <- c(rep(1,4), rep(2,4), rep(3,4),
                      rep(4,4), rep(5,4), rep(6,6))
nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                     psi_start_vec=rep(0, max(prm$beta_1_track) + max(prm$beta_x_track)))
(psi_hat_vec <- nr_out[[1]])
(VCV <- df %>% calculate_variance(prm, psi_hat_vec, trt_models))




#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           A vs C
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##




# df_CvA <- df_CvA %>% left_join(df_baseline, by="id")
# df_CvA <- df_CvA %>% mutate(sex=as.factor(sex),
#                             race=factor(race, labels = c("NZ maori/pacific" = 1,
#                                                          "Aus Indigenous" = 2,
#                                                          "Asian" = 3,
#                                                          "White/other" = 4)),
#                             firstyear = as.factor(format(firstdate,"%Y")),
#                             smoke = as.factor(smoke)
# )
# df_CvA <- df_CvA %>% filter(!is.na(sercreat)) %>% 
#   filter(!latereferral=="") %>%
#   filter(smoke!=4) %>%
#   filter(!is.na(bmi_cat))

##### 6 "year" model, censored at end of period 25
prm <- list()
prm$sim_label <- "(const_1)"
prm$t_a_vec <- c(seq(0, 90*25, 90))
prm$beta_1_track <- rep(1, 26)
prm$beta_x_track <- rep(0, 26)
prm$psi_lab <- c("psi_1")
prm$censor_date <- as.Date("2003-10-01") + 27*90
prm$censor <- T

prm$beta_1_track <- c(rep(1,4), rep(2,4), rep(3,4),
                      rep(4,4), rep(5,4), rep(6,6))

# df_CvA <- df_CvA %>% filter(firstdate < prm$censor_date) %>%
#              mutate(x = 0) %>%
#             mutate(C_i = difftime(prm$censor_date, firstdate,  "days"))


#full treatment model
prm$trt_mod_list <- sapply(0:25, function(k) {
  mdl_k <- c("1", "age", "sex", "race", "smoke", "firstyear",
             "sercreat", "latereferral", "cad", "lung", "pvd",
             "diabetes", "bmi_cat", "at_home")
  if(k > 0){ mdl_k <- c(mdl_k, paste0("a_", k-1))}
  return(mdl_k)
})

df_CvA_aug <- df_CVC_vs_AVFG_wide %>% df_augment(prm=prm)
fit_trt_out <- df_CvA_aug %>% fit_treatment_models(prm=prm)
df <- fit_trt_out[[1]]
trt_models <- fit_trt_out[[2]]

nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                     psi_start_vec=rep(0, max(prm$beta_1_track) + max(prm$beta_x_track)))
(psi_hat_vec <- nr_out[[1]])
(VCV <- df %>% calculate_variance(prm, psi_hat_vec, trt_models))


#reduced treatment model
prm$trt_mod_list <- sapply(0:25, function(k) {
  mdl_k <- c("1", #"age", 
             "sex", "race", #"smoke", #"firstyear",
             "sercreat", "latereferral", "cad", "lung", "pvd",
             "diabetes", "bmi_cat", "at_home")
  if(k > 0){ mdl_k <- c(mdl_k, paste0("a_", k-1))}
  return(mdl_k)
})

fit_trt_out <- df_CvA %>% fit_treatment_models(prm=prm)
df <- fit_trt_out[[1]]
trt_models <- fit_trt_out[[2]]

nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                     psi_start_vec=rep(0, max(prm$beta_1_track) + max(prm$beta_x_track)))
(psi_hat_vec <- nr_out[[1]])
(VCV <- df %>% calculate_variance(prm, psi_hat_vec, trt_models))




#baseline model then only past treatment


prm$trt_mod_list <- sapply(0:25, function(k) {
  if(k==0){
    mdl_k <- c("1", "age", "sex", "race", "smoke", "firstyear",
               "sercreat", "latereferral", "cad", "lung", "pvd",
               "diabetes", "bmi_cat", "at_home")
  } else {  
    mdl_k <- c("1", paste0("a_", k-1))
  }
  return(mdl_k)
})

fit_trt_out <- df_CvA %>% fit_treatment_models(prm=prm)
df <- fit_trt_out[[1]]
trt_models <- fit_trt_out[[2]]

nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                     psi_start_vec=rep(0, max(prm$beta_1_track) + max(prm$beta_x_track)))
(psi_hat_vec <- nr_out[[1]])
(VCV <- df %>% calculate_variance(prm, psi_hat_vec, trt_models))


#full treatment model with varying windows
prm$trt_mod_list <- sapply(0:25, function(k) {
  mdl_k <- c("1", "age", "sex", "race", "smoke", "firstyear",
             "sercreat", "latereferral", "cad", "lung", "pvd",
             "diabetes", "bmi_cat", "at_home")
             #paste0("HvF_", k))
  
  if(k > 0){ mdl_k <- c(mdl_k, paste0("a_", k-1))}
  return(mdl_k)
})

prm$beta_1_track <- c(rep(1,5), rep(2,4), rep(3,4),
                      rep(4,4), rep(5,9))

df_CvA_aug <- df_CVC_vs_AVFG_wide %>% 
                df_augment(prm=prm) %>%
                filter(firstdate < prm$censor_date)
fit_trt_out <- df_CvA_aug %>% fit_treatment_models(prm=prm)
df <- fit_trt_out[[1]]
trt_models <- fit_trt_out[[2]]



nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                     psi_start_vec=rep(0, max(prm$beta_1_track) + max(prm$beta_x_track)))
(psi_hat_vec <- nr_out[[1]])
(VCV <- df %>% calculate_variance(prm, psi_hat_vec, trt_models))

prm$beta_1_track <- c(rep(1,4), rep(2,4), rep(3,4),
                      rep(4,4), rep(5,10))

nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                     psi_start_vec=rep(0, max(prm$beta_1_track) + max(prm$beta_x_track)))
(psi_hat_vec <- nr_out[[1]])
(VCV <- df %>% calculate_variance(prm, psi_hat_vec, trt_models))


MA_out_5 <- lapply(0:4,
          function(k){
            prm$beta_1_track <- c(rep(1,4+k), rep(2,4), rep(3,4),
                                  rep(4,4), rep(5,10-k))
            
            nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                                 psi_start_vec=rep(0, max(prm$beta_1_track) + max(prm$beta_x_track)))
            (psi_hat_vec <- nr_out[[1]])
            VCV <- NA
            #(VCV <- df %>% calculate_variance(prm, psi_hat_vec, trt_models))
            
            return(list(psi_hat_vec=psi_hat_vec,
                   VCV=VCV))
          })

MA_out_4 <- lapply(0:8,
                 function(k){
                   prm$beta_1_track <- c(rep(1,4+k), rep(2,4), rep(3,4),
                                         rep(4,14-k))
                   
                   nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                                        psi_start_vec=rep(0, max(prm$beta_1_track) + max(prm$beta_x_track)))
                   (psi_hat_vec <- nr_out[[1]])
                   #VCV <- NA
                   (VCV <- df %>% calculate_variance(prm, psi_hat_vec, trt_models))
                   
                   return(list(psi_hat_vec=psi_hat_vec,
                               VCV=VCV))
                 })



# Seeing if I can estimate shorter periods by starting at a previous psi_hat
offset=-1
prm$beta_1_track <- c(rep(1,4+offset), rep(2,4), rep(3,4),
                                          rep(4,14-offset))

nr_out_m1 <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                     psi_start_vec=c(-0.8602233, -0.4174207, -0.3215344, -0.1973486),
                                     print_results=T
                                     )

offset=-2
prm$beta_1_track <- c(rep(1,4+offset), rep(2,4), rep(3,4),
                      rep(4,14-offset))
nr_out_m2 <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                     psi_start_vec=nr_out_m1[[1]],
                                     print_results=T
)

offset=-3
prm$beta_1_track <- c(rep(1,4+offset), rep(2,4), rep(3,4),
                      rep(4,14-offset))
nr_out_m3 <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                        psi_start_vec=nr_out_m2[[1]],
                                        print_results=T
)

# Building on previous results to extend size of psi_hat

# Reducing the 1st period
beta_1_start <- c(1, 2, rep(3,4), rep(4,4))
prm$beta_1_track <- c(beta_1_start, rep(max(beta_1_start)+1, 26-length(beta_1_start)))
psi_hat_vec <- nr_out_m3[[1]]
psi_start <- c(psi_hat_vec[1], psi_hat_vec[1], psi_hat_vec[2:length(psi_hat_vec)])

nr_out_1234444etc <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                                psi_start_vec=psi_start,
                                                print_results=T)

beta_1_start <- c(1, 2, 3, rep(4,4), rep(5,4))
prm$beta_1_track <- c(beta_1_start, rep(max(beta_1_start)+1, 26-length(beta_1_start)))
psi_hat_vec <- nr_out_123333etc[[1]]
psi_start <- c(psi_hat_vec[1:2], psi_hat_vec[2], psi_hat_vec[3:length(psi_hat_vec)])

nr_out_1234444etc <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                               psi_start_vec=psi_start,
                                               print_results=T)

beta_1_start <- c(1, 2, 3, 4, rep(5,4), rep(6,4))
prm$beta_1_track <- c(beta_1_start, rep(max(beta_1_start)+1, 26-length(beta_1_start)))
psi_hat_vec <- nr_out_1234444etc[[1]]
psi_start <- c(psi_hat_vec[1:3], psi_hat_vec[3], psi_hat_vec[4:length(psi_hat_vec)])

nr_out_12345555etc <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                                psi_start_vec=psi_start,
                                                print_results=T)
psi_hat_7 <- nr_out_12345555etc[[1]]
(VCV_7 <- df %>% calculate_variance(prm, psi_hat_vec, trt_models))
psi_hat_7 %>% rbind(sapply(1:7, function(k){sqrt(VCV_7[k,k])}))


beta_1_start <- c(1:5, rep(6,3), rep(7,4))
prm$beta_1_track <- c(beta_1_start, rep(max(beta_1_start)+1, 26-length(beta_1_start)))
psi_hat_vec <- nr_out_7[[1]]
psi_start <- c(psi_hat_vec[1:4], psi_hat_vec[4], psi_hat_vec[5:length(psi_hat_vec)])
nr_out_8 <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                                 psi_start_vec=psi_start,
                                                 print_results=T)
psi_hat_7 <- nr_out_12345555etc[[1]]
(VCV_7 <- df %>% calculate_variance(prm, psi_hat_vec, trt_models))
psi_hat_7 %>% rbind(sapply(1:7, function(k){sqrt(VCV_7[k,k])}))


beta_1_start <- c(1:5, rep(6,3), rep(7,4))
prm$beta_1_track <- c(beta_1_start, rep(max(beta_1_start)+1, 26-length(beta_1_start)))
psi_hat_vec <- nr_out_7[[1]]
psi_start <- c(psi_hat_vec[1:4], psi_hat_vec[4], psi_hat_vec[5:length(psi_hat_vec)])
nr_out_8 <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                       psi_start_vec=psi_start,
                                       print_results=T)
psi_hat_7 <- nr_out_12345555etc[[1]]
(VCV_7 <- df %>% calculate_variance(prm, psi_hat_vec, trt_models))
psi_hat_7 %>% rbind(sapply(1:7, function(k){sqrt(VCV_7[k,k])}))




nr_out_8_11 <- lapply(1:4, 
      function(i){
        k_max <- 4+i
        beta_1_start <- c(1:k_max, rep(k_max+1, 4-i), rep(k_max+2, 4))
        prm$beta_1_track <- c(beta_1_start, rep(max(beta_1_start)+1, 26-length(beta_1_start)))
        psi_start <- c(psi_hat_7[1:4], rep(psi_hat_7[4], i), psi_hat_7[4:length(psi_hat_7)])
        nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                             psi_start_vec=psi_start,
                                             print_results=T)
        return(nr_out[[1]])
})

