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
library(labelled)
library(DescTools)

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

MA_4_1st <- tibble(period=3:11,
                   m1_psi_1=-sapply(1:9,
                          function(k){
                            MA_out_4[[k]]$psi_hat_vec[[1]]
                          }),
                   stdev=sapply(1:9,
                          function(k){
                            MA_out_4[[k]]$VCV[1,1] %>% sqrt()
                          })
                   )

MA_4_1st %>%
  ggplot(aes(x=period, y=psi_1)) +
  geom_line() +
  geom_errorbar(aes(ymin=psi_1-2*stdev, ymax=psi_1+2*stdev),
                width=.2, position=position_dodge(0.05)) +
  labs(title="Period 1 to x with 2*se CI's", 
       y = "-1*psi on the log scale", 
       x = "Final period of group")



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
psi_hat_8 <- nr_out_8[[1]]
#(VCV_7 <- df %>% calculate_variance(prm, psi_hat_vec, trt_models))
#psi_hat_7 %>% rbind(sapply(1:7, function(k){sqrt(VCV_7[k,k])}))



beta_1_start <- c(1:6, rep(6,3), rep(7,4))
prm$beta_1_track <- c(beta_1_start, rep(max(beta_1_start)+1, 26-length(beta_1_start)))
psi_hat_vec <- nr_out_7[[1]]
psi_start <- c(psi_hat_vec[1:4], psi_hat_vec[4], psi_hat_vec[5:length(psi_hat_vec)])
nr_out_8 <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                       psi_start_vec=psi_start,
                                       print_results=T)
psi_hat_8 <- nr_out_8[[1]]
#(VCV_7 <- df %>% calculate_variance(prm, psi_hat_vec, trt_models))
#psi_hat_7 %>% rbind(sapply(1:7, function(k){sqrt(VCV_7[k,k])}))


nr_out_9_12 <- lapply(1:4, 
      function(i){
        k_max <- 5+i
        beta_1_start <- c(1:k_max, rep(k_max+1, 5-i), rep(k_max+2, 4))
        prm$beta_1_track <- c(beta_1_start, rep(max(beta_1_start)+1, 26-length(beta_1_start)))
        psi_start <- c(psi_hat_8[1:5], rep(psi_hat_8[5], i), psi_hat_8[6:length(psi_hat_8)])
        nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                             psi_start_vec=psi_start,
                                             print_results=T)
        return(nr_out[[1]])
})
psi_hat_12 <- nr_out_9_12[[4]]


psi_hat_old <- psi_hat_12
nr_out_13_16 <- lapply(1:4, 
      function(i){
        k_max <- 10+i
        beta_1_start <- c(1:k_max, rep(k_max+1, 5-i))
        prm$beta_1_track <- c(beta_1_start, rep(max(beta_1_start)+1, 26-length(beta_1_start)))
        psi_start <- c(psi_hat_old[1:10], rep(psi_hat_old[10], i), psi_hat_old[11:length(psi_hat_old)])
        nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                             psi_start_vec=psi_start,
                                             print_results=T)
        return(nr_out[[1]])
      })

psi_hat_13 <- nr_out_13_16[[1]]
beta_1_start <- c(1:11, rep(12, 4))
prm$beta_1_track <- c(beta_1_start, rep(max(beta_1_start)+1, 26-length(beta_1_start)))
(VCV_13 <- df %>% calculate_variance(prm, psi_hat_13, trt_models))
#psi_hat_7 %>% rbind(sapply(1:7, function(k){sqrt(VCV_7[k,k])}))
psi_VCV_13 <- psi_hat_13 %>% rbind(sapply(1:13, function(k){VCV_13[k,k]}))
prm$beta_1_track

stdev_13 <- sapply(1:13, function(k){VCV_13[k,k] %>% sqrt()})

df_13 <- tibble(period=1:13, psi=-psi_hat_13, 
                stdev=)
df_13 %>%
  ggplot(aes(x=period, y=psi)) +
  geom_line() +
  geom_errorbar(aes(ymin=psi-2*stdev, ymax=psi+2*stdev),
                width=.2, position=position_dodge(0.05)) +
  labs(title="Periods 1-11 are single periods, then groups with 2*se CI's", 
       y = "psi on the log scale", 
       x = "Period (1-11 single, then groups)")
                                   


#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           Building a model that is a sequence of
#           single periods followed by one big period
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##

prm$beta_1_track <- rep(1,26)
nr_out_inc_list <- list()
for(i in 0:13){
  if(i > 0){
    prm$beta_1_track <- c(1:i, rep(i+1, 26-i))
  } else {
    prm$beta_1_track <- rep(1,26)
  }
  
  if(i == 0){
    psi_start <- c(0)
  } else if(i == 1){
    psi_start <- c(psi_hat_old, psi_hat_old)
  } else {
    psi_start <- c(psi_hat_old[1:(i-1)], 
                   psi_hat_old[i-1],
                   psi_hat_old[i:length(psi_hat_old)])

  }
  
  nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                       psi_start_vec=psi_start,
                                       print_results=T)
  psi_hat_old <- nr_out[[1]]
  nr_out_inc_list <- nr_out_inc_list %>% append(., list(psi_hat_old))
}


nr_out_piece_list <- list()
#for(i in 0:13){
for(i in 14:20){
  if(i > 0){
    prm$beta_1_track <- c(1:i, rep(i+1, 26-i))
  } else {
    prm$beta_1_track <- rep(1,26)
  }
  
  if(i == 0){
    psi_start <- c(0)
  } else if(i == 1){
    psi_start <- c(psi_hat_old, psi_hat_old)
  } else {
    psi_start <- c(psi_hat_old[1:(i-1)], 
                   psi_hat_old[i-1],
                   psi_hat_old[i:length(psi_hat_old)])
    
  }
  
  nr_out <- df %>% newton_raphson_piece(prm=prm, tol=0.01,
                                       psi_start_vec=psi_start,
                                       print_results=T)
  psi_hat_old <- nr_out[[1]]
  nr_out_piece_list <- nr_out_piece_list %>% append(., list(psi_hat_old))
}

df_indiv <- tibble(period=0:25)
for(i in 1:20){
  psi_hat <- nr_out_piece_list[[i]]
  
  psi_by_period <- sapply(1:26, 
                    function(k){
                      if(k < length(psi_hat)){
                        return(-psi_hat[k])
                      } else {
                        return(-psi_hat[length(psi_hat)])
                      }
                    })
  
  #ste_by_period <- sapply(1:9, 
  #                        function(k){
  #                          ste_vec[[prm$beta_1_track[[k]]]]  
  #                        })
  psi_col <- as.name(paste0("itr_", i))
  #ste_col <- as.name(paste0("itr_", i, "_ste"))
  df_indiv <- df_indiv %>% add_column("{psi_col}" := psi_by_period)
  
  #nr_out_MA_list <- nr_out_MA_list %>% append(., list(psi_hat_old))
}
df_indiv %>% pivot_longer(cols = !period) %>% 
  ggplot(aes(x=period, y=value, group=name)) +
  geom_point() + geom_line()

# Last period that worked went up to period 8 as individuals.
# This broke at period 9 so testing if combining periods will
# allow it to work.
psi_hat_8 <- nr_out_inc_list[[9]]
beta_1_start <- c(1:8, rep(9,2))
prm$beta_1_track <- c(beta_1_start, rep(10, 26-length(beta_1_start)))
psi_start <- c(psi_hat_8[1:8], psi_hat_8[[8]], psi_hat_8[[9]])

nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                     psi_start_vec=psi_start,
                                     print_results=T)



psi_hat_8 <- nr_out_inc_list[[8]]
beta_1_start <- c(1:6, rep(7,2), 8:9)
prm$beta_1_track <- c(beta_1_start, rep(10, 26-length(beta_1_start)))
psi_start <- c(psi_hat_8[1:6], psi_hat_8[[7]], rep(psi_hat_8[[8]],3))

nr_out_9 <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                     psi_start_vec=psi_start,
                                     print_results=T)


psi_hat_9 <- nr_out_9[[1]]
#psi_hat_8 <- nr_out_inc_list[[8]]
beta_1_start <- c(1:6, rep(7,2), 8:10, rep(11,2), 12:14)
prm$beta_1_track <- c(beta_1_start, rep(15, 26-length(beta_1_start)))
psi_start <- c(psi_hat_9[1:9], rep(psi_hat_9[[10]],6))

nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                     psi_start_vec=psi_start,
                                     print_results=T)

# nr_out_13_16 <- lapply(1:4, 
#                        function(i){
#                          k_max <- 10+i
#                          beta_1_start <- c(1:k_max, rep(k_max+1, 5-i))
#                          prm$beta_1_track <- c(beta_1_start, rep(max(beta_1_start)+1, 26-length(beta_1_start)))
#                          psi_start <- c(psi_hat_old[1:10], rep(psi_hat_old[10], i), psi_hat_old[11:length(psi_hat_old)])
#                          nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
#                                                               psi_start_vec=psi_start,
#                                                               print_results=T)
#                          return(nr_out[[1]])
#                        })






#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           Building models with a "moving" average of 4 periods
#           First 4-4-x then 1-4-3-x, 2-4-2-x etc
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##

prm$beta_1_track <- rep(1,26)
nr_out_MA_list <- list()
df_MA <- tibble(period=0:8)
for(i in 0:4){
  if(i==0){
    prm$beta_1_track <- c(rep(1,4), rep(2,4), rep(3, 18))
  } else if(i < 4){
    prm$beta_1_track <- c(rep(1,i), rep(2, 4), rep(3, 4-i), rep(4, 18))
  } else {
    prm$beta_1_track <- c(rep(1,4), rep(2,4), rep(3, 18))
  }

  if(i == 0){
    psi_start <- rep(0, 3)
  } else if(i == 1){
    psi_start <- c(psi_hat_old[[1]], psi_hat_old[[2]], psi_hat_old[[2]],
                   psi_hat_old[[3]])
  } else if (i %in% c(2,3)) {
    psi_start <- psi_hat_old
  } else{
    psi_start <- c(psi_hat_old[[1]], psi_hat_old[[2]], psi_hat_old[[4]])
  }

  nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                       psi_start_vec=psi_start,
                                       print_results=T)
  psi_hat_old <- nr_out[[1]]
  VCV <- df %>% calculate_variance(prm, psi_hat_old, trt_models)
  
  ste_vec <- sapply(1:(dim(VCV)[[1]]), 
                     function(k){VCV[k,k] %>% sqrt()})
  
  psi_by_period <- sapply(1:9, 
                    function(k){
                      -psi_hat_old[[prm$beta_1_track[[k]]]]  
                    })
  ste_by_period <- sapply(1:9, 
                          function(k){
                            ste_vec[[prm$beta_1_track[[k]]]]  
                          })
  psi_col <- as.name(paste0("itr_", i))
  ste_col <- as.name(paste0("itr_", i, "_ste"))
  df_MA <- df_MA %>% add_column("{psi_col}" := psi_by_period,
                                "{ste_col}" := ste_by_period)
  
  #nr_out_MA_list <- nr_out_MA_list %>% append(., list(psi_hat_old))
}


ggplot(df_MA, aes(x=period)) +
  geom_point(aes(y=itr_0, size=4), color="red") +
  geom_line(aes(y=itr_0), color="red") +
  geom_errorbar(aes(ymin=itr_0-2*itr_0_ste, ymax=itr_0+2*itr_0_ste),
                width=.2, position=position_dodge(0.05), color="red") + 
  geom_point(aes(y=itr_1, size=4), color="darkred") +
  geom_line(aes(y=itr_1), color="darkred") +
  geom_errorbar(aes(ymin=itr_1-2*itr_1_ste, ymax=itr_1+2*itr_1_ste),
                width=.2, position=position_dodge(0.05), color="darkred") + 
  geom_point(aes(y=itr_2, size=4), color="blue") +
  geom_line(aes(y=itr_2), color="blue") +
  geom_errorbar(aes(ymin=itr_2-2*itr_2_ste, ymax=itr_2+2*itr_2_ste),
                width=.2, position=position_dodge(0.05), color="blue") + 
  geom_point(aes(y=itr_3, size=4), color="lightblue") +
  geom_line(aes(y=itr_3), color="lightblue") +
  geom_errorbar(aes(ymin=itr_3-2*itr_3_ste, ymax=itr_3+2*itr_3_ste),
                width=.2, position=position_dodge(0.05), color="lightblue") + 
  geom_point(aes(y=itr_4, size=4), color="green") + 
  geom_line(aes(y=itr_4), color="green") +
  geom_errorbar(aes(ymin=itr_4-2*itr_4_ste, ymax=itr_4+2*itr_4_ste),
                width=.2, position=position_dodge(0.05), color="green")






#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#         Lyle's block models (tm)
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##

#nr_out_MA_list <- list()
df_Lyle_psis <- tibble(psi=1:5)
for(i in 1:8){
  prm$beta_1_track <- c(rep(1, i), rep(2,4), rep(3,4), rep(4,4), rep(5, 26-i-3*4))
  
  if(i == 1){
    psi_start <- c(df_MA$itr_0[[1]], rep(df_MA$itr_0[[5]],3),
                   df_MA$itr_0[[9]])
  } else {
    psi_start <- psi_hat_old
  }
  
  #nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
  #                                     psi_start_vec=psi_start,
  #                                     print_results=T)
  
  
  nr_out_piecewise <- df %>% newton_raphson_piece(prm=prm, tol=0.01,
                                       psi_start_vec=psi_start,
                                       max_sub_iter = 20,
                                       print_results=T)
  psi_hat <- nr_out_piecewise[[1]]
  #VCV <- df %>% calculate_variance(prm, psi_hat_old, trt_models)
  
  #ste_vec <- sapply(1:(dim(VCV)[[1]]), 
  #                  function(k){VCV[k,k] %>% sqrt()})
  
  #psi_by_period <- sapply(1:9, 
  #                        function(k){
  #                          -psi_hat_old[[prm$beta_1_track[[k]]]]  
  #                        })
  #ste_by_period <- sapply(1:9, 
  #                        function(k){
  #                          ste_vec[[prm$beta_1_track[[k]]]]  
  #                        })
  psi_col <- as.name(paste0("itr_", i))
  #ste_col <- as.name(paste0("itr_", i, "_ste"))
  #df_MA <- df_MA %>% add_column("{psi_col}" := psi_by_period,
  #"{ste_col}" := ste_by_period)
  
  df_Lyle_psis <- df_Lyle_psis %>% add_column("{psi_col}" := psi_hat)
  
  psi_hat_old <- psi_hat
  #nr_out_MA_list <- nr_out_MA_list %>% append(., list(psi_hat_old))
}

df_lyle <- tibble(period=0:25)
matrix_lyle <- df_Lyle_psis %>% as.matrix()
for(i in 1:8){
  prm$beta_1_track <- c(rep(1, i), rep(2,4), rep(3,4), rep(4,4), rep(5, 26-i-3*4))
  psi_hat <- matrix_lyle[1:5, i+1]
  
  psi_by_period <- sapply(1:26, 
                    function(k){
                      return(-psi_hat[prm$beta_1_track[[k]]])
                    })
  psi_col <- as.name(paste0("itr_", i))
  df_lyle <- df_lyle %>% add_column("{psi_col}" := psi_by_period)
}

df_lyle_2 %>% pivot_longer(cols = !period) %>% 
  ggplot(aes(x=period, y=value, group=name)) +
  geom_point() + geom_line()


