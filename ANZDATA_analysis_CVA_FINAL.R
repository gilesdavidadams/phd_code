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
#df_all <- df_all %>% filter(usemarker == 1)
df_filtered <- df_all %>%
                filter(totaldialdur > 90,
                       !is.na(bmi),
                       bmi > 15,
                       bmi <= 50,
                       !is.na(sercreat),
                       !is.na(rxhomeva),
                       !is.na(rxhomeday90),
                       !is.na(vacategory90)) %>%
                      rename(t0 = "_t0", t="_t", 
                            smoke=cig, 
                            latereferral=referral3) %>%
                      filter(!latereferral=="") %>%
                      filter(smoke!=4)


#collect baseline variables
df_baseline <- df_filtered %>%
  mutate(cad=BEcoronary,
         lung=BElung,
         pvd=BEpvd,
         diabetes=diabB,
         at_home=rxhome,
         race=racereduced) %>%
  group_by(id) %>%
  select(id, firstdate, finaldate,
         age, sex, race, 
         sercreat, smoke,
         latereferral, cad, lung, pvd,
         diabetes, bmi, at_home,
         usemarker, censordate) %>%
  filter(row_number()==1) %>%
  ungroup()

df_base_summarise <- df_filtered %>%
  mutate(t_period=ifelse(t %% 90==0, t-1, t))
df_base_summarise <- df_base_summarise %>%
  group_by(id) %>% 
  summarise(lastperiod=max((t_period %/% 90)-1),
            censorind=max(censorind))
df_baseline <- df_baseline %>% 
  left_join(df_base_summarise, by="id")

df_baseline$bmi_cat <- cut(df_baseline$bmi, c(0, 20, 25, 30, 100))

df_hasPD <- df_filtered %>% 
              mutate(isPD = ifelse(t0 > 90 & rxhomeva==6, 1, 0)) %>%
              group_by(id) %>%
              summarise(hasPD = max(isPD))
              
df_filtered_PD <- df_filtered %>% 
                    left_join(df_hasPD, by="id") %>%
                    filter(hasPD==0) %>%
                    select(-c(hasPD))
              

#choose variables/columns that we want to keep
df_reduced <- df_filtered_PD %>% 
  select(c(id, sequno, firstdate, finaldate, deathdate, deadind,
           survtime, rxhomeva, t, t0,
           age, sex, racereduced, sercreat, smoke,
           latereferral, rxhome, censordate, censorind,
           lu, pv
  ))
  #select(-c(sequno, firstdate, finaldate, deathdate, survtime)) %>%
  



#split the dataset at the period cut-off points
df_split <- survSplit(Surv(t0, t, deadind)~., df_reduced, cut=seq(90,3650,90))
#create a variable indicating which period each row belongs to
df_split <- df_split %>% mutate(period= (t0 %/% 90)-1)

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

for(i in 0:max(df_split$rxhomeva, na.rm=TRUE)){
  if(i==0){
    progressbar <- txtProgressBar(min=0, max=max(df_split$rxhomeva, na.rm=TRUE),
                                  style=3, width=50, char="=")}
  else{
    setTxtProgressBar(progressbar, i) }
  
  rxhomeva_sum_cols <- c(rxhomeva_sum_cols, paste0("rxhomeva_sum_", i))
  col_name <- as.name(paste0("rxhomeva_", i))
  df_split_grp <- df_split %>% group_by(id, period) %>%
    summarise("rxhomeva_sum_{i}":= sum({{col_name}}))
  
  df_split <- df_split %>% 
    left_join(df_split_grp,
              by=c("id", "period")) %>%
    select(-c({{col_name}}))
  
  if(i==max(df_split$rxhomeva, na.rm=TRUE)){close(progressbar)}
}


df_split <- df_split %>% 
  mutate(rxhomeva_max = 
           colnames(df_split[, rxhomeva_sum_cols])[max.col(df_split[, rxhomeva_sum_cols])]) %>%
  select(-all_of(rxhomeva_sum_cols)) %>%
  mutate(rxhomeva_90=as.integer(substring(rxhomeva_max, 14)))



#Calculates ti for each individual
df_split <- df_split %>% 
  group_by(id) %>%
  mutate(ti=max(t)) %>%
  ungroup()

#Calculates which treatment the patient is on the longest
# within each period
df_90_original <- df_split %>%
  group_by(id, period) %>% 
  summarise(rxhomeva_90=max(rxhomeva_90),
            deadind=max(deadind),
            ti=max(ti),
            tk=sum(t-t0)) %>%
  ungroup()



# df_90_wide <- df_90_original %>% select(-c(tk)) %>%
#           pivot_wider(names_from = period,
#                       names_prefix = "rx_",
#                       values_from = rxhomeva_90)


df_90 <- df_90_original %>% 
  group_by(id) %>% 
  mutate(dies=max(deadind)) %>% 
  ungroup() %>%
  select(-c(deadind, tk)) %>%
  rename(ai = rxhomeva_90)



df_90_NA <- df_90 %>%
  mutate(ai_NA = as.numeric(is.na(ai))) %>%
  group_by(id) %>%
  summarise(has_NA = max(ai_NA))

df_90_home <- df_90 %>%
  mutate(ai_home = ifelse(ai %in% c(0,1,2), 1, 0)) %>%
  group_by(id) %>%
  summarise(has_home = max(ai_home))

# df_90_incPD <- df_90 %>%
#   mutate(is.PD = (ai == 6)) %>%
#   group_by(id) %>%
#   summarise(has_PD = max(as.numeric(is.PD)))

# df_baseline <- df_baseline %>%
#   left_join(df_90_NA, by="id") %>%
#   left_join(df_90_incPD, by="id")




df_90_filtered <- df_90 %>%
  left_join(df_90_NA, by="id") %>%
  filter(has_NA != 1) %>%
  left_join(df_90_home, by="id") %>%
  filter(has_home != 1) %>%
  rename(ai_finest=ai) %>%
  filter(ti > 90)



# df_90_PD <- df_90 
# 
# df_90_hasPD <- df_90 %>% group_by(id) %>%
#   summarise(has.PD = max(as.numeric(is.PD)))
# 
# df_90 <- df_90 %>%
#   left_join(df_90_hasPD, by="id") %>%
#   filter(ti>90)

#df_90_wide <- df_90 %>% select(-c(ai_NA, has_NA, is.PD, has.PD)) %>%
#              pivot_wider(names_from = period,
#                    names_prefix = "a_",
#                    values_from = ai_finest)

# df_90_home <- df_90 %>% mutate(HvF = ifelse(ai_finest %in% c(0,1,2), 1, 0)) %>%
#   select(c(id, period, HvF)) %>%
#   pivot_wider(names_from = period,
#               names_prefix = "HvF_",
#               values_from = HvF)

# df_CVC_vs_AVFG_long <- df_90 %>%
#                         filter(has.PD != 1) %>%
#                         mutate(ai = ifelse(ai_finest %in% c(2, 5), 1, 0))



df_CVC_vs_AVFG_long <- df_90_filtered %>%
  # filter(has.PD != 1) %>%
  #mutate(ai = ifelse(ai_finest %in% c(2, 5), 1, 0))
  mutate(ai = ifelse(ai_finest %in% c(2, 5), 0, 1))

df_CVC_vs_AVFG_wide <- df_CVC_vs_AVFG_long %>% 
  select(c(id, period, ti, dies, ai)) %>%
  pivot_wider(names_from = period,
              names_prefix = "a_",
              values_from = ai)



#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           C vs A
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##

df_cva <- df_CVC_vs_AVFG_wide
df_cva <- df_cva %>% left_join(df_baseline, by="id")
df_cva <- df_cva %>% mutate(sex=as.factor(sex),
                            race=factor(race, labels = c("NZ maori/pacific" = 1,
                                                         "Aus Indigenous" = 2,
                                                         "Asian" = 3,
                                                         "White/other" = 4)),
                            firstyear = as.factor(format(firstdate,"%Y")),
                            smoke = as.factor(smoke),
                            x=0
)
# df_cva <- df_cva %>% filter(!is.na(sercreat)) %>% 
#   filter(!latereferral=="") %>%
#   filter(smoke!=4) %>%
#   filter(!is.na(bmi_cat)) %>%
#   filter(usemarker==1)




prm_27 <- list()
prm_27$sim_label <- "(const_1)"
# prm_27$t_a_vec <- c(seq(90, 90*28, 90))
prm_27$psi_lab <- c("psi_1")
prm_27$censor_date <- df_all$finaldate %>% max()
prm_27$censor <- T


prm_27_full <- prm_27
prm_27_full$t_a_vec <- c(seq(90, 90*27, 90))
prm_27_full$beta_1_track <- c(1:16, rep(17, 27-16))
prm_27_full$beta_x_track <- rep(0, 27)

final_time <- prm_27_full$t_a_vec[length(prm_27_full$t_a_vec)]

df_cva_27 <- df_cva %>% #filter(firstdate < prm_27$censor_date) %>%
  filter(!(censordate < prm_27$censor_date) | (is.na(censordate))) %>%
  mutate(C_i = as.integer(difftime(prm_27$censor_date, firstdate,  "days"))) %>%
  mutate(C_i = pmax(C_i, final_time )) %>%
  mutate(ti = pmin(ti, C_i),
         iscensored = ifelse(ti==C_i, 1, 0))

prm_27_full$censor_max <- df_cva_27 %>% select(C_i) %>% max()




#Full treatment model

prm_27_full$trt_mod_list_full <- lapply(0:26, function(k) {
  mdl_k <- c("1", "age", "sex", "race", "smoke", "firstyear",
             "sercreat", "latereferral", "cad", "lung", "pvd",
             "diabetes", "bmi_cat"#, "at_home"
  )
  return(mdl_k)
})
#df_cva_27_full <- df_cva_27 %>% filter(firstdate < prm_27_full$censor_date) %>%
#          mutate(C_i = difftime(prm_27_full$censor_date, firstdate,  "days"))
fit_trt_out <- df_cva_27 %>% fit_treatment_models(prm=prm_27_full)
df_cva_27_full <- fit_trt_out[[1]]
prm_27_full$trt_models <- fit_trt_out[[2]]


nr_out_27 <- df_cva_27_full %>%
  newton_raphson_piece(prm=prm_27_full,
                       psi_max = 4,
                       max_iter=200)

psi_hat_vec_27 <- nr_out_27[[1]]
psi_hat_vec <- psi_hat_vec_27

VCV <- df_cva_27_full %>% calculate_variance(prm=prm_27_full,
                                             psi_hat_vec=psi_hat_vec)
ste <- sapply(1:dim(VCV)[[1]],
              function(k){
                return(sqrt(VCV[k,k]))
              })
var_vec <- ste^2

tibble(psi=psi_hat_vec, ste=ste, VCV=ste^2)

df_psi_27 <- lapply(1:(prm_27_full$beta_1_track %>% length()),
                    function(k){
                      beta_value <- prm_27_full$beta_1_track[k]
                      psi_now <- psi_hat_vec_27[beta_value]
                      ste_now <- ste[beta_value]
                      
                      tibble(period=k-1, psi=psi_now, ste=ste_now)
                    }) %>% Reduce(rbind, .)
df_psi_27 %>% ggplot(aes(x=period, y=psi)) + 
  geom_point() +
  geom_errorbar(aes(ymin=psi-2*ste, ymax=psi+2*ste),
                width=.2, position=position_dodge(0.05))

# INVERSE VARIANCE WEIGHTED AVERAGING

ivw_psi_hat <- sapply(1:26,
                      function(j){
                        if(j==1){
                          weights <- var_vec[1:2]
                          psi_hat_part <- psi_hat_vec[1:2]
                        } else if (j==16){
                          weights <- var_vec[15:16]
                          psi_hat_part <- psi_hat_vec[15:16]
                        } else if (j >= 17) {
                          psi_hat_part <- psi_hat_vec[16]
                          weights <- 1
                        } else {
                          weights <- var_vec[(j-1):(j+1)]
                          psi_hat_part <- psi_hat_vec[(j-1):(j+1)]
                        }
                        return(sum(psi_hat_part/weights)/(sum(1/weights)))
                      })

df_ivw <- lapply(1:(prm_27_full$beta_1_track %>% length()),
                 function(k){
                   psi_now <- ivw_psi_hat[prm_27_full$beta_1_track[k]]
                   
                   tibble(period=k-1, psi=psi_now)
                 }) %>% Reduce(rbind, .)

df_ivw %>% ggplot(aes(x=period, y=psi)) + geom_point()







df_ivw_5 <- lapply(1:26,
                   function(j){
                     if(j==1){
                       items <- 1:3
                       # weights <- var_vec[1:3]
                       # psi_hat_part <- psi_hat_vec[1:3]
                     } else if(j==2) {
                       items <- 1:4
                       # weights <- var_vec[1:4]
                       # psi_hat_part <- psi_hat_vec[1:4]
                     } else if (j==15){
                       items <- 13:16
                       # weights <- var_vec[13:16]
                       # psi_hat_part <- psi_hat_vec[13:16]
                     } else if (j==16){
                       items <- 14:16
                       # weights <- var_vec[14:16]
                       # psi_hat_part <- psi_hat_vec[14:16]
                     } else if (j >= 17) {
                       items <- 17
                       # psi_hat_part <- psi_hat_vec[16]
                       # weights <- 1
                     } else {
                       items <- (j-2):(j+2)
                       # weights <- var_vec[(j-2):(j+2)]
                       # psi_hat_part <- psi_hat_vec[(j-2):(j+2)]
                     }
                     weights <- var_vec[items]
                     psi_hat_part <- psi_hat_vec[items]
                     #var_parts <- var_vec[items]
                     
                     psi_avg <- sum(psi_hat_part/weights)/(sum(1/weights))
                     var_avg <- sum(outer(1/sqrt(weights), 1/sqrt(weights)))/((sum(1/weights))^2)
                     # (sum(1/weights))^(-1)*(sum(outer(1/sqrt(weights, 1/weights)))
                     return(tibble(period=j, psi=-psi_avg, var=var_avg, ste=sqrt(var_avg)))
                   }) %>% Reduce(rbind, .)

# df_ivw_5 <- lapply(1:(prm_27_full$beta_1_track %>% length()),
#                  function(k){
#                    psi_now <- ivw_psi_hat_5[prm_27_full$beta_1_track[k]]
#                    
#                    tibble(period=k-1, psi=-psi_now)
#                  }) %>% Reduce(rbind, .)

df_ivw_5 %>% ggplot(aes(x=period, y=psi)) + 
  geom_point() + 
  geom_errorbar(aes(ymin=psi-2*ste, ymax=psi+2*ste),
                width=.2, position=position_dodge(0.05))






n=15
df_two_stage_list <- tibble()
for(n in 1:15){
  #for(n in 15:20){
  prm_temp <- prm_20_full
  prm_temp$beta_1_track <- c(rep(1, n), rep(2, 20-n))
  
  nr_out_20 <- df_cva_20_full %>%
    newton_raphson_piece(prm=prm_20_full,
                         psi_max = 4,
                         max_iter=200)
  
  psi_hat_vec <- nr_out[[1]]
  
  VCV <- df_cva_20_full %>% calculate_variance(prm=prm_20_full,
                                               psi_hat_vec=psi_hat_vec)
  #trt_models=trt_models_25_full)
  
  ste <- sapply(1:dim(VCV)[[1]],
                function(k){
                  return(sqrt(VCV[k,k]))
                })
  
  df_two_stage_list <- list(df_two_stage_list,
                            tibble(psi_1 = psi_hat_vec[1],
                                   psi_2 = psi_hat_vec[2],
                                   ste_1 = ste[1],
                                   ste_2 = ste[2]))
  
}






fit_trt_out <- df_CvA_aug %>% fit_treatment_models(prm=prm)
df <- fit_trt_out[[1]]
trt_models <- fit_trt_out[[2]]


psi_extender_list <- list()

lapply(0:15, function(i){
  if(i > 0){
    prm_25_full$beta_1_track <- c(1:i, rep(i+1, 26-i))
  } else {
    prm_25_full$beta_1_track <- rep(1,26)
  }
  
  nr_out <- df_cva_25_full %>% newton_raphson_piece(prm=prm_25_full)
  psi_hat_vec <- nr_out[[1]]
  
  tibble(psi_hat_vec)
}) %>% Reduce(rbind, .)




psi_indiv_list <- list()
ste_indiv_list <- list()
for(i in 0:8){
  #for(i in 7:10){
  #psi_indiv_list <- lapply(1:10,
  #function(i){
  if(i > 0){
    prm_25_full$beta_1_track <- c(1:i, rep(i+1, 26-i))
  } else {
    prm_25_full$beta_1_track <- rep(1,26)
  }
  
  # if(i == 0){
  #   psi_start_vec <- c(0)
  # } else if(i == 1){
  #   psi_start_vec <- psi_indiv_list[[i]]
  #   psi_start_vec <- c(psi_start_vec, psi_start_vec)
  # } else {
  #   psi_start_vec <- psi_indiv_list[[i]]
  #   psi_start_vec <- c(psi_start_vec[1:(i-1)],
  #                    psi_start_vec[i],
  #                    psi_start_vec[i:length(psi_start_vec)])
  # }
  
  nr_out <- df_cva_25_full %>% newton_raphson_piece(prm=prm_25_full, #tol=0.01,
                                                    #psi_start_vec=psi_start_vec,
                                                    psi_max_vec = c(rep(4,7), rep(3,4)),
                                                    max_iter=200)
  psi_hat_vec <- nr_out[[1]]
  # return(tibble(psi_hat=psi_hat_vec))
  #})
  psi_indiv_list <- psi_indiv_list %>% append(., list(psi_hat_vec))
  
  # VCV <- df_cva_25_full %>% calculate_variance(prm=prm_25_full,
  #                                      psi_hat_vec=psi_hat_vec)
  #                                      #trt_models=trt_models_25_full)
  
  # ste <- sapply(1:ifelse(i==0, 1, dim(VCV)[[1]]),
  #               function(k){
  #                 return(ifelse(i==0, sqrt(VCV), sqrt(VCV[k,k])))
  #               })
  # ste_indiv_list <- ste_indiv_list %>% append(., list(ste))
}



df_indiv <- tibble(period=0:25)
ste_indiv <- tibble(period=0:25)
for(i in 1:length(psi_indiv_list)){
  psi_hat <- psi_indiv_list[[i]]
  ste_hat <- ste_indiv_list[[i]]
  
  psi_by_period <- sapply(1:26, 
                          function(k){
                            if(k < length(psi_hat)){
                              return(-psi_hat[k])
                            } else {
                              return(-psi_hat[length(psi_hat)])
                            }
                          })
  
  ste_by_period <- sapply(1:26, 
                          function(k){
                            if(k < length(ste_hat)){
                              return(ste_hat[k])
                            } else {
                              return(ste_hat[length(ste_hat)])
                            }
                          })
  
  psi_col <- as.name(paste0("itr_", i-1))
  #ste_col <- as.name(paste0("itr_", i, "_ste"))
  df_indiv <- df_indiv %>% add_column("{psi_col}" := psi_by_period)
  ste_indiv <- ste_indiv %>% add_column("{psi_col}" := ste_by_period)
  #nr_out_MA_list <- nr_out_MA_list %>% append(., list(psi_hat_old))
}


df_two_stage_list <- tibble()
for(n in 1:15){
  #for(n in 15:20){
  prm_temp <- prm_25_full
  prm_temp$beta_1_track <- c(rep(1, n), rep(2, 25-n))
  
  nr_out <- df_cva_25_full %>%
    newton_raphson_piece(prm=prm_temp,
                         psi_max = 4,
                         max_iter=200)
  
  psi_hat_vec <- nr_out[[1]]
  
  VCV <- df_cva_25_full %>% calculate_variance(prm=prm_temp,
                                               psi_hat_vec=psi_hat_vec)
  #trt_models=trt_models_25_full)
  
  ste <- sapply(1:dim(VCV)[[1]],
                function(k){
                  return(sqrt(VCV[k,k]))
                })
  
  df_two_stage_list <- list(df_two_stage_list,
                            tibble(psi_1 = psi_hat_vec[1],
                                   psi_2 = psi_hat_vec[2],
                                   ste_1 = ste[1],
                                   ste_2 = ste[2]))
  
}




ste_indiv_list <- list()
for(i in 0:8){
  #for(i in 7:10){
  #psi_indiv_list <- lapply(1:10,
  #function(i){
  if(i > 0){
    prm_25_full$beta_1_track <- c(1:i, rep(i+1, 26-i))
  } else {
    prm_25_full$beta_1_track <- rep(1,26)
  }
  
  nr_out <- df_cva_25_full %>% newton_raphson_piece(prm=prm_25_full, #tol=0.01,
                                                    #psi_start_vec=psi_start_vec,
                                                    psi_max_vec = c(rep(4,7), rep(3,4)),
                                                    max_iter=200)
  psi_hat_vec <- nr_out[[1]]
  # return(tibble(psi_hat=psi_hat_vec))
  #})
  psi_indiv_list <- psi_indiv_list %>% append(., list(psi_hat_vec))
  
}





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
  VCV <- df %>% calculate_variance(prm, psi_hat_old)
  
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
df_Lyle_stes <- tibble(ste=1:5)
for(i in 0:8){
  
  if(i == 0){
    prm_25_full$beta_1_track <- c(rep(1,4), rep(2,4), rep(3,4), rep(4, 26-i-3*4))
  } else {
    prm_25_full$beta_1_track <- c(rep(1, i), rep(2,4), rep(3,4), rep(4,4), rep(5, 26-i-3*4))
  }
  
  if(i == 0){ psi_start <- NA}
  else {
    psi_start <- psi_hat_old
  }
  
  #nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
  #                                     psi_start_vec=psi_start,
  #                                     print_results=T)
  
  
  nr_out_piecewise <- df_cva_25_full %>% 
    newton_raphson_piece(prm=prm_25_full)#, tol=0.01,
  #psi_start_vec=psi_start,
  #max_sub_iter = 20,
  #print_results=T)
  psi_hat <- nr_out_piecewise[[1]]
  
  VCV <- df %>% calculate_variance(prm, psi_hat)#, trt_models)
  ste_vec <- sapply(1:(dim(VCV)[[1]]), 
                    function(k){VCV[k,k] %>% sqrt()})
  
  
  if(i == 0){
    psi_hat <- c(psi_hat, NA)
    #ste_vec <- c(ste_vec, NA)
  }  
  
  psi_col <- as.name(paste0("itr_", i))
  
  df_Lyle_psis <- df_Lyle_psis %>% add_column("{psi_col}" := psi_hat)
  # df_Lyle_stes <- df_Lyle_stes %>% add_column("{psi_col}" := ste_vec)
  
  psi_hat_old <- psi_hat
}



df_lyle <- tibble(period=0:25)
df_lye_stes <- tibble(period=0:25)
matrix_lyle_stes <- df_Lyle_stes %>% as.matrix()
matrix_lyle <- df_Lyle_psis %>% as.matrix()
for(i in 0:8){
  if(i == 0){
    prm$beta_1_track <- c(rep(1,4), rep(2,4), rep(3,4), rep(4, 26-i-3*4))
  } else {
    prm$beta_1_track <- c(rep(1, i), rep(2,4), rep(3,4), rep(4,4), rep(5, 26-i-3*4))
  }
  #prm$beta_1_track <- c(rep(1, i), rep(2,4), rep(3,4), rep(4,4), rep(5, 26-i-3*4))
  psi_hat <- matrix_lyle[1:5, i+2]
  ste_hat <- matrix_lyle_stes[1:5, i+2]
  
  psi_by_period <- sapply(1:26, 
                          function(k){
                            return(-psi_hat[prm$beta_1_track[[k]]])
                          })
  ste_by_period <- sapply(1:26, 
                          function(k){
                            return(ste_hat[prm$beta_1_track[[k]]])
                          })
  psi_col <- as.name(paste0("itr_", i))
  df_lyle <- df_lyle %>% add_column("{psi_col}" := psi_by_period)
  df_lye_stes <- df_lye_stes %>% add_column("{psi_col}" := ste_by_period)
}

matrix_lyle_by_period <- df_lyle %>% select(-c(period)) %>% as.matrix()

MA_design <- c(rep(1,1), rep(0,8), #1
               rep(1,2), rep(0,7),
               rep(1,3), rep(0,6),
               rep(1,4), rep(0,5),
               rep(1/2,1), rep(1,3), rep(1/2,1), rep(0,4), #5
               rep(1/2,2), rep(1,2), rep(1/2,2), rep(0,3),
               rep(1/2,3), rep(1,1), rep(1/2,3), rep(0,2),
               rep(1/2,4), rep(0,0), rep(1/2,4), rep(0,1),
               rep(1/3,1), rep(1/2,3), rep(1/3), rep(1/2, 3), rep(1/3,1),
               rep(1/3,1), rep(1/2,3), rep(1/3), rep(1/2, 3), rep(1/3,1), #10
               rep(1/3,1), rep(1/2,3), rep(1/3), rep(1/2, 3), rep(1/3,1),
               rep(1/3,1), rep(1/2,3), rep(1/3), rep(1/2, 3), rep(1/3,1),
               rep(0,1), rep(1/2,8),
               rep(0,2), rep(1/2,3), rep(1,1), rep(1/2,3),
               rep(0,3), rep(1/2,2), rep(1,2), rep(1/2,2), #15
               rep(0,4), rep(1/2,1), rep(1,3), rep(1/2,1),
               rep(0,5), rep(1,4),
               rep(0,6), rep(1,3),
               rep(0,7), rep(1,2), 
               rep(0,8), rep(1,1)) %>% matrix(ncol=9, byrow=T) %>% t()

MA_converter <- sapply(1:20, function(k){MA_design[,k] / (sum(MA_design[,k]))})

psi_MA_smoothed <- sapply(1:dim(MA_converter)[[2]],
                          function(k){
                            return(matrix_lyle_by_period[k,]%*%MA_converter[,k])
                          })

df_lyle_2 %>% pivot_longer(cols = !period) %>% 
  ggplot(aes(x=period, y=value, group=name)) +
  geom_point() + geom_line()


