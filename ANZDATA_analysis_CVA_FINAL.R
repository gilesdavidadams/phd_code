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
library(gtsummary)
library(haven)
library(survival)
library(stats)
library(labelled)
library(janitor)
library(labelled)
library(DescTools)
library(cowplot)
library(xtable)
library(MASS, exclude=c("select"))
source("G-estimation_FINAL.R")


options(dplyr.summarise.inform = FALSE)




#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           DATA CLEANING
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##


#read in data
df_all <- read_dta("C:/Users/gdadams/OneDrive - The University of Melbourne/ANZREQ-194 Stata datasets/mergedJK_20131018_FINAL.dta")
df_filtered <- df_all %>%
                filter(totaldialdur > 90,
                       !is.na(bmi),
                       bmi > 15,
                       bmi <= 50,
                       !is.na(sercreat),

                       !is.na(vacategory90)) %>%
                      rename(t0 = "_t0", t="_t"
                      ) %>%
  select(-c(race))

#collect baseline variables
df_baseline <- df_filtered %>%
  rename(cad=coronB,
         lung=luB,
         pvd=pvB,
         diabetes=diabB,
         at_home=rxhome,
         race=racereduced,
         latereferral=referral3,
         smoke=cig) %>%
  group_by(id) %>%
  filter(row_number()==1) %>%
  ungroup() %>% 
  mutate(sex=as.factor(sex),
         race=factor(race, levels=c(1:4),
                     labels = c("NZ maori/pacific",
                                "Aus Indigenous",
                                "Asian",
                                "White/other")),
         cad=factor(cad, levels=c(0,1), labels=c("No/Unknown", "Yes/Suspected")),
         lung=factor(lung, levels=c(0,1), labels=c("No/Unknown", "Yes/Suspected")),
         pvd=factor(pvd, levels=c(0,1), labels=c("No/Unknown", "Yes/Suspected")),
         diabetes=factor(diabetes, levels=c(0,1, 2), labels=c("No/Unknown", "Type 1", "Type 2")),
         firstyear = as.factor(format(firstdate,"%Y")),
         smoke=ifelse(smoke %in% c(3,4), 0, smoke),
         smoke=factor(smoke, levels=c(0,1,2), labels=c("Never/Unknown", "Current", "Former")),
         latereferral=ifelse(latereferral %in% c("U", "", "N"), 0, 1),
         latereferral=factor(latereferral, levels=c(0,1), labels=c("No/Unknown", "Yes")),
         primary=ifelse((disease %>% as.integer) %in% c(800, 801, 802, 803), 1,
                        ifelse((disease %>% as.integer) %in% c(100, 110:112, 121, 122,
                                                               seq(130, 190, 10), 182,
                                                               191), 2,
                        ifelse((disease %>% as.integer) %in% c(300, 302), 3, 4))),
         primary=factor(primary, levels=c(1,2,3,4), labels=c("Diabetes", "Glomerulonephritis",
                                                             "Hypertension", "Other"))
         ) %>%
  select(id, age, sex, race, 
         sercreat, smoke,
         latereferral, cad, lung, pvd,
         diabetes, bmi, censordate, firstdate, firstyear, primary)

df_baseline$bmi_cat <- cut(df_baseline$bmi, c(15, 20, 25, 30, 50))

df_base_summarise <- df_filtered %>%
  mutate(t_period=ifelse(t %% 90==0, t-1, t)) %>%
  group_by(id) %>% 
  summarise(lastperiod=max((t_period %/% 90)-1),
            censorind=max(censorind))

df_baseline <- df_baseline %>% 
  left_join(df_base_summarise, by="id")



df_modalityNA <- df_filtered %>% mutate(isNA = is.na(rxhomeva)) %>%
                    group_by(id) %>%
                    summarise(hasNA = max(isNA))

df_filtered_NA <- df_filtered %>% 
                    left_join(df_modalityNA, by="id") %>%
                    filter(hasNA==0) %>%
                    select(-c(hasNA))


df_hasPD <- df_filtered_NA %>% 
  mutate(isPD = ifelse(t > 90 & rxhomeva==6, 1, 0)) %>%
  group_by(id) %>%
  summarise(hasPD = max(isPD))
              
df_filtered_PD <- df_filtered_NA %>% 
                    left_join(df_hasPD, by="id") %>%
                    filter(hasPD==0)

df_hasHome <- df_filtered_PD %>% 
  mutate(isPD = ifelse(t > 90 & rxhomeva==6, 1, 0)) %>%
  group_by(id) %>%
  summarise(hasPD = max(isPD))
              

#choose variables/columns that we want to keep
df_reduced <- df_filtered_PD %>% 
  select(c(id, sequno, firstdate, finaldate, deathdate, deadind,
           survtime, rxhomeva, t, t0, censordate
  ))
  

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
  rxhomeva_sum_cols <- c(rxhomeva_sum_cols, paste0("rxhomeva_sum_", i))
  col_name <- as.name(paste0("rxhomeva_", i))
  df_split_grp <- df_split %>% group_by(id, period) %>%
    summarise("rxhomeva_sum_{i}":= sum({{col_name}}))
  
  df_split <- df_split %>% 
    left_join(df_split_grp,
              by=c("id", "period")) %>%
    select(-c({{col_name}}))
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



df_90_filtered <- df_90 %>%
  left_join(df_90_NA, by="id") %>%
  filter(has_NA == 0) %>%
  left_join(df_90_home, by="id") %>%
  filter(has_home == 0) %>%
  rename(ai_finest=ai) %>%
  filter(ti > 90)


df_CVC_vs_AVFG_long <- df_90_filtered %>%
  mutate(ai = ifelse(ai_finest %in% c(2, 5), 1L, 0L))
 
df_cva <- df_CVC_vs_AVFG_long %>% 
  select(c(id, period, ti, dies, ai)) %>%
  pivot_wider(names_from = period,
              names_prefix = "a_",
              values_from = ai) %>%
  left_join(df_baseline, by="id") %>%
  mutate(x=0)


prm_base <- list()
prm_base$sim_label <- "(const_1)"
prm_base$psi_lab <- c("psi_1")
prm_base$censor_date <- df_all$finaldate %>% max()
prm_base$censor <- T


#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#          EXPLORATORY DATA ANALYSIS
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##

#Count of time spent on AVF/AVG/CVC H or F

df_reduced %>%
  select(id, rxhomeva, t, t0) %>%
  filter(t > 90) %>%
  mutate(totaltime = t-t0) %>%
  group_by(rxhomeva) %>%
  summarise(time_count=sum(totaltime)) #%>%
  # select(time_count) %>% sum


#Count of time spent on AVF/AVG/CVC facility only after filtering
df_split %>%
    filter(id %in% df_90_filtered$id) %>% 
    select(id, rxhomeva, t, t0) %>%
    filter(t > 90) %>%
    mutate(totaltime = t-t0) %>%
    group_by(rxhomeva) %>%
    summarise(time_count=sum(totaltime))



# df_reduced  %>% 
  left_join(df_reduced %>% 
              mutate(homeorpd = ifelse(rxhomeva %in% c(0,1,2,6),1,0)) %>% 
              group_by(id) %>% 
              summarise(has_homepd = max(homeorpd)),
            by="id") %>% 
  filter(has_homepd==0) %>% 
  select(id, rxhomeva, t, t0) %>% 
  mutate(avfg = ifelse(rxhomeva %in% c(3,4), 1, 0), totaltime = t-t0) %>% 
  group_by(avfg) %>% 
  summarise(sum(totaltime))





df_count <- df_cva
period_stats <- lapply(-1:32, function(k){
  a_k <- as.name(paste0("a_", k))
  df_now <- df_cva %>% 
              filter(ti >= 90*(k+1))
  n_now <- df_now %>% nrow()
  a_1_now <- df_now %>% select({a_k}) %>% sum(na.rm=T)
  a_0_now <- n_now - a_1_now
  
  if(k==-1){
    switchTo0 <- NA
    switchTo1 <- NA
    switch_tot <- NA
    dies_tot <- NA
  } else {
    a_prev <- as.name(paste0("a_", k-1))
    df_now <- df_now %>% mutate(switchTo1 = ifelse({{a_prev}}==0 & {{a_k}}==1, 1, 0),
                                switchTo0 = ifelse({{a_prev}}==1 & {{a_k}}==0, 1, 0),
                                dies_now = ifelse((ti < 90*(k+2)) & dies==1, 1, 0))
    switchTo0 <- df_now %>% select(switchTo0) %>% sum(na.rm=T) %>% as.integer()
    switchTo1 <- df_now %>% select(switchTo1) %>% sum(na.rm=T) %>% as.integer()
    switch_tot <- switchTo0 + switchTo1 %>% as.integer()
    dies_tot <- df_now %>% select(dies_now) %>% sum(na.rm=T) %>% as.integer()
  }
  return(tibble(k=k, n=n_now, a_0=a_0_now, a_1=a_1_now, deaths=dies_tot, 
                switch_tot=switch_tot, switchTo0=switchTo0, switchTo1=switchTo1))
}) %>% Reduce(rbind, .) %>% 
        mutate(d = format(round(deaths/n*100, 1)), s = format(round(switch_tot/n*100, 1)))
period_stats %>% xtable(digts=0) %>% print(include.rownames=F, 
                         format.args = list(big.mark = ","))

period_stats %>% select(k, n, deaths) %>% filter(k >= 0) %>%
  ggplot(aes(x=k, y=(deaths/n)*100)) +
  geom_point() + geom_smooth(method="lm", formula = (y ~ 1 + x + x^2), se=F) +
  xlab("Period") + ylab("% Deaths") + theme_bw()



 


tk_max <- 25
pk_max <- 16
prm_ai <- prm_base
prm_ai$t_a_vec <- c(seq(90, 90*(tk_max+1), 90))
prm_ai$beta_1_track <- c(1:pk_max, rep(pk_max+1, (tk_max+1)-pk_max))
prm_ai$beta_x_track <- rep(0, tk_max+1)

final_time <- prm_ai$t_a_vec[length(prm_ai$t_a_vec)] + 90

df_cva_ai <- df_cva %>% #filter(firstdate < prm_26$censor_date) %>%
  filter(!(censordate < prm_ai$censor_date) | (is.na(censordate))) %>%
  mutate(C_i = as.integer(difftime(prm_ai$censor_date, firstdate,  "days"))) %>%
  mutate(C_i = pmin(C_i, final_time )) %>%
  mutate(ti = pmin(ti, C_i),
         iscensored = ifelse(ti==C_i, 1, 0))

prm_ai$censor_max <- df_cva_ai %>% select(C_i) %>% max()

prm_ai$trt_mod_list_full <- lapply(0:tk_max, function(k) {
    mdl_k <- c("1", "age", "sex", "race", "smoke", #"firstyear",
             "sercreat", "latereferral", "cad", "lung", "pvd", 
             "diabetes",  "primary", "bmi_cat"
    )
  return(mdl_k)
})

fit_trt_out <- df_cva_ai %>% fit_treatment_models(prm=prm_ai)
df_cva_ai_trt <- fit_trt_out[[1]]
prm_ai$trt_models <- fit_trt_out[[2]]

#for(j in 1:length(prm_ai$trt_models)){
# z_df <-   
tibs_full <- lapply(1:(tk_max+1),
      function(j){
          mod_values <- prm_ai$trt_models[[j]] %>% 
                              summary %>%
                              .$coefficients
          
          cov_df <- tibble(Covariate=mod_values %>% rownames,
                           Estimate=mod_values[,1],
                           'Std. Error'=mod_values[,2]
          )
          # %>%
          #                     as_tibble(rownames="Covariate") %>%
          #                     select(Covariate, Estimate, 'Std. Error')
          for(i in 1:(dim(cov_df)[[1]])){
            if(i==1){
              tibs <- cov_df[i,] %>% 
                      select(Covariate, Estimate) %>%
                      add_row(Covariate=cov_df[i,]$Covariate,
                              Estimate=cov_df[i,]$'Std. Error')
            } else {
              tibs <- tibs %>%
                      add_row(Covariate=cov_df[i,]$Covariate,
                              Estimate=cov_df[i,]$Estimate) %>%
                      add_row(Covariate=cov_df[i,]$Covariate,
                              Estimate=cov_df[i,]$'Std. Error')      
            }
          }
          tibs <- tibs %>% mutate(
                          est_str=Estimate %>% 
                                  #format(scientific=F) %>%
                                  round(digits=4) %>%
                                  format(scientific=F),
                          est_str=ifelse(row_number() %% 2 == 1,
                                         est_str,
                                         paste0("(", gsub(" ", "", est_str), ")")),
                          Covariate=ifelse(row_number() %% 2 ==1,
                                          Covariate,
                                          "")
                                  )
          
          if(j > 1) {
            tibs <- tibs %>% select(-c(Covariate))
          }
          new_name <- as.name(paste0("p_", j-1))
          tibs <- tibs %>% rename("{new_name}" := est_str)
          tibs <- tibs %>% select(-c(Estimate))
          return(tibs)
      }) %>% Reduce(cbind, .)

tibs_full %>% xtable() %>% print(include.rownames=F)



covariates_table <- lapply(1:(length(prm_ai$trt_models)-1),
      function(j){
        new_name <- as.name(paste0("p_", j))
        zvalues <- prm_ai$trt_models[[j]] %>% 
                    summary %>%
                    .$coefficients %>%
                    as_tibble(rownames="Covariate") %>%
                    select(Covariate, Estimate, 'Std. Error') %>%
                    rename("{new_name}" := 'z value')
        zvalues <- zvalues %>% mutate(
                    "{new_name}" := format(round({{new_name}}, digits=1))
        )
      if(j==1){
        return(zvalues)
      } else {
        return(zvalues %>% select(-c(Covariate)))
      }

}) %>% Reduce(cbind, .) 

covariates_table %>% xtable %>% print(include.rownames=F, 
                         format.args = list(big.mark = ","))





#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           C vs A  26
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
pk_max <- 23
prm_26_full <- prm_base
prm_26_full$t_a_vec <- c(seq(90, 90*26, 90))
prm_26_full$beta_1_track <- c(1:pk_max, rep(pk_max+1, 26-pk_max))
prm_26_full$beta_x_track <- rep(0, 26)

final_time <- prm_26_full$t_a_vec[length(prm_26_full$t_a_vec)] + 90

df_cva_26 <- df_cva %>%
  filter(!(censordate < prm_26_full$censor_date) | (is.na(censordate))) %>%
  mutate(C_i = as.integer(difftime(prm_26_full$censor_date, firstdate,  "days"))) %>%
  mutate(C_i = pmin(C_i, final_time )) %>%
  mutate(ti = pmin(ti, C_i),
         iscensored = ifelse(ti==C_i, 1, 0))

prm_26_full$censor_max <- df_cva_26 %>% select(C_i) %>% max()




#Full treatment model
prm_26_full$trt_mod_list_full <- lapply(0:(26-1), function(k) {
  mdl_k <- c("1", "age", "sex", "race", "smoke",
             "sercreat", "latereferral", "cad", "lung", "pvd",
             "diabetes", "bmi_cat", "primary"
  )
  return(mdl_k)
})


fit_trt_out <- df_cva_26 %>% fit_treatment_models(prm=prm_26_full)
df_cva_26_full <- fit_trt_out[[1]]
prm_26_full$trt_models <- fit_trt_out[[2]]


nr_out_26 <- df_cva_26_full %>%
  newton_raphson_grad_large(prm=prm_26_full,
                       psi_max = 4,
                       max_iter=200)

psi_hat_vec_26 <- nr_out_26[[1]]
psi_hat_vec <- psi_hat_vec_26

VCV <- df_cva_26_full %>% calculate_variance(prm=prm_26_full,
                                             psi_hat_vec=psi_hat_vec)



ste <- sapply(1:dim(VCV)[[1]],
              function(k){
                return(sqrt(VCV[k,k]))
              })

var_vec <- ste^2

tibble(psi=psi_hat_vec, ste=ste, VCV=ste^2)
tibble_cva <- tibble(psi=psi_hat_vec, ste=ste, VCV=ste^2)

df_psi_26 <- lapply(1:(prm_26_full$beta_1_track %>% length()),
                    function(k){
                      beta_value <- prm_26_full$beta_1_track[k]
                      psi_now <- psi_hat_vec_26[beta_value]
                      ste_now <- ste[beta_value]
                      
                      tibble(period=((k-1) %>% as.integer()), psi=psi_now, ste=ste_now, epsi=exp(psi))
                    }) %>% Reduce(rbind, .)

# INVERSE VARIANCE WEIGHTED AVERAGING

ivw_3 <- lapply(1:26,
                   function(j){
                     if(j==1){
                       items <- c(1,2)
                     } else if (j==pk_max+1) {
                       items <- c(pk_max, pk_max+1)
                     } else if (j > pk_max+1){
                       items <- c(pk_max+1)
                     } else {
                       items <- (j-1):(j+1)
                     }
                     weights <- var_vec[items]
                     psi_hat_part <- psi_hat_vec[items]
                     
                     psi_avg <- sum(psi_hat_part/weights)/(sum(1/weights))
                     var_avg <- mean(weights)

                     return(tibble(period=j, psi_avg=psi_avg, var_avg=var_avg, ste=sqrt(var_avg)))
                   })

df_ivw_3 <- lapply(1:(prm_26_full$beta_1_track %>% length()),
                 function(k){
                   psi_now <- ivw_3[[k]]$psi_avg
                   ste_now <- ivw_3[[k]]$var_avg %>% sqrt()
                   tibble(period=k-1, psi_ivw3=psi_now, ste_ivw3=ste_now)
                 }) %>% Reduce(rbind, .)

ivw_5 <- lapply(1:26,
                   function(j){
                     if(j==1){
                       items <- 1:3
                     } else if(j==2) {
                       items <- 1:4
                     } else if (j==pk_max) {
                       items <- c((pk_max-2):(pk_max), rep(pk_max+1, 2))
                     } else if (j==pk_max+1){
                       items <- c(pk_max-1, pk_max, rep(pk_max+1, 3))
                     } else if (j==pk_max+2){
                       items <- c(pk_max, rep(pk_max+1, 4))
                     } else if (j > pk_max+2){
                       items <- c(pk_max+1)
                     } else {
                       items <- (j-2):(j+2)
                     }
                     weights <- var_vec[items]
                     psi_hat_part <- psi_hat_vec[items]
                     
                     psi_avg <- sum(psi_hat_part/weights)/(sum(1/weights))
                     var_avg <- mean(weights)

                     return(tibble(period=j, psi_avg=psi_avg, var_avg=var_avg, ste=sqrt(var_avg)))
                   })

df_ivw_5 <- lapply(1:(prm_26_full$beta_1_track %>% length()),
                 function(k){
                   psi_now <- ivw_5[[k]]$psi_avg
                   var_now <- ivw_5[[k]]$var_avg
                   ste_now <- var_now %>% sqrt()
                   tibble(period=k-1, psi_ivw5=psi_now, ste_ivw5=ste_now)
                 }) %>% Reduce(rbind, .)




df_psi_all <- df_psi_26 %>% 
  left_join(df_ivw_3, by="period") %>%
  left_join(df_ivw_5, by="period") %>%
  mutate(period=as.integer(period))

df_psi_all %>%
  xtable(digts=0) %>% print(include.rownames=F, 
                         format.args = list(big.mark = ","))




y_lim_low <- -4
y_lim_high <- 4
p_plot_max <- 26

df_psi_all %>% filter(period <= p_plot_max) %>%
  ggplot(aes(x=period, y=psi)) + 
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=psi-2*ste, ymax=psi+2*ste),
                width=.3, linewidth=0.9,
                position=position_dodge(0.05)) + 
  ylim(y_lim_low, y_lim_high) + 
  ylab("psi") + theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", color="blue")

df_psi_all %>% filter(period <= p_plot_max) %>%
        ggplot(aes(x=period, y=psi_ivw3)) +
        geom_point(size = 2) +
        geom_errorbar(aes(ymin=psi_ivw3-2*ste_ivw3, ymax=psi_ivw3+2*ste_ivw3),
                      width=.3, linewidth=0.9,
                      position=position_dodge(0.05)) +
        ylim(y_lim_low, y_lim_high) + theme_bw() +
        geom_hline(yintercept=0, linetype="dashed", color="blue") +
        ylab("IVW3")


df_psi_all %>% filter(period <= p_plot_max) %>%
  ggplot(aes(x=period, y=psi_ivw5)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=psi_ivw5-2*ste_ivw5, ymax=psi_ivw5+2*ste_ivw5),
                width=.3, linewidth=0.9,
                position=position_dodge(0.05)) +
  ylim(y_lim_low, y_lim_high) + theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", color="blue") +
  ylab("IVW5")