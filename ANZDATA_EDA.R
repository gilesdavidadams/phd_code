# Counting participants in slices

df_dialysis <- read_dta("C:/Users/gdadams/OneDrive - The University of Melbourne/ANZREQ-194 Stata datasets/ANZREQ-194 dialysis.dta")

df_all %>% distinct(id, .keep_all=T) %>% nrow()

df_all %>% distinct(id, .keep_all=T) %>%
  filter(totaldialdur > 90) %>% nrow()

df_all %>% distinct(id, .keep_all=T) %>%
  filter(totaldialdur > 90,
         !is.na(bmi),
         bmi > 15,
         bmi <= 50) %>% nrow()

df_all %>% distinct(id, .keep_all=T) %>%
  filter(totaldialdur > 90,
         !is.na(bmi),
         bmi > 15,
         bmi <= 50,
         !is.na(sercreat)) %>% nrow()

# df_all %>% distinct(id, .keep_all=T) %>%
#   filter(totaldialdur > 90,
#          !is.na(bmi),
#          bmi > 15,
#          bmi <= 50,
#          !is.na(sercreat),
#          !is.na(firstacc)) %>% nrow()

df_all %>% mutate(rx90.na = is.na(rxhomeday90)) %>% 
  distinct(id, .keep_all=T) %>%
  filter(totaldialdur > 90,
         !is.na(bmi),
         bmi > 15,
         bmi <= 50,
         !is.na(sercreat),
         rx90.na == 0) %>% 
         nrow

df_all %>%
  filter(({{t0}} <= 90) & ({{t}} > 90)) %>% #%>%
  #distinct(id, .keep_all=T) %>%
  filter(totaldialdur > 90,
         !is.na(bmi),
         bmi > 15,
         bmi <= 50,
         !is.na(sercreat),
         !is.na(rxhomeva)) %>% select(id, rxhomeva) %>% nrow

df_all %>% 
  filter(totaldialdur > 90,
         !is.na(bmi),
         bmi > 15,
         bmi <= 50,
         !is.na(sercreat),
         !is.na(rxhomeva),
         !is.na(rxhomeday90),
         !is.na(vacategory90)) %>%
          distinct(id, .keep_all=T) %>% nrow


df_all %>% 
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
          filter(smoke!=4) %>%
  distinct(id, .keep_all=T) %>% nrow






# Exploratory analysis of treatment models

#First considering what lag to have on past treatment

mdl_1 <- with(df %>% filter(.$ti > prm$t_a_vec[[2]]), 
            glm(a_1 ~ 1 + a_0, family=binomial))


mdl_2a <- with(df %>% filter(.$ti > prm$t_a_vec[[2+1]]), 
              glm(a_2 ~ 1 + a_1, family=binomial))
mdl_2b <- with(df %>% filter(.$ti > prm$t_a_vec[[2+1]]), 
               glm(a_2 ~ 1 + a_1 + a_0, family=binomial))
mdl_2c <- with(df %>% filter(.$ti > prm$t_a_vec[[2+1]]), 
               glm(a_2 ~ 1 + a_1 + a_0 + a_1:a_0, family=binomial))
summary(mdl_2a)
summary(mdl_2b)
summary(mdl_2c)



mdl_3a <- with(df %>% filter(.$ti > prm$t_a_vec[[3+1]]), 
               glm(a_3 ~ 1 + a_2, family=binomial))
mdl_3b <- with(df %>% filter(.$ti > prm$t_a_vec[[3+1]]), 
               glm(a_3 ~ 1 + a_2 + a_1, family=binomial))
mdl_3c <- with(df %>% filter(.$ti > prm$t_a_vec[[3+1]]), 
               glm(a_3 ~ 1 + a_2 + a_1 + a_0, family=binomial))
summary(mdl_3a)
summary(mdl_3b)
summary(mdl_3c)




mdl_4a <- with(df %>% filter(.$ti > prm$t_a_vec[[4+1]]), 
               glm(a_4 ~ 1 + a_3, family=binomial))
mdl_4b <- with(df %>% filter(.$ti > prm$t_a_vec[[4+1]]), 
               glm(a_4 ~ 1 + a_3 + a_2, family=binomial))
mdl_4c <- with(df %>% filter(.$ti > prm$t_a_vec[[4+1]]), 
               glm(a_4 ~ 1 + a_3 + a_2 + a_1, family=binomial))
mdl_4d <- with(df %>% filter(.$ti > prm$t_a_vec[[4+1]]), 
               glm(a_4 ~ 1 + a_3 + a_2 + a_1 + a_0, family=binomial))
summary(mdl_4a)
summary(mdl_4b)
summary(mdl_4c)
summary(mdl_4d)






#Now considering baseline covariates
mdl_0 <- with(df %>% filter(.$ti > prm$t_a_vec[[0+1]]), 
              glm(a_1 ~ 1 + a_0 + age + sex + race, family=binomial))
summary(mdl_0)







#Looking at number of treatment changes
changes_vec <- sapply(1:33, df=df, prm=prm, 
                      FUN= function(k, df, prm){
  col_1 <- as.name(paste0("a_", k-1))
  col_2 <- as.name(paste0("a_", k))
  df_k <- df %>% filter(.$ti > prm$t_a_vec[[k+1]]) %>%
                  mutate(a_change = as.integer({{col_1}} == {{col_2}}))
  return(sum(df_k$a_change))
})







model_list <- lapply(0:32, function(k){
  mdl_formula <- paste0("a_", k, " ~ ", 
                        paste0(prm$trt_mod_list[[k+1]], collapse = " + "))
  
  model_temp <- with(df %>% filter(.$ti > prm$t_a_vec[[k+1]]), 
                     glm(  as.formula(mdl_formula),
                           family=binomial))
  
})





q <-  sapply(1:32, 
      function(k){
        var1 <- as.name(paste0("a_", k-1))
        var2 <- as.name(paste0("a_", k))
        df_q <- df %>% filter(.$ti > prm$t_a_vec[[k+1]]) %>% 
                mutate(a_diff = abs({{var2}} - {{var1}})) %>% select(id, {{var1}}, {{var2}}, a_diff)
        return(df_q$a_diff %>% sum())
      })






ggplot(df_0, aes(logit0, age))+
  geom_point(size = 0.5, alpha = 0.5)



p0_S <- tibble(psi0=numeric(), S0=numeric())

q <- lapply(1:10, function(psi0){
          S <- df %>% calculate_score(., psi_hat_vec=c(psi0, 0), prm=prm)
          
          return(tibble(psi0=psi0, S0=S[[1]]))
}) %>% Reduce(add_row, .)




#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           TESTING MODELS FOR BASELINE COVARIATES
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##


df_CvA <- df_CVC_vs_AVFG_wide
df_CvA <- df_CvA %>% left_join(df_baseline, by="id")
df_CvA <- df_CvA %>% mutate(sex=as.factor(sex),
                    race=factor(race, labels = c("NZ maori/pacific" = 1,
                                                 "Aus Indigenous" = 2,
                                                 "Asian" = 3,
                                                 "White/other" = 4)),
                    firstyear = as.factor(format(firstdate,"%Y")),
                    smoke = as.factor(smoke)
)
df_CvA <- df_CvA %>% filter(!is.na(sercreat)) %>% 
                    filter(!latereferral=="") %>%
                    filter(smoke!=4) %>%
                    filter(!is.na(bmi_cat))

model_out <- function(k){
      mdl_formula <- paste0("a_", k, " ~ 1 + age + sex + race + smoke + firstyear +
                sercreat + latereferral + cad + lung + pvd + diabetes +
                bmi_cat + at_home")
      
      if(k > 0){
        mdl_formula <- paste0(mdl_formula, " + a_", k-1)
      }
      
      with(df %>% filter(.$ti > prm$t_a_vec[[k+1]]), 
         glm(as.formula(mdl_formula),
              family=binomial))
}

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
  mdl_k <- c("1", "age", "sex", "race", "smoke", #"firstyear",
             "sercreat", "latereferral", "cad", "lung", "pvd",
             "diabetes", "bmi_cat", "at_home")
  if(k > 0){ mdl_k <- c(mdl_k, paste0("a_", k-1))}
  return(mdl_k)
})
prm$censor <- T

df <- df_CvA %>% filter(firstdate < prm$censor_date)

fit_trt_out <- df %>% fit_treatment_models(prm=prm)
df <- fit_trt_out[[1]]
trt_models <- fit_trt_out[[2]]
df <- df %>% mutate(x = 0) %>%
  mutate(C_i = difftime(prm$censor_date, firstdate,  "days"))

##### EDA looking at how the coefficients change over time
trt_mdls <- trt_models
trt_coeffs <- lapply(1:length(trt_mdls), 
                     function(k){
                        trt_mdls[[k]]$coefficients %>% data.frame()
                     })
trt_coeffs[[1]] <- trt_coeffs[[1]] %>% rbind(a_k=0)
coeffs_df <- trt_coeffs %>% Reduce(cbind, .)


##### Runnings various models
nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                     psi_start_vec=rep(0, max(prm$beta_1_track) + max(prm$beta_x_track)))
(psi_hat_vec <- nr_out[[1]])



prm$beta_1_track <- (c(rep(1,4), rep(2, 26-4)))
nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                     psi_start_vec=rep(0, max(prm$beta_1_track) + max(prm$beta_x_track)))
(psi_hat_vec <- nr_out[[1]])


prm$beta_1_track <- c(rep(1,4), rep(2,4), rep(3,4),
                      rep(4,4), rep(5,4), rep(6,6))
nr_out <- df %>% newton_raphson_grad(prm=prm, tol=0.01,
                                     psi_start_vec=rep(0, max(prm$beta_1_track) + max(prm$beta_x_track)))
(psi_hat <- nr_out[[1]])
(VCV <- df %>% calculate_variance(prm, psi_hat_vec, trt_models))




#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           COUNTING TREATMENT SWITCHES FOR CvA
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##


switch_counter <- lapply(1:24, function(k){
  var_now <- as.name(paste0("a_", k))
  var_prev <- as.name(paste0("a_", k-1))
  
  df <- df_cva_25_full %>% 
        mutate(lastperiod = pmin(lastperiod, (as.integer(C_i) %/% 90)-1)) %>% 
        select({{var_now}}, {{var_prev}}, lastperiod) %>% 
        filter(lastperiod >= k) %>%    
        mutate(switch_0 = ifelse({{var_now}}!={{var_prev}} & {{var_now}}==0, 1, 0),
               switch_1 = ifelse({{var_now}}!={{var_prev}} & {{var_now}}==1, 1, 0),
               dies_0 = ifelse(lastperiod==k & {{var_now}}==0, 1, 0),
               dies_1 = ifelse(lastperiod==k & {{var_now}}==1, 1, 0)
               )
  
  switch_0_amount <- df %>% filter(switch_0==1) %>% nrow()
  switch_1_amount <- df %>% filter(switch_1==1) %>% nrow()
  
  switch_0_deaths <- df %>% filter(switch_0==1) %>% 
                  select(dies_0) %>% sum(na.rm=T)
  switch_1_deaths <- df %>% filter(switch_1==1) %>% 
                      select(dies_1) %>% sum(na.rm=T)
  
  deaths_0 <- df %>% filter(dies_0==1) %>% nrow()
  deaths_1 <- df %>% filter(dies_1==1) %>% nrow()
  total_deaths_check <- df %>% filter(lastperiod==k) %>% nrow()
  
  n_alive_start <- df %>% nrow()
  n_alive_end <- (df %>% nrow()) - total_deaths_check 
  
  return(tibble(period=k,
                swh_0=switch_0_amount,
                sww_1=switch_1_amount,
                sw_tot=switch_0_amount+switch_1_amount,
                
                sw_0_dth = switch_0_deaths,
                sw_1_dth = switch_1_deaths,
                sw_d_tot = switch_0_deaths + switch_1_deaths,
                
                dths_0 = deaths_0,
                dths_1 = deaths_1,
                dths_tot = deaths_0 + deaths_1,
                dths_chk = total_deaths_check,
                n_alv_strt = n_alive_start,
                n_alv_end = n_alive_end
             ))
}) %>% Reduce(rbind, .) %>% print(n=25)

# df_check <- df_cva_25_full %>%
#               mutate(dies_now=ifelse(lastperiod==k, 1, 0)
#                      
#               filter(lastperiod==k)




#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           COUNTING TREATMENT SWITCHES FOR PDvO
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##


switch_counter_PD <- lapply(1:25, function(k){
  var_now <- as.name(paste0("a_", k))
  var_prev <- as.name(paste0("a_", k-1))
  
  df <- df_PD_25 %>% 
    select({{var_now}}, {{var_prev}}, lastperiod) %>% 
    filter(lastperiod >= k) %>%    
    mutate(switch_0 = ifelse({{var_now}}!={{var_prev}} & {{var_now}}==0, 1, 0),
           switch_1 = ifelse({{var_now}}!={{var_prev}} & {{var_now}}==1, 1, 0),
           dies_0 = ifelse(lastperiod==k & {{var_now}}==0, 1, 0),
           dies_1 = ifelse(lastperiod==k & {{var_now}}==1, 1, 0)
    )
  
  switch_0_amount <- df %>% filter(switch_0==1) %>% nrow()
  switch_1_amount <- df %>% filter(switch_1==1) %>% nrow()
  
  switch_0_deaths <- df %>% filter(switch_0==1) %>% 
    select(dies_0) %>% sum(na.rm=T)
  switch_1_deaths <- df %>% filter(switch_1==1) %>% 
    select(dies_1) %>% sum(na.rm=T)
  
  deaths_0 <- df %>% filter(dies_0==1) %>% nrow()
  deaths_1 <- df %>% filter(dies_1==1) %>% nrow()
  total_deaths_check <- df %>% filter(lastperiod==k) %>% nrow()
  
  n_alive_start <- df %>% nrow()
  n_alive_end <- (df %>% nrow()) - total_deaths_check 
  
  return(tibble(period=k,
                swh_0=switch_0_amount,
                sww_1=switch_1_amount,
                sw_tot=switch_0_amount+switch_1_amount,
                
                sw_0_dth = switch_0_deaths,
                sw_1_dth = switch_1_deaths,
                sw_d_tot = switch_0_deaths + switch_1_deaths,
                
                dths_0 = deaths_0,
                dths_1 = deaths_1,
                dths_tot = deaths_0 + deaths_1,
                dths_chk = total_deaths_check,
                n_alv_strt = n_alive_start,
                n_alv_end = n_alive_end
  ))
}) %>% Reduce(rbind, .) %>% print(n=25)

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           BUILDING MODELS FOR MA
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##


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



#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##
#           BUILDING MODELS FOR PERIOD SPECIFIC PSIS
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   ##


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

################################
### INVERSE VARIANCE WEIGHTING
#################################
var_indiv <- ste_indiv^2 %>% select(-c(period))
df_indiv_less <- df_indiv %>% select(-c(period))

ivwa <- lapply(1:7, 
  function(k){
    sum_inv_var <- (1/var_indiv[k, (k+1):7]) %>% sum()  
    mean_inv_weight <- ((df_indiv[k, (k+1):7] / var_indiv[k, (k+1):7]) %>% sum()) /  sum_inv_var
        tibble(period=k-1, 
            psi_avg = mean_inv_weight,
            var_avg = 6/sum_inv_var,
            ste_avg = sqrt(6/sum_inv_var))
  }) %>% Reduce(rbind, .) #%>% ggplot(aes(x=period, y=psi_avg)) + geom_point() #+ ylim(0,2)

ivwa 
#%>% Reduce(rbind, .) %>% ggplot(aes(x=period, y=psi_avg)) + geom_point()






psi_parts_vec <- c(1, 5, 10, 15, 20)
d=2

test_out <- tibble()

test_out <- lapply(psi_parts_vec,
       function(psi_part){
         tibble_row <- tibble(psi=psi_part)
         for (a_curr in (0:(length(prm$t_a_vec)-1))){
           if (prm$beta_1_track[[a_curr+1]] == d) {
             
             df_temp <- df %>% calculate_tau_k(prm=prm, 
                                               psi_hat_vec=c(0, -psi_part),
                                               a_k=a_curr)
             
             a_k <- as.name(paste0("a_", a_curr))
             fit_k <- as.name(paste0("fit_", a_curr))
             
             df_temp <- df_temp %>% mutate(
               S_i := ifelse((is.na({{a_k}}) | is.na({{fit_k}}) | is.na(tau_k)), 0,
                             ({{a_k}} - {{fit_k}})*tau_k)
             )
             
             S_inc <- df_temp %>% select(S_i) %>% sum()
             
             tibble_row <- tibble_row %>% add_column("{a_k}":= S_inc)
           }
           
         }
         return(tibble_row)
       }) %>% Reduce(rbind, .)

test_out %>% pivot_longer(cols = !psi) %>% 
  ggplot(aes(x=psi, y=value, group=name)) +
  geom_point() + geom_line()














