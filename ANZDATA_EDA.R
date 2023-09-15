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
