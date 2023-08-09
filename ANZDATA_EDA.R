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
