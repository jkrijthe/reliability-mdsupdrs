library(tidyverse,warn.conflicts = FALSE)
library(broom)
library(parallel)
source("R/dlmTrend.R")

set.seed(1)
load("bootstrap-data.RData")


# Make selection of data
df_datasets <- bind_rows(Yearly=df_yearly_long,
                         All=df_long, 
                         Factors=df_factors_long,
                         .id = "Criterion")%>%
  group_by(Measurement,Criterion) %>% 
  nest(.key = "Dataset")

df_models <- tibble(
  f = rep("dlmTrend_mle",1),
  params = list(list(GG=1))
)

# Function to apply dlm model to (long) data.frame with observations
apply_model <- function(df,df_models) {
  mat_in <- df %>%
    select(PATNO,Month,Score) %>% 
    #distinct(PATNO,Month,.keep_all=TRUE) %>% 
    spread(., Month, Score) %>%
    select(-PATNO) %>%
    as.matrix(.)

  months <- df$Month %>% unique %>% sort
  intervals <- months-c(-3,months[-length(months)]) # Assume start is 3 months before first measurement

  df_models  %>% 
    mutate(Model = invoke_map(f,params,Y=mat_in,weight=intervals))
}

print(system.time(models <- df_datasets %>%
              mutate(Model = map(Dataset, apply_model, df_models=df_models)) %>%
              unnest(Model,.drop=FALSE) %>%
              mutate(params=as.character(params))))


# Bootstrapping
apply_bootstrap <- function(df,df_models,reps=1000) {
  mat_in <- df %>%
    select(PATNO,Month,Score) %>% 
    #distinct(PATNO,Month,.keep_all=TRUE) %>% 
    spread(., Month, Score) %>%
    select(-PATNO) %>%
    as.matrix(.)
  
  months <- df$Month %>% unique %>% sort
  intervals <- months-c(-3,months[-length(months)]) # Assume start is 3 months before first measurement
  mclapply(1:reps, function(i) {
    print(i);
    tryCatch(df_models  %>% 
               mutate(rep=i, Model = invoke_map(f,params,Y=mat_in[sample(nrow(mat_in),replace=TRUE),],weight=intervals)),
             error=function(x){data.frame()}) %>% bind_rows
  },mc.cores=30)
}


models_bootstrapped <- df_datasets %>%
  #{map(.$Dataset, apply_bootstrap, df_models=df_models)})
  mutate(Model = map(Dataset, apply_bootstrap, df_models=df_models))
save(df_models, models_bootstrapped, file="dlm_estimates_ppmi.RData")
