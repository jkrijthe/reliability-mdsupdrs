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

# system.time(models <- df_datasets %>% 
#               mutate(Model = map(Dataset, apply_model, df_models=df_models)) %>% 
#               unnest(Model,.drop=FALSE) %>% 
#               mutate(params=as.character(params)))


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


# Bootstrap for comparing stricter criterion to strict criterion
apply_bootstrap_compare <- function(df,df_models,reps=100) {
  patnos <- df$PATNO %>% unique
  months <- df$Month %>% unique %>% sort
  intervals <- months-c(-3,months[-length(months)]) # Assume start is 3 months before first measurement
  

  
  mclapply(1:reps, function(i) {
    print(i);
    selected_patnos <- sample(patnos,replace=TRUE)
    
    mat_in <- df %>%
      filter(PATNO %in% selected_patnos) %>% 
      select(PATNO,Month,Score) %>% 
      spread(., Month, Score) %>%
      select(-PATNO) %>%
      as.matrix(.)
    
    mat_stricter <- df %>%
      filter(PATNO %in% selected_patnos, strictly==TRUE) %>%
      select(PATNO,Month,Score) %>% 
      spread(., Month, Score) %>%
      select(-PATNO) %>%
      as.matrix(.)
    
    tryCatch(
      bind_rows(
        Strict=df_models  %>% mutate(rep=i, Model = invoke_map(f,params,Y=mat_in,weight=intervals)),
        Stricter=df_models  %>% mutate(rep=i, Model = invoke_map(f,params,Y=mat_stricter,weight=intervals)),
        .id="Level"),
      error=function(x) { data.frame() }
    ) %>% bind_rows
  },mc.cores=30)
}

df_all_measurements <- df_yearly_long %>% filter(Measurement=="MDS_UPDRS_III_off") %>% select(-Measurement)
models_bootstrapped_compare <- apply_bootstrap_compare(df_all_measurements,df_models)

save(models_bootstrapped_compare, file="dlm_estimates_ppmi-compare.RData")