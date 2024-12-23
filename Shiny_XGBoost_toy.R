library(xgboost)
library(tidyverse)

## Helper function: calculate permutation p-value
xgb_permu_p_val = function(x_df, label, bst, obs_imp_df, nrounds, num_permutations, early_stopping = TRUE){
  sample_size = nrow(x_df)
  tmp = data.frame(Feature = colnames(x_df))
  obs_imp_df = left_join(tmp, obs_imp_df, by = "Feature")
  obs_gain_imp = ifelse(is.na(obs_imp_df$Gain), 0, obs_imp_df$Gain)
  obs_cover_imp = ifelse(is.na(obs_imp_df$Cover), 0, obs_imp_df$Cover)
  
  extreme_gain_cnt = 0
  extreme_cover_cnt = 0
  for(i in 1:num_permutations){
    shuffle_label = sample(label, size = sample_size, replace = FALSE)
    shuffle_xgb_df = xgb.DMatrix(data = x_df, label = shuffle_label)
    if(early_stopping){
      shuffle_bst = xgboost(data = shuffle_xgb_df, max.depth = bst$params$max_depth, nrounds = bst$best_iteration, objective = bst$params$objective, verbose = 0)
    }
    else{
      shuffle_bst = xgboost(data = shuffle_xgb_df, max.depth = bst$params$max_depth, nrounds = nrounds, objective = bst$params$objective, verbose = 0)
    }
    shuffle_imp_df = xgb.importance(model = shuffle_bst)
    shuffle_imp_df = left_join(tmp, shuffle_imp_df, by = "Feature")
    shuffle_gain_imp = ifelse(is.na(shuffle_imp_df$Gain), 0, shuffle_imp_df$Gain)
    shuffle_cover_imp = ifelse(is.na(shuffle_imp_df$Cover), 0, shuffle_imp_df$Cover)
    
    extreme_gain_cnt = extreme_gain_cnt + as.integer(shuffle_gain_imp >= obs_gain_imp)
    extreme_cover_cnt = extreme_cover_cnt + as.integer(shuffle_cover_imp >= obs_cover_imp)
  }
  gain_res = data.frame(Feature = colnames(x_df), Gain = obs_gain_imp, Gain_pval = (extreme_gain_cnt + 1) / (num_permutations + 1))
  cover_res = data.frame(Feature = colnames(x_df), Cover = obs_cover_imp, Cover_pval = (extreme_cover_cnt + 1) / (num_permutations + 1))
  return(list(gain = gain_res, cover = cover_res))
}

## Input
## 1. df: the input data
## 2. outcome_var: string, the name of outcome variable in the data;
## time_var, status_var: string, the name of time and status variables for survival data
## 3. dt_type: lm/glm/cox
## 4. max.depth: max depth of each tree in XGBoost
## 5. nrounds: number of iterations in XGBoost
## 6. early_stopping_rounds: number of rounds for early stopping
## 7. npermu: number of iterations for permutation test
## 8. alpha: type-1 error for permutation test
## 9. seed: random seed for replication
## -----------------------------------------------------------------------------
## Output: a dataframe including the following columns
## xnames: names of explanatory variables
## runtime: running time of the variable selection process
## importance: variable importance (the higher, the more important)
## pval: p-value of each variable based on permutation test
## selection: whether the algorithm selects each variable (1: selected; 0: not selected)
## -----------------------------------------------------------------------------
## Notes
## To test your code, you may reduce max.depth, nrounds, and/or npermu
var_select_XGBoost = function(df, outcome_var=NULL, time_var=NULL, status_var=NULL,
                           dt_type="lm", max.depth = 6, nrounds = 500,
                           early_stopping_rounds = 10, npermu = 100,
                           alpha = 0.05, seed=1024){
  set.seed(seed)
  if(dt_type == "lm"){
    label = df[, outcome_var]
    x_df = as.matrix(df[, colnames(df) != outcome_var])
    xgb_df = xgb.DMatrix(data = x_df, label = label)
    
    tstart = proc.time()
    if(!is.null(early_stopping_rounds)){
      bst = xgboost(data = xgb_df, max.depth = max.depth, nrounds = nrounds,
                    objective = "reg:squarederror", 
                    early_stopping_rounds = early_stopping_rounds, verbose = 0)
      obs_imp_df = xgb.importance(model = bst)
      
      res = xgb_permu_p_val(x_df, label, bst, obs_imp_df, nrounds = nrounds,
                            num_permutations = npermu, early_stopping = TRUE)
    }
    else{
      bst = xgboost(data = xgb_df, max.depth = max.depth, nrounds = nrounds,
                    objective = "reg:squarederror", 
                    verbose = 0)
      obs_imp_df = xgb.importance(model = bst)
      
      res = xgb_permu_p_val(x_df, label, bst, obs_imp_df, nrounds = nrounds,
                            num_permutations = npermu, early_stopping = FALSE)
    }
    runtime = (proc.time() - tstart)[[3]]
    
    # the important variables are those whose pvalue < alpha
    # selection=1 indicates that a variable is import
    res$gain$selection = ifelse(res$gain$Gain_pval < alpha, 1, 0)
    
    selection_res = data.frame(xnames=res$gain$Feature,
                               runtime=runtime,
                               importance=res$gain$Gain,
                               pval=res$gain$Gain_pval,
                               selection=res$gain$selection,
                               row.names = NULL)
  }
  else if (dt_type == "glm"){
    label = df[, outcome_var]
    x_df = as.matrix(df[, colnames(df) != outcome_var])
    xgb_df = xgb.DMatrix(data = x_df, label = label)
    
    tstart = proc.time()
    if(!is.null(early_stopping_rounds)){
      bst = xgboost(data = xgb_df, max.depth = max.depth, nrounds = nrounds,
                    objective = "binary:logistic", 
                  early_stopping_rounds = early_stopping_rounds, verbose = 0)
      obs_imp_df = xgb.importance(model = bst)
      
      res = xgb_permu_p_val(x_df, label, bst, obs_imp_df, nrounds = nrounds,
                            num_permutations = npermu, early_stopping = TRUE)
    }
    else{
      bst = xgboost(data = xgb_df, max.depth = max.depth, nrounds = nrounds,
                    objective = "binary:logistic", 
                    verbose = 0)
      obs_imp_df = xgb.importance(model = bst)
      
      res = xgb_permu_p_val(x_df, label, bst, obs_imp_df, nrounds = nrounds,
                            num_permutations = npermu, early_stopping = FALSE)
      
    }
    runtime = (proc.time() - tstart)[[3]]
    
    # the important variables are those whose pvalue < alpha
    # selection=1 indicates that a variable is import
    res$gain$selection = ifelse(res$gain$Gain_pval < alpha, 1, 0)
    
    selection_res = data.frame(xnames=res$gain$Feature,
                               runtime=runtime,
                               importance=res$gain$Gain,
                               pval=res$gain$Gain_pval,
                               selection=res$gain$selection,
                               row.names = NULL)
  }
  else{
    label = ifelse(df[, status_var] == 1, df[, time_var], -df[, time_var])
    x_df = as.matrix(df[, !colnames(df) %in% c(status_var, time_var)])
    xgb_df = xgb.DMatrix(data = x_df, label = label)
    
    tstart = proc.time()
    if(!is.null(early_stopping_rounds)){
      bst = xgboost(data = xgb_df, max.depth = max.depth, nrounds = nrounds,
                    objective = "survival:cox", 
                    early_stopping_rounds = early_stopping_rounds, verbose = 0)
      obs_imp_df = xgb.importance(model = bst)
      
      res = xgb_permu_p_val(x_df, label, bst, obs_imp_df, nrounds = nrounds,
                            num_permutations = npermu, early_stopping = TRUE)
    }
    else{
      bst = xgboost(data = xgb_df, max.depth = max.depth, nrounds = nrounds,
                    objective = "survival:cox", 
                    verbose = 0)
      obs_imp_df = xgb.importance(model = bst)
      
      res = xgb_permu_p_val(x_df, label, bst, obs_imp_df, nrounds = nrounds,
                            num_permutations = npermu, early_stopping = FALSE)
    }
    runtime = (proc.time() - tstart)[[3]]
    
    # the important variables are those whose pvalue < alpha
    # selection=1 indicates that a variable is import
    res$gain$selection = ifelse(res$gain$Gain_pval < alpha, 1, 0)
    
    selection_res = data.frame(xnames=res$gain$Feature,
                               runtime=runtime,
                               importance=res$gain$Gain,
                               pval=res$gain$Gain_pval,
                               selection=res$gain$selection,
                               row.names = NULL)
  }
  return(selection_res)
}