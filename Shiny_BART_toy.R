library(bartMachine)
library(tidyverse)

## Input
## 1. df: the input data
## 2. outcome_var: string, the name of outcome variable in the data
## 3. dt_type: lm/glm (BART does not support survival outcomes due to the limitation of the R package)
## 4. nburn: number of burn in rounds of MCMC
## 5. niter: number of rounds for sampling from the posterior
## 6. npermu: number of rounds for the permutation test 
## 7. alpha: type-1 error for permutation test based variable selection
## 8. ntree: number of trees in the algorithm
## 9. seed: random seed for replication
## -----------------------------------------------------------------------------
## Output: a dataframe including the following columns
## xnames: names of explanatory variables
## runtime: running time of the variable selection process
## importance: variable importance (scale of (0,1), the higher, the more important)
## pval: p-value of each variable based on permutation test
##   note that there are 3 different methods to calculate the p-value, including
##   local, global max, and global SE. So here we calculate all 3 types of p-values
## selection: whether the algorithm selects each variable (1: selected; 0: not selected)
##   We have 3 different types of p-values, so we can obtain 3 different selection results
## -----------------------------------------------------------------------------
## Notes
## It is normal that this method is slow
## To test your code, you may reduce nburn, niter, and/or npermu
var_select_BART = function(df, outcome_var=NULL,
                           dt_type="lm", nburn = 2000, niter = 5000, npermu = 100,
                           alpha = 0.05, ntree = 20, seed=1024){
  set.seed(seed)
  df_x = df[, names(df)[which(names(df) != outcome_var)]]
  df_y = df[, outcome_var]
  if(dt_type == "glm"){
    df_y = as.factor(df_y)
  }
  
  tstart = proc.time()
  # train a BART model
  invisible(capture.output(bart <- bartMachine(X = df_x, y = df_y, num_burn_in = nburn, 
                                               num_iterations_after_burn_in = niter, 
                                               num_trees = ntree, run_in_sample = FALSE,
                                               verbose = FALSE)))
  tmp = data.frame(feature = colnames(bart$X))
  
  invisible(capture.output(bart_var_sel <- var_selection_by_permute(bart, num_reps_for_avg = 10, 
                                                                    num_permute_samples = npermu,
                                                                    num_trees_for_permute = ntree, 
                                                                    alpha = alpha, plot = FALSE)))
  runtime = (proc.time() - tstart)[3]
  
  # variable importance, by the descending order
  bart_vip = bart_var_sel$var_true_props_avg
  permute_mat = bart_var_sel$permute_mat
  bart_vip_df = data.frame(feature = names(bart_vip), vip = bart_vip)
  tmp = left_join(tmp, bart_vip_df, by = "feature")
  
  # calculate permutation p-values
  # first type of p-value: local method
  pval_local_df = data.frame(feature = names(bart_vip), pval_local = sapply(1:ncol(permute_mat), function(j){
    mean(permute_mat[, j] >= bart_vip[j])
  }))
  tmp = left_join(tmp, pval_local_df, by = "feature")
  # selection=1 indicates that a variable is import
  tmp$selection_local = 0
  tmp$selection_local[bart_var_sel$important_vars_local_col_nums] = 1
  
  # second type of p-value: global max method
  vip_max = apply(permute_mat, MARGIN = 1, max)
  pval_global_max_df = data.frame(feature = names(bart_vip), pval_global_max = sapply(bart_vip, function(x){
    mean(vip_max >= x)
  }))
  tmp = left_join(tmp, pval_global_max_df, by = "feature")
  # selection=1 indicates that a variable is import
  tmp$selection_global_max = 0
  tmp$selection_global_max[bart_var_sel$important_vars_global_max_col_nums] = 1
  
  # third type of p-value: global SE method
  mk = colMeans(permute_mat)
  sk = apply(permute_mat, MARGIN = 2, sd)
  C_star = (bart_vip - mk)/sk
  pval_global_se_df = data.frame(feature = names(bart_vip), pval_global_se = sapply(C_star, function(C){
    1 - mean(sapply(1:nrow(permute_mat), function(x){
      all(permute_mat[x, ] <= mk + sk * C)
    }))
  }))
  tmp = left_join(tmp, pval_global_se_df, by = "feature")
  # selection=1 indicates that a variable is import
  tmp$selection_global_se = 0
  tmp$selection_global_se[bart_var_sel$important_vars_global_se_col_nums] = 1
  
  selection_res = data.frame(xnames=tmp$feature,
                             runtime=runtime,
                             importance=tmp$vip,
                             pval_local=tmp$pval_local,
                             pval_global_max=tmp$pval_global_max,
                             pval_global_se=tmp$pval_global_se,
                             selection_local=tmp$selection_local,
                             selection_global_max=tmp$selection_global_max,
                             selection_global_se=tmp$selection_global_se,
                             row.names = NULL)
  
  return(selection_res)
}