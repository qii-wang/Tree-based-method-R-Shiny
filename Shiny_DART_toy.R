library(BART)

## Input
## 1. df: the input data
## 2. outcome_var: string, the name of outcome variable in the data;
## time_var, status_var: string, the name of time and status variables for survival data
## 3. dt_type: lm/glm/cox
## 4. nskip: number of burn in rounds of MCMC
## 5. ndpost: number of rounds for sampling from the posterior
## 6. keepevery_cox: for the survival outcome, it is recommended to only keep
## the sample from posterior every keepevery rounds
## 7. thresh_dart: the threshold for variable selection
## 8. ntree: number of trees in the algorithm
## 9. seed: random seed for replication
## -----------------------------------------------------------------------------
## Output: a dataframe including the following columns
## xnames: names of explanatory variables
## runtime: running time of the variable selection process
## importance: variable importance (scale of (0,1), the higher, the more important)
## selection: whether the algorithm selects each variable (1: selected; 0: not selected)
## -----------------------------------------------------------------------------
## Notes
## It is normal that dealing with survival outcomes is slow
## To test your code, you may reduce nskip_cox, ndpost_cox, and/or keepevery_cox
var_select_DART = function(df, outcome_var=NULL, time_var=NULL, status_var=NULL,
                           dt_type="lm", nskip_lm = 2000, ndpost_lm = 5000,
                           nskip_glm = 2000, ndpost_glm = 5000,
                           nskip_cox = 2000, ndpost_cox = 5000, keepevery_cox = 10, 
                           thresh_dart = 0.5, ntree = 20, seed=1024){
  set.seed(seed)
  if(dt_type == "lm"){
    df_x = df[, names(df)[which(names(df) != outcome_var)]]
    tstart = proc.time()
    invisible(capture.output(dart <- wbart(x.train = df_x, y.train = df[, outcome_var],
                                           sparse = TRUE, ntree = ntree, 
                                           printevery = 1e5,
                                           ndpost = ndpost_lm, nskip = nskip_lm)))
    runtime = (proc.time() - tstart)[[3]]
    
    mpvip_dart = apply(dart$varcount, MARGIN = 2, function(x){
      mean(x>0)
    })
    
    # the important variables are those whose MPVIP > thresh_dart
    # selection_dart=1 indicates that a variable is import
    selection_dart = as.integer(mpvip_dart > thresh_dart)
    selection_res = data.frame(xnames=names(df_x),
                               runtime=runtime,
                               importance=mpvip_dart,
                               selection=selection_dart,
                               row.names = NULL)
    
  }
  
  else if (dt_type == "glm"){
    df_x = df[, names(df)[which(names(df) != outcome_var)]]
    
    tstart = proc.time()
    invisible(capture.output(dart <- pbart(x.train = df_x, y.train = df[, outcome_var],
                                           sparse = TRUE, ntree = ntree, 
                                           printevery = 1e5,
                                           ndpost = ndpost_glm, nskip = nskip_glm)))
    runtime = (proc.time() - tstart)[[3]]
    
    mpvip_dart = apply(dart$varcount, MARGIN = 2, function(x){
      mean(x>0)
    })
    
    # the important variables are those whose MPVIP > thresh_dart
    # selection_dart=1 indicates that a variable is import
    selection_dart = as.integer(mpvip_dart > thresh_dart)
    
    selection_res = data.frame(xnames=names(df_x),
                               runtime=runtime,
                               importance=mpvip_dart,
                               selection=selection_dart,
                               row.names = NULL)
  }
  
  else{
    df_x = df[, names(df)[which(!(names(df) %in% c(time_var, status_var)))]]
    
    tstart = proc.time()
    invisible(capture.output(dart <- surv.bart(x.train = df_x, times = df[, time_var],
                                               delta = df[, status_var], sparse = TRUE,
                                               ntree = ntree, 
                                               printevery = 1e5,
                                               ndpost = ndpost_cox, keepevery = keepevery_cox,
                                               nskip = nskip_cox)))
    runtime = (proc.time() - tstart)[[3]]
    
    mpvip_dart = apply(dart$varcount[, -1], MARGIN = 2, function(x){
      mean(x>0)
    })
    
    # the important variables are those whose MPVIP > thresh_dart
    # selection_dart=1 indicates that a variable is import
    selection_dart = as.integer(mpvip_dart > thresh_dart)
    
    selection_res = data.frame(xnames=names(df_x),
                               runtime=runtime,
                               importance=mpvip_dart,
                               selection=selection_dart,
                               row.names = NULL)
  }
  return(selection_res)
}