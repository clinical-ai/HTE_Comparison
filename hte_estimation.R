hte_estimation = function(scenario, seed) {
  
  library(glmnet)
  library(BART)
  library(grf)
  # library(causalLearning)  # Package for causal boosting and causal MARS -- under development
  wd = getwd()
  setwd(wd)
  source('dgp.R')
  
  #-----------------------------------------------------------------------------------------------------
  # Simulate and format data
  #-----------------------------------------------------------------------------------------------------
  df = dgp(scenario, seed)
  df = df[complete.cases(df),]
  set.seed(seed)
  cap = 10000
  if (nrow(df)>cap) { 
    df = df[sample(1:nrow(df), cap, replace = F),]
  }
  zName = c('high_vol', 'endo_revasc', 'apixaban', 'NPH')[scenario]
  xNames = setdiff(names(df), c(zName, 'y_fact', 'tau', 'pi_fact'))
  n_samp = nrow(df)
  
  #-----------------------------------------------------------------------------------------------------
  # Null estimator
  #-----------------------------------------------------------------------------------------------------
  ATE_naive = mean(df$pi_fact[df[,zName]==1]) - mean(df$pi_fact[df[,zName]==0])
  CATE_true = df$tau
  ATE_true = mean(CATE_true)
  ATE_bias_null = abs( (ATE_naive - ATE_true) / ATE_true * 100 )
  CATE_RMSE_null = sqrt(mean((0 - CATE_true)^2))
  
  #-----------------------------------------------------------------------------------------------------
  # Regularized Logistic Regression estimator (LASSO + RIDGE)
  #-----------------------------------------------------------------------------------------------------
  ptm = proc.time()
  ## Selection of marginal effects with LASSO
  x = data.matrix(cbind(df[,xNames], Z = df[,zName]))
  y = as.factor(df$y_fact)
  penal = rep(1, ncol(x))
  penal[which(colnames(x) %in% c('Z'))] = 0    # impose Exposure variable
  LASSO_fit = cv.glmnet(x, y, family = "binomial", alpha = 1,
                        nfolds = 10, type.measure = "deviance",
                        penalty.factor = penal)
  coef.lasso = coef(LASSO_fit)
  print( impvarNames = colnames(x)[(which(coef.lasso != 0)[-1]-1)] )
  
  if (length(impvarNames)==1) {
    CATE_RMSE_RLR = NA
  } else { 
    ## Interaction model with RIDGE
    fml = as.formula(' ~ .^2')
    x_reduced =  x[,which(colnames(x) %in% impvarNames)]
    x_interac = cbind( model.matrix(fml, data.frame(x_reduced))[,-1])
    GLMNET_fit = cv.glmnet(x_interac, y, family = "binomial", alpha = 0,
                           nfolds = 10, type.measure = "deviance")
    ## CATE estimates
    x_Z0 = x_Z1 = x
    x_Z0[,'Z'] = 0
    x_Z1[,'Z'] = 1
    x_reduced =  x_Z0[,which(colnames(x_Z0) %in% impvarNames)]
    x_interac_Z0 = cbind( model.matrix(fml, data.frame(x_reduced))[,-1])
    x_reduced =  x_Z1[,which(colnames(x_Z1) %in% impvarNames)]
    x_interac_Z1 = cbind( model.matrix(fml, data.frame(x_reduced))[,-1])
    pi_hat_Z0 = plogis(predict(GLMNET_fit, x_interac_Z0))
    pi_hat_Z1 = plogis(predict(GLMNET_fit, x_interac_Z1))
    CATE_hat = as.vector(pi_hat_Z1 - pi_hat_Z0)
    ATE_hat = mean(CATE_hat)
    ATE_bias_RLR = abs( (ATE_hat - ATE_true) / ATE_true * 100 )
    CATE_RMSE_RLR = sqrt(mean((CATE_hat - CATE_true)^2))
  }
  
  ## Extract run time
  rtm_RLR = (proc.time() - ptm)['elapsed']
  
  ## Cleaning
  rm(x,x_Z0,x_Z1,x_reduced,x_interac,x_interac_Z0,x_interac_Z1,y)
  
  #-----------------------------------------------------------------------------------------------------
  # BART (with and without PS-adjustment)
  #-----------------------------------------------------------------------------------------------------
  ## Hyperparameters -- defaults except the following MCMC parameters:
  nskip = 500
  ndpost = 1000
  
  ## BART model for Z
  ptm = proc.time()
  X = df[,xNames]
  y = df[,zName]
  BART_fit_Z = pbart(X, y,
                     nskip = nskip, ndpost = ndpost)
  PS = pnorm(BART_fit_Z$yhat.train.mean)
  rtm_BART_Z = (proc.time() - ptm)['elapsed']
  
  ## BART model for Y
  ptm = proc.time()
  X_fact = cbind(df[,xNames], Z = df[,zName])
  X_count = X_fact
  X_count$Z = abs(X_fact$Z - 1)
  y = df$y_fact
  BART_fit_Y = pbart(x.train = X_fact, y.train = y, x.test = X_count,
                     nskip = nskip, ndpost = ndpost)
  rtm_BART_Y = (proc.time() - ptm)['elapsed']
  
  ## BART model for Y (PS adjusted)
  ptm = proc.time()
  X_fact = cbind(df[,xNames], Z = df[,zName], PS)
  X_count = X_fact
  X_count$Z = abs(X_fact$Z - 1)
  y = df$y_fact
  psBART_fit_Y = pbart(x.train = X_fact, y.train = y, x.test = X_count,
                       nskip = nskip, ndpost = ndpost)
  rtm_BARTps_Y = (proc.time() - ptm)['elapsed']
  
  ## CATE estimates with BART
  ptm = proc.time()
  pi_hat_Z0 = pnorm(BART_fit_Y$yhat.train)
  pi_hat_Z0[,X_fact$Z==1] = pnorm(BART_fit_Y$yhat.test)[ ,X_fact$Z==1]
  pi_hat_Z1 = pnorm(BART_fit_Y$yhat.train)
  pi_hat_Z1[,X_fact$Z==0] = pnorm(BART_fit_Y$yhat.test)[ ,X_fact$Z==0]
  CATE_post = pi_hat_Z1 - pi_hat_Z0
  ATE_post  = apply(CATE_post, 1, mean)
  ATE_bias_post = abs( (ATE_post - ATE_true) / ATE_true * 100 )
  ATE_bias_BART0 = mean(ATE_bias_post)
  ATE_CI = quantile(ATE_post, c(0.025, 0.975))
  ATE_IL_BART0 = ATE_CI[2] - ATE_CI[1]
  ATE_coverage_BART0 = as.numeric(ifelse(ATE_true>=ATE_CI[1] & ATE_true<=ATE_CI[2], 1, 0))
  CATE_hat = apply(CATE_post, 2, mean)
  CATE_RMSE_BART0 = sqrt(mean((CATE_hat - CATE_true)^2))
  rtm_BART_pred = (proc.time() - ptm)['elapsed']
  
  ## CATE estimates with PS-BART
  ptm = proc.time()
  pi_hat_Z0 = pnorm(psBART_fit_Y$yhat.train)
  pi_hat_Z0[,X_fact$Z==1] = pnorm(psBART_fit_Y$yhat.test)[ ,X_fact$Z==1]
  pi_hat_Z1 = pnorm(psBART_fit_Y$yhat.train)
  pi_hat_Z1[,X_fact$Z==0] = pnorm(psBART_fit_Y$yhat.test)[ ,X_fact$Z==0]
  CATE_post = pi_hat_Z1 - pi_hat_Z0
  ATE_post  = apply(CATE_post, 1, mean)
  ATE_bias_post = abs( (ATE_post - ATE_true) / ATE_true * 100 )
  ATE_bias_BART1 = mean(ATE_bias_post)
  ATE_CI = quantile(ATE_post, c(0.025, 0.975))
  ATE_IL_BART1 = ATE_CI[2] - ATE_CI[1]
  ATE_coverage_BART1 = as.numeric(ifelse(ATE_true>=ATE_CI[1] & ATE_true<=ATE_CI[2], 1, 0))
  CATE_hat = apply(CATE_post, 2, mean)
  CATE_RMSE_BART1 = sqrt(mean((CATE_hat - CATE_true)^2))
  rtm_BARTps_pred = (proc.time() - ptm)['elapsed']
  
  ## Extract run times
  rtm_BART0 = rtm_BART_Y + rtm_BART_pred
  rtm_BART1 = rtm_BART_Z + rtm_BARTps_Y + rtm_BARTps_pred
  rm(rtm_BART_Y,rtm_BART_pred,rtm_BART_Z,rtm_BARTps_Y,rtm_BARTps_pred)
  
  ## Cleaning
  rm(X,X_fact,X_count,y)
  
  #-----------------------------------------------------------------------------------------------------
  # Causal Forest
  #-----------------------------------------------------------------------------------------------------
  ptm = proc.time()
  
  X = as.matrix(df[,xNames])
  W = df[,zName]
  Y = df$y_fact
  CF_fit = causal_forest(X, Y, W, seed = seed)
  CATE_hat = as.vector(predict(CF_fit)$predictions)
  CATE_RMSE_CF = sqrt(mean((CATE_hat - CATE_true)^2))
  ATE_estimates = estimate_average_effect(CF_fit, target.sample = 'all', method = 'AIPW')
  ATE_hat = ATE_estimates[1]
  ATE_bias_CF = as.numeric( abs( (ATE_hat - ATE_true) / ATE_true * 100 ) )
  ATE_CI = c(ATE_hat-1.96*ATE_estimates[2], ATE_hat+1.96*ATE_estimates[2])
  ATE_IL_CF = ATE_CI[2] - ATE_CI[1]
  ATE_coverage_CF = as.numeric(ifelse(ATE_true>=ATE_CI[1] & ATE_true<=ATE_CI[2], 1, 0))
  
  ## Extract run time
  rtm_CF = (proc.time() - ptm)['elapsed']
  
  ## Cleaning
  rm(X,W,Y)
  
  #-----------------------------------------------------------------------------------------------------
  # Causal Boosting and Causal MARS (with and without PS-adjustment).
  # The following code was written for a rough beta version of the 'causalLearning' package and might have 
  # to be modified for the published version.
  #-----------------------------------------------------------------------------------------------------
  # X = as.matrix(df[,xNames])
  # Z = df[,zName]
  # Y = df$y_fact
  # 
  # ## CB
  # ptm = proc.time()
  # ps_adjus = FALSE
  # set.seed(seed)
  # CB_fit = cv.causalBoosting(X, Z, Y,
  #                            num.trees = 300, maxleaves = 4, splitSpread = 0.2,
  #                            type.measure = 'effect', nfolds = 3,
  #                            propensity = ps_adjus)
  # CATE_hat = predict(CB_fit, X)
  # ATE_hat = mean(CATE_hat)
  # ATE_bias_CB0 = as.numeric( abs( (ATE_hat - ATE_true) / ATE_true * 100 ) )
  # CATE_RMSE_CB0 = sqrt(mean((CATE_hat - CATE_true)^2))
  # rtm_CB0 = (proc.time() - ptm)['elapsed']
  # 
  # ## BCM
  # ptm = proc.time()
  # set.seed(seed)
  # BCM_fit = bagged.causalMARS(X, Z, Y, propensity = ps_adjus)
  # CATE_hat = predict(BCM_fit, X)
  # ATE_hat = mean(CATE_hat)
  # ATE_bias_BCM0 = as.numeric( abs( (ATE_hat - ATE_true) / ATE_true * 100 ) )
  # CATE_RMSE_BCM0 = sqrt(mean((CATE_hat - CATE_true)^2))
  # rtm_BCM0 = (proc.time() - ptm)['elapsed']
  # 
  # ## PS-CB
  # ptm = proc.time()
  # ps_adjus = TRUE
  # stratum = stratify(PS, Z)$stratum
  # set.seed(seed)
  # CB_fit = cv.causalBoosting(X, Z, Y,
  #                            num.trees = 300, maxleaves = 4, splitSpread = 0.2,
  #                            type.measure = 'effect', nfolds = 3,
  #                            propensity = ps_adjus, stratum = stratum)
  # CATE_hat = predict(CB_fit, X)
  # ATE_hat = mean(CATE_hat)
  # ATE_bias_CB1 = as.numeric( abs( (ATE_hat - ATE_true) / ATE_true * 100 ) )
  # CATE_RMSE_CB1 = sqrt(mean((CATE_hat - CATE_true)^2))
  # rtm_CB1 = (proc.time() - ptm)['elapsed']
  # 
  # ## PS-BCM
  # ptm = proc.time()
  # set.seed(seed)
  # BCM_fit = bagged.causalMARS(X, Z, Y, propensity = ps_adjus, stratum = stratum)
  # CATE_hat = predict(BCM_fit, X)
  # ATE_hat = mean(CATE_hat)
  # ATE_bias_BCM1 = as.numeric( abs( (ATE_hat - ATE_true) / ATE_true * 100 ) )
  # CATE_RMSE_BCM1 = sqrt(mean((CATE_hat - CATE_true)^2))
  # rtm_BCM1 = (proc.time() - ptm)['elapsed']
  
  #-----------------------------------------------------------------------------------------------------
  # Save results
  #-----------------------------------------------------------------------------------------------------
  out = list(seed=seed, scenario=scenario, n_samp=n_samp,
             ATE_bias_null=ATE_bias_null, ATE_bias_RLR=ATE_bias_RLR, ATE_bias_BART0=ATE_bias_BART0,
             ATE_bias_BART1=ATE_bias_BART1, ATE_bias_CF=ATE_bias_CF,
             # ATE_bias_CB0=ATE_bias_CB0, ATE_bias_BCM0=ATE_bias_BCM0,
             # ATE_bias_CB1=ATE_bias_CB1, ATE_bias_BCM1=ATE_bias_BCM1,
             ATE_coverage_BART0=ATE_coverage_BART0, ATE_coverage_BART1=ATE_coverage_BART1, ATE_coverage_CF=ATE_coverage_CF,
             ATE_IL_BART0=ATE_IL_BART0, ATE_IL_BART1=ATE_IL_BART1, ATE_IL_CF=ATE_IL_CF,
             CATE_RMSE_null=CATE_RMSE_null, CATE_RMSE_RLR=CATE_RMSE_RLR,
             CATE_RMSE_BART0=CATE_RMSE_BART0, CATE_RMSE_BART1=CATE_RMSE_BART1,
             CATE_RMSE_CF=CATE_RMSE_CF,
             # CATE_RMSE_CB0=CATE_RMSE_CB0, CATE_RMSE_BCM0=CATE_RMSE_BCM0,
             # CATE_RMSE_CB1=CATE_RMSE_CB1, CATE_RMSE_BCM1=CATE_RMSE_BCM1,
             rtm_RLR=rtm_RLR, rtm_BART0=rtm_BART0, rtm_BART1=rtm_BART1, rtm_CF=rtm_CF
             # rtm_CB0=rtm_CB0, rtm_BCM0=rtm_BCM0, rtm_CB1=rtm_CB1, rtm_BCM1=rtm_BCM1
  )
  return(out)
}
