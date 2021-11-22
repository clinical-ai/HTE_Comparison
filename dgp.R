dgp = function (scenario, seed) {

  ## Load empirical data and model parameters
  dfName = c('tavr', 'cli', 'nvaf', 't2dm')[scenario]
  zName = c('high_vol', 'endo_revasc', 'apixaban', 'NPH')[scenario]
  yName = c('rdm30d', 'rdm30d', 'stroke1y', 'hypoglyc1y')[scenario]
  load(file = paste(dfName, '_outcome_model.RData', sep=''))
  df = obj[['df']]
  po = obj[['po']]
  xNames = setdiff(names(df), yName)
  rm(obj)
  
  ## Define treatment effect funciton on linear predictor scale
  t_tau = qlogis(po[,2]) - qlogis(po[,1])

  ## Specify marginal effect based on real study
  t_tau_2 = t_tau
  if(scenario==1) {
    t_tau_2 = t_tau + (log(0.87)*2-mean(t_tau))
  }
  if(scenario==3) {
    t_tau_2 = t_tau + (log(0.51)-mean(t_tau))
  }
  if(scenario==4) {
    t_tau_2 = t_tau + (log(1)-mean(t_tau))
  }
  po[,2] = plogis(qlogis(po[,1]) + t_tau_2)
  tau = po[,2] - po[,1]
  
  ## Sample factual outcome
  pi_fact = ifelse(df[,zName]==0, po[,1], po[,2])
  set.seed(seed)
  y_fact = rbinom(nrow(df), 1, pi_fact)
  
  ## Bind empirical and simulated data
  df_sim = cbind(df[,xNames], y_fact, pi_fact, tau)
  
  return(df_sim) 
}