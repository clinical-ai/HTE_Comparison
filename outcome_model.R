rm(list = ls())
setwd('/home')

## Load dataset
scenario = 1
dfName = c('tavr', 'cli', 'nvaf', 't2dm')[scenario]
yName = c('rdm30d', 'rdm30d', 'stroke1y', 'hypoglyc1y')[scenario]
zName = c('high_vol', 'endo_revasc', 'apixaban', 'NPH')[scenario]
data = read.csv(paste(dfName, '.csv', sep=''), header = T)
xNames = setdiff(names(data), yName)
data[,yName] <- as.factor(data[,yName])
unique(data[,yName])

## Deep learning with H2O
library(h2o)
localH2O = h2o.init(nthreads=-1)
df = as.h2o(data)
splits = h2o.splitFrame(df, c(0.8), seed=5)
train  = h2o.assign(splits[[1]], "train.hex") # 80%
valid  = h2o.assign(splits[[2]], "valid.hex") # 20%

## Grid search of optimal hyperparameters
hyper_params = list(
  activation=c("Rectifier","Tanh","Maxout","RectifierWithDropout","TanhWithDropout","MaxoutWithDropout"),
  hidden=list(rep(50,2),rep(50,3),rep(100,2),rep(100,3),c(50,100)),
  input_dropout_ratio=c(0.2,0.5),
  l2=c(0,1e-3)
)
hyper_params
grid = h2o.grid(
  algorithm="deeplearning",
  grid_id="dl_grid",
  training_frame=train,
  validation_frame=valid,
  x=xNames,
  y=yName,
  epochs=1000,
  stopping_metric="logloss",
  stopping_tolerance=1e-3,        ## stop when misclassification does not improve by >= 0.1% for 2 scoring events
  hyper_params=hyper_params
)
grid = h2o.getGrid("dl_grid",sort_by="logloss",decreasing=FALSE)
grid

## Fit best model to the training dataset (update the hyperparemeter values based on the output of 'grid')
dff.fit = h2o.deeplearning(x = xNames, y = yName, training_frame = train, validation_frame = valid,
                            hidden = c(50, 100),
                            activation = 'RectifierWithDropout',
                            input_dropout_ratio = 0.2,
                            l2 = 0.00,
                            epochs = 10000,
                            score_training_samples = 0,
                            stopping_metric = "logloss",
                            stopping_tolerance = 1e-3,
                            nfolds = 0
                            )
h2o.performance(dff.fit, valid = T)

## Get estimated probabilities for the potential outcomes
df_Z0 = df_Z1 = df
df_Z0[,zName] = 0
df_Z1[,zName] = 1
po = array(0, c(nrow(data), 2))
colnames(po) = c('Z0', 'Z1')
po[,1] = as.vector(predict(dff.fit, df_Z0)[,3])
po[,2] = as.vector(predict(dff.fit, df_Z1)[,3])
h2o.shutdown(prompt=FALSE)

## Export results
fName = paste(dfName, '_outcome_model.RData', sep='')
obj = list(df=data, po=po)
save(obj, file = fName)

