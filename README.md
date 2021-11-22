This file documents the contents of the project ./HeterogeneousTreatmentEffects.
The following files are present.

outome_model.R: code used to model the observed binary outcome with a
deep feed-forward neural network and generate initial estimates of the
conditional outcome probabilities given control and treatment (Step 1-2 of the
data generating process in Section 3.2 of the manuscript).

dgp.R: user-written function to generate a dataset (i.e. simulate treatment
effects and outcomes according to Step 2-4 of the data generating process).

hte_estimation.R: user-written function to estimate the conditional average
treatment effects in a simulated dataset and evaluate the performance of the
estimators described in Section 3.3. The R package for causal boosting and
causal MARS is not provided as it is still under development. 