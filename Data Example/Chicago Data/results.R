hhh4_models <- list.files("results/hhh4", pattern = "hhh4", full.names = TRUE)
mgcv_models <- list.files("results/mgcv", pattern = "mgcv", full.names = TRUE)

hhh4_models <- hhh4_models[!grepl("predict", hhh4_models)]
mgcv_models <- mgcv_models[!grepl("predict", mgcv_models)]


hhh4_models_predict <- list.files("results/hhh4", pattern = "hhh4", full.names = TRUE)
mgcv_models_predict <- list.files("results/mgcv", pattern = "mgcv", full.names = TRUE)

hhh4_models_predict <- hhh4_models_predict[grepl("predict", hhh4_models_predict)]
mgcv_models_predict <- mgcv_models_predict[grepl("predict", mgcv_models_predict)]



### HHH4 Model
fitting_time_hhh4 <- numeric(length(hhh4_models))
aic_hhh4 <- numeric(length(hhh4_models))
bic_hhh4 <- numeric(length(hhh4_models))
mse_hhh4 <- numeric(length(hhh4_models))
mae_hhh4 <- numeric(length(hhh4_models))
mse_hhh4_train <- numeric(length(hhh4_models))
mae_hhh4_train <- numeric(length(hhh4_models))
mse_hhh4_test <- numeric(length(hhh4_models))
mae_hhh4_test <- numeric(length(hhh4_models))

deviance_explained <- numeric(length(hhh4_models))
deviance_explained_train <- numeric(length(hhh4_models))
deviance_explained_test <- numeric(length(hhh4_models))

for(i in seq_along(hhh4_models)){
  results <- readRDS(hhh4_models[i])
  results_predict <- readRDS(hhh4_models_predict[i])
  fitting_time_hhh4[i] <- results$fitting_time$toc - results$fitting_time$tic
  aic_hhh4[i] <- AIC(results$model)
  bic_hhh4[i] <- BIC(results$model)
  mse_hhh4[i] <- mean(residuals(results$model, type = "response")^2)
  mae_hhh4[i] <- mean(abs(residuals(results$model, type = "response")))
  
  mse_hhh4_train[i] <- mean(residuals(results_predict$model, type = "response")^2)
  mae_hhh4_train[i] <- mean(abs(residuals(results_predict$model, type = "response")))
  mse_hhh4_test[i] <- mean((results_predict$prediction - results_predict$truth)^2)
  mae_hhh4_test[i] <- mean(abs(results_predict$prediction - results_predict$truth))
  
  # TODO
  deviance_train_fit <- -2 * sum(dpois(results$model$stsObj@observed[results$model$control$subset,], results$model$fitted.values, log = TRUE) - dpois(results$model$stsObj@observed[results$model$control$subset,], results$model$stsObj@observed[results$model$control$subset,], log = TRUE))
  deviance_train <- -2 * sum(dpois(results$model$stsObj@observed[results$model$control$subset,], mean(results$model$stsObj@observed[results$model$control$subset,]), log = TRUE) - dpois(results$model$stsObj@observed[results$model$control$subset,], results$model$stsObj@observed[results$model$control$subset,], log = TRUE))
  deviance_explained[i] <- 1 - deviance_train_fit / deviance_train
  
  deviance_train_fit <- -2 * sum(dpois(results_predict$model$stsObj@observed[results_predict$model$control$subset,], results_predict$model$fitted.values, log = TRUE) - dpois(results_predict$model$stsObj@observed[results_predict$model$control$subset,], results_predict$model$stsObj@observed[results_predict$model$control$subset,], log = TRUE))
  deviance_train <- -2 * sum(dpois(results_predict$model$stsObj@observed[results_predict$model$control$subset,], mean(results_predict$model$stsObj@observed[results_predict$model$control$subset,]), log = TRUE) - dpois(results_predict$model$stsObj@observed[results_predict$model$control$subset,], results_predict$model$stsObj@observed[results_predict$model$control$subset,], log = TRUE))
  deviance_explained_train[i] <- 1 - deviance_train_fit / deviance_train
  
  deviance_test_fit <- -2 * sum(dpois(results_predict$truth, results_predict$prediction, log = TRUE) - dpois(results_predict$truth, results_predict$truth, log = TRUE))
  deviance_test <- -2 * sum(dpois(results_predict$truth, mean(results_predict$truth), log = TRUE) - dpois(results_predict$truth, results_predict$truth, log = TRUE))
  deviance_explained_test[i] <- 1 - deviance_test_fit / deviance_test
  
}

hhh4_results <- data.frame(fitting_time = fitting_time_hhh4,
           aic = aic_hhh4, bic = bic_hhh4, 
           mse = mse_hhh4, deviance_explained = deviance_explained,
           mse_train = mse_hhh4_train, deviance_explained_train = deviance_explained_train,
           mse_test = mse_hhh4_test, deviance_explained_test = deviance_explained_test)



### MGCV Model
fitting_time_mgcv <- numeric(length(mgcv_models))
aic_mgcv <- numeric(length(mgcv_models))
bic_mgcv <- numeric(length(mgcv_models))
mse_mgcv <- numeric(length(mgcv_models))
mae_mgcv <- numeric(length(mgcv_models))
mse_mgcv_train <- numeric(length(mgcv_models))
mae_mgcv_train <- numeric(length(mgcv_models))
mse_mgcv_test <- numeric(length(mgcv_models))
mae_mgcv_test <- numeric(length(mgcv_models))

deviance_explained <- numeric(length(mgcv_models))
deviance_explained_train <- numeric(length(mgcv_models))
deviance_explained_test <- numeric(length(mgcv_models))


for(i in seq_along(mgcv_models)){
  results <- readRDS(mgcv_models[i])
  results_predict <- readRDS(mgcv_models_predict[i])
  fitting_time_mgcv[i] <- results$fitting_time
  aic_mgcv[i] <- AIC(results$model)
  bic_mgcv[i] <- BIC(results$model)
  mse_mgcv[i] <- mean(residuals(results$model, type = "response")^2)
  mae_mgcv[i] <- mean(abs(residuals(results$model, type = "response")))
  
  mse_mgcv_train[i] <- mean(residuals(results_predict$model, type = "response")^2)
  mae_mgcv_train[i] <- mean(abs(residuals(results_predict$model, type = "response")))
  mse_mgcv_test[i] <- mean((results_predict$prediction - results_predict$truth)^2)
  mae_mgcv_test[i] <- mean(abs(results_predict$prediction - results_predict$truth))
  
  deviance_train_fit <- -2 * sum(dpois(results$model$y, results$model$fitted.values, log = TRUE) - dpois(results$model$y, results$model$y, log = TRUE))
  deviance_train <- -2 * sum(dpois(results$model$y, mean(results$model$y), log = TRUE) - dpois(results$model$y, results$model$y, log = TRUE))
  deviance_explained[i] <- 1 - deviance_train_fit / deviance_train
  
  deviance_train_fit <- -2 * sum(dpois(results_predict$model$y, results_predict$model$fitted.values, log = TRUE) - dpois(results_predict$model$y, results_predict$model$y, log = TRUE))
  deviance_train <- -2 * sum(dpois(results_predict$model$y, mean(results_predict$model$y), log = TRUE) - dpois(results_predict$model$y, results_predict$model$y, log = TRUE))
  deviance_explained_train[i] <- 1 - deviance_train_fit / deviance_train
  
  deviance_test_fit <- -2 * sum(dpois(results_predict$truth, results_predict$prediction, log = TRUE) - dpois(results_predict$truth, results_predict$truth, log = TRUE))
  deviance_test <- -2 * sum(dpois(results_predict$truth, mean(results_predict$truth), log = TRUE) - dpois(results_predict$truth, results_predict$truth, log = TRUE))
  deviance_explained_test[i] <- 1 - deviance_test_fit / deviance_test
}


mgcv_results <- data.frame(fitting_time = fitting_time_mgcv,
           aic = aic_mgcv, bic = bic_mgcv,
           mse = mse_mgcv, deviance_explained = deviance_explained,
           mse_train = mse_mgcv_train, deviance_explained_train = deviance_explained_train,
           mse_test = mse_mgcv_test, deviance_explained_test = deviance_explained_test)


mgcv_results[order(mgcv_results$mse_test),]

# Univariate results
uni_log <- readRDS("results/univariate/result_univariate_log.rds")
uni_linear <- readRDS("results/univariate/result_univariate_linear.rds")

# PSTARMA
# readRDS("results/pstarma/results.rds")



