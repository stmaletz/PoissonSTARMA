
x <- list.files("models_new/")


training <- list.files("models_new/", full.names = TRUE)
training <- training[!grepl("predict", training)]
test <- list.files("models_new/", full.names = TRUE, pattern = "predict")

fitting_times <- numeric(length(training))
mse <- numeric(length(training))
mse_training <- numeric(length(test))
mse_test <- numeric(length(test))

deviance_explained <- numeric(length(training))
deviance_explained_training <- numeric(length(training))
deviance_explained_test <- numeric(length(test))

variance_explained <- numeric(length(training))
variance_explained_training <- numeric(length(training))
variance_explained_test <- numeric(length(test))




for(i in seq_along(training)){
  model <- readRDS(training[i])
  fitting_times[i] <- model$fitting_time
  mse[i] <- mean(residuals(model$model, type = "response")^2)
  
  deviance_train_fit <- -2 * sum(dpois(model$model$y, model$model$fitted.values, log = TRUE) - dpois(model$model$y, model$model$y, log = TRUE))
  deviance_train <- -2 * sum(dpois(model$model$y, mean(model$model$y), log = TRUE) - dpois(model$model$y, model$model$y, log = TRUE))
  deviance_explained[i] <- 1 - deviance_train_fit / deviance_train
  
  variance_explained[i] <- 1 - mean(residuals(model$model, type = "response")^2) / var(c(model$model$y))
  
  
  
}

for(i in seq_along(test)){
  model <- readRDS(test[i])
  mse_training[i] <- mean(residuals(model$model, type = "response")^2)
  mse_test[i] <- mean((model$prediction - model$truth)^2)
  
  deviance_train_fit <- -2 * sum(dpois(model$model$y, model$model$fitted.values, log = TRUE) - dpois(model$model$y, model$model$y, log = TRUE))
  deviance_train <- -2 * sum(dpois(model$model$y, mean(model$model$y), log = TRUE) - dpois(model$model$y, model$model$y, log = TRUE))
  deviance_explained_training[i] <- 1 - deviance_train_fit / deviance_train
  
  
  deviance_test_fit <- -2* sum(dpois(model$truth, model$prediction, log = TRUE) - dpois(model$truth, model$truth, log = TRUE))
  deviance_test <- -2 * sum(dpois(model$truth, mean(model$truth), log = TRUE) - dpois(model$truth, model$truth, log = TRUE))
  deviance_explained_test[i] <- 1 - deviance_test_fit / deviance_test
  
  variance_explained_training[i] <- 1 - mean(residuals(model$model, type = "response")^2) / var(c(model$model$y))
  variance_explained_test[i] <- 1 - mean((model$prediction - model$truth)^2) / var(c(model$truth))
}

test <- gsub("predict_", "", test)
all(test == training)

infos <- strsplit(training, "_|.rds")
infos <- data.frame(trend = as.logical(sapply(infos, "[[", 5)),
                  space = as.numeric(sapply(infos, "[[", 7)),
                  time = as.numeric(sapply(infos, "[[", 9)))



results <- data.frame(trend = infos$trend, space = infos$space, time = infos$time,
  fitting_time = fitting_times, mse = mse, mse_training = mse_training, mse_test = mse_test,
  deviance_explained = deviance_explained, deviance_train = deviance_explained_training,
  deviance_test = deviance_explained_test,
  var_explained = variance_explained, var_train = variance_explained_training,
  var_test = variance_explained_test)

results[order(results$space, results$time, results$trend),]
results[order(results$mse_test),]


res_pstarma <- readRDS("results.rds")



training <- list.files("hhh4_models/", full.names = TRUE)
training <- training[!grepl("predict", training)]
test <- list.files("hhh4_models//", full.names = TRUE, pattern = "predict")



fitting_times <- numeric(length(training))
mse <- numeric(length(training))
mse_training <- numeric(length(test))
mse_test <- numeric(length(test))

deviance_explained <- numeric(length(training))
deviance_explained_training <- numeric(length(training))
deviance_explained_test <- numeric(length(test))

variance_explained <- numeric(length(training))
variance_explained_training <- numeric(length(training))
variance_explained_test <- numeric(length(test))




for(i in seq_along(training)){
  model <- readRDS(training[i])
  fitting_times[i] <- model$fitting_time$toc - model$fitting_time$tic
  mse[i] <- mean(residuals(model$model, type = "response")^2)
  
  deviance_train_fit <- sum(dpois(model$model$stsObj@observed[model$model$control$subset,], model$model$fitted.values, log = TRUE) - dpois(model$model$stsObj@observed[model$model$control$subset,], model$model$stsObj@observed[model$model$control$subset,], log = TRUE))
  deviance_train <- sum(dpois(model$model$stsObj@observed[model$model$control$subset,], mean(model$model$stsObj@observed[model$model$control$subset,]), log = TRUE) - dpois(model$model$stsObj@observed[model$model$control$subset,], model$model$stsObj@observed[model$model$control$subset,], log = TRUE))
  deviance_explained[i] <- 1 - deviance_train_fit / deviance_train
  
  variance_explained[i] <- 1 - mean(residuals(model$model, type = "response")^2) / var(c(model$model$stsObj@observed))
  
}

for(i in seq_along(test)){
  model <- readRDS(test[i])
  mse_training[i] <- mean(residuals(model$model, type = "response")^2)
  mse_test[i] <- mean((model$prediction - model$truth)^2)
  
  deviance_train_fit <- sum(dpois(model$model$stsObj@observed[model$model$control$subset,], model$model$fitted.values, log = TRUE) - dpois(model$model$stsObj@observed[model$model$control$subset,], model$model$stsObj@observed[model$model$control$subset,], log = TRUE))
  deviance_train <- sum(dpois(model$model$stsObj@observed[model$model$control$subset,], mean(model$model$stsObj@observed[model$model$control$subset,]), log = TRUE) - dpois(model$model$stsObj@observed[model$model$control$subset,], model$model$stsObj@observed[model$model$control$subset,], log = TRUE))
  deviance_explained_training[i] <- 1 - deviance_train_fit / deviance_train
  
  
  deviance_test_fit <- sum(dpois(model$truth, model$prediction, log = TRUE) - dpois(model$truth, model$truth, log = TRUE))
  deviance_test <- sum(dpois(model$truth, mean(model$truth), log = TRUE) - dpois(model$truth, model$truth, log = TRUE))
  deviance_explained_test[i] <- 1 - deviance_test_fit / deviance_test
  
  variance_explained_training[i] <- 1 - mean(residuals(model$model, type = "response")^2) / var(c(model$model$stsObj@observed))
  variance_explained_test[i] <- 1 - mean((model$prediction - model$truth)^2) / var(c(model$truth))
}

training <- gsub("///", "//", training)
test <- gsub("predict_", "", test)
test <- gsub("///", "//", test)
all(test == training)
results_hhh4 <- data.frame(fitting_time = fitting_times, mse = mse, mse_training = mse_training, mse_test = mse_test,
                           deviance_explained = deviance_explained, deviance_train = deviance_explained_training,
                           deviance_test = deviance_explained_test,
                           var_explained = variance_explained, var_train = variance_explained_training,
                           var_test = variance_explained_test)

