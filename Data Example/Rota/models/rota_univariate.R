library("glmSTARMA")
library("tscount")
library("rjson")
library("spdep")
library("tictoc")

## Rota Fit
rota <- read.csv2("../data/rotavirus.csv")
rota <- t(as.matrix(rota[-1]))
rota <- rota[-413,]
rownames(rota) <- gsub("X", "", rownames(rota))
rota <- rota[match(shape$RKI_ID, row.names(rota)),]

covariates <- cbind(season_cos = cos(2 * pi / 52 * seq(ncol(rota))), 
                    season_sin = sin(2 * pi / 52 * seq(ncol(rota))),
                    vaccine = 1 * (seq(ncol(rota)) >= 654))

covariates_train <- covariates[seq(835),]
covariates_test <- covariates[-seq(835),]



# Fit models

qic_ar <- numeric(8)
aic_ar <- numeric(8)
bic_ar <- numeric(8)
fitting_time_ar <- numeric(8)
mse_ar <- numeric(8)
mse_train_ar <- numeric(8)
mse_test_ar <- numeric(8)

deviance_ar <- numeric(8)
deviance_train_ar <- numeric(8)
deviance_test_ar <- numeric(8)

parameters_ar <- vector("list", 8)
uni_fit <- vector(mode = "list", nrow(rota))
uni_fit_train <- vector(mode = "list", nrow(rota))
test_pred <- matrix(NA, ncol = 68, nrow = nrow(rota))

for(i in 1:8){
  cat("Index", i, "\n")
  tic()
  for(j in seq(nrow(rota))){
    uni_fit[[j]] <- tsglm(rota[j,], model = list(past_obs = 1:i), xreg = covariates, link = "log")
  }
  time <- toc()
  #qic_ar[i] <- fit$qic
  #aic_ar[i] <- fit$aic
  #bic_ar[i] <- fit$bic
  fitting_time_ar[i] <- time$toc - time$tic
  
  
  mse_ar[i] <- mean(sapply(uni_fit, function(x) mean(residuals(x, type = "response")^2)))
  
  ll1 <- sum(sapply(uni_fit, function(x){
    sum(dpois(x$ts[-seq(i)], x$fitted.values[-seq(i)], log = TRUE) - dpois(x$ts[-seq(i)], x$ts[-seq(i)], log = TRUE))
  }))
  ll2 <- sum(dpois(rota[, -seq(i)], mean(rota), log = TRUE) - dpois(rota[, -seq(i)], rota[, -seq(i)], log = TRUE))

  deviance_ar[i] <- 1 - ll1 / ll2
  
  # parameters_ar[[i]] <- summary(fit)
  
  for(j in seq(nrow(rota))){
    uni_fit_train[[j]] <- tsglm(rota[j, 1:835], model = list(past_obs = 1:i), xreg = covariates_train, link = "log")
    test_pred[j,] <- predict(uni_fit_train[[j]], newobs = rota[j, 836:903], n.ahead = 68, newxreg = covariates_test)$pred
  }
  
  mse_train_ar[i] <- mean(sapply(uni_fit_train, function(x) mean(residuals(x, type = "response")^2)))
  mse_test_ar[i] <- mean((test_pred - rota[, 836:903])^2)
  
  ll1 <- sum(sapply(uni_fit_train, function(x){
    sum(dpois(x$ts[-seq(i)], x$fitted.values[-seq(i)], log = TRUE) - dpois(x$ts[-seq(i)], x$ts[-seq(i)], log = TRUE))
  }))
  ll2 <- sum(dpois(rota[, -c(seq(i), 836:903)], mean(rota[, -c(836:903)]), log = TRUE) - dpois(rota[, -c(seq(i), 836:903)], rota[, -c(seq(i), 836:903)], log = TRUE))
  deviance_train_ar[i] <- 1 - ll1 / ll2
  
  ll1 <- sum(dpois(rota[, 836:903], test_pred, log = TRUE) - dpois(rota[, 836:903], rota[, 836:903], log = TRUE))
  ll2 <- sum(dpois(rota[, 836:903], mean(rota[, 836:903]), log = TRUE) - dpois(rota[, 836:903], rota[, 836:903], log = TRUE))
  deviance_test_ar[i] <- 1 - ll1 / ll2
}


ar_models <- list(
  #qic = qic_ar,
  #aic = aic_ar,
  #bic = bic_ar,
  fitting_time = fitting_time_ar,
  mse = mse_ar,
  mse_train = mse_train_ar,
  mse_test = mse_test_ar,
  deviance = deviance_ar,
  deviance_train = deviance_train_ar,
  deviance_test = deviance_test_ar
  # parameters = parameters_ar
)

saveRDS(ar_models, "../results/univariate/rota_results_uni.rds")