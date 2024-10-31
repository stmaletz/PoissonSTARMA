library("glmSTARMA")
library("rjson")
library("spdep")
library("tictoc")

## Erste beiden Ziffern geben Bundesland an.
## 05 ist NRW


shape <- read_sf("../data/germany_county_shapes.json")
nb <- poly2nb(shape, row.names = shape$RKI_ID)
W1 <- nb2mat(nb, style="W")
row.names(W1) <- colnames(W1) <- shape$RKI_ID

W2 <- 1 * (W1 > 0)
W2 <- W2 %*% W2
diag(W2) <- 0
W2[W1 > 0] <- 0
W2 <- W2 / colSums(W2)



## Read in features:
ddr_feature <- read.csv("../data/region_features.csv")
ddr_feature$RKI_ID <- ifelse(nchar(as.character(ddr_feature$RKI_ID)) == 4, paste0("0", ddr_feature$RKI_ID), as.character(ddr_feature$RKI_ID))
ddr_feature <- ddr_feature[match(shape$RKI_ID, ddr_feature$RKI_ID),]
rownames(ddr_feature) <- ddr_feature$RKI_ID
ddr_feature <- ddr_feature$is_east

population <- read.csv("../data/germany_population_data.csv")
population <- aggregate(population ~ county + year, data = population, FUN = sum)
population <- reshape(population, timevar = "year", idvar = "county", direction = "wide")
population <- population[match(shape$RKI_NameDE, population$county),]
population <- population$population.2001



## Rota Fit
rota <- read.csv2("../data/rotavirus.csv")
rota <- t(as.matrix(rota[-1]))
rota <- rota[-413,]
rownames(rota) <- gsub("X", "", rownames(rota))
rota <- rota[match(shape$RKI_ID, row.names(rota)),]


covariates_log <- list()
covariates_log$is_east <- TimeConstant(ddr_feature)
covariates_log$population <- TimeConstant(log(population / 100000))
covariates_log$season_cos <- SpatialConstant(cos(2 * pi / 52 * seq(ncol(rota))))
covariates_log$season_sin <- SpatialConstant(sin(2 * pi / 52 * seq(ncol(rota))))
covariates_log$vaccine <- SpatialConstant(1 * (seq(ncol(rota)) >= 654))


covariates_log_train <- list()
covariates_log_train$is_east <- TimeConstant(ddr_feature)
covariates_log_train$population <- TimeConstant(log(population / 100000))
covariates_log_train$season_cos <- SpatialConstant(cos(2 * pi / 52 * seq(835)))
covariates_log_train$season_sin <- SpatialConstant(sin(2 * pi / 52 * seq(835)))
covariates_log_train$vaccine <- SpatialConstant(1 * (seq(835) >= 654))


covariates_log_test <- list()
covariates_log_test$is_east <- TimeConstant(ddr_feature)
covariates_log_test$population <- TimeConstant(log(population / 100000))
covariates_log_test$season_cos <- SpatialConstant(cos(2 * pi / 52 * (836:903)))
covariates_log_test$season_sin <- SpatialConstant(sin(2 * pi / 52 * (836:903)))
covariates_log_test$vaccine <- SpatialConstant(1 * (836:903 >= 654))




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

for(i in 1:8){
  tic()
  fit <- glmstarma(rota, model = list(past_obs = rep(2, i)), wlist = list(diag(412), W1, W2),
                   covariates = covariates_log, 
                   family = vpoisson("log"), control = list(method = "nloptr", maxit = 10000))
  time <- toc()
  qic_ar[i] <- fit$qic
  aic_ar[i] <- fit$aic
  bic_ar[i] <- fit$bic
  fitting_time_ar[i] <- time$toc - time$tic
  mse_ar[i] <- mean(residuals(fit)^2)
  
  ll1 <- sum(dpois(fit$ts[, -seq(fit$max_time_lag)], fit$fitted.values[, -seq(fit$max_time_lag)], log = TRUE) - dpois(fit$ts[, -seq(fit$max_time_lag)], fit$ts[, -seq(fit$max_time_lag)], log = TRUE))
  ll2 <- sum(dpois(fit$ts[, -seq(fit$max_time_lag)], mean(fit$ts), log = TRUE) - dpois(fit$ts[, -seq(fit$max_time_lag)], fit$ts[, -seq(fit$max_time_lag)], log = TRUE))
  deviance_ar[i] <- 1 - ll1 / ll2
  
  parameters_ar[[i]] <- summary(fit)
  
  fit_train <- glmstarma(rota[,1:835], model = list(past_obs = rep(2, i)), wlist = list(diag(412), W1, W2),
                         covariates = covariates_log_train, 
                         family = vpoisson("log"), control = list(method = "nloptr", maxit = 10000))
  predictions <- predict(fit_train, newobs = rota[, 836:903], newxreg = covariates_log_test, n.ahead = 68)
  
  mse_train_ar[i] <- mean(residuals(fit_train)^2)
  mse_test_ar[i] <- mean((predictions - rota[, 836:903])^2)
  
  ll1 <- sum(dpois(fit_train$ts[, -seq(fit_train$max_time_lag)], fit_train$fitted.values[, -seq(fit_train$max_time_lag)], log = TRUE) - dpois(fit_train$ts[, -seq(fit_train$max_time_lag)], fit_train$ts[, -seq(fit_train$max_time_lag)], log = TRUE))
  ll2 <- sum(dpois(fit_train$ts[, -seq(fit_train$max_time_lag)], mean(fit_train$ts), log = TRUE) - dpois(fit_train$ts[, -seq(fit_train$max_time_lag)], fit_train$ts[, -seq(fit_train$max_time_lag)], log = TRUE))
  deviance_train_ar[i] <- 1 - ll1 / ll2
  
  ll1 <- sum(dpois(rota[, 836:903], predictions, log = TRUE) - dpois(rota[, 836:903], rota[, 836:903], log = TRUE))
  ll2 <- sum(dpois(rota[, 836:903], mean(rota[, 836:903]), log = TRUE) - dpois(rota[, 836:903], rota[, 836:903], log = TRUE))
  deviance_test_ar[i] <- 1 - ll1 / ll2
}


ar_models <- list(
  qic = qic_ar,
  aic = aic_ar,
  bic = bic_ar,
  fitting_time = fitting_time_ar,
  mse = mse_ar,
  mse_train = mse_train_ar,
  mse_test = mse_test_ar,
  deviance = deviance_ar,
  deviance_train = deviance_train_ar,
  deviance_test = deviance_test_ar,
  parameters = parameters_ar
)

saveRDS(ar_models, "../results/pstarma/rota_results_ar.rds")




qic_arma <- numeric(8)
aic_arma <- numeric(8)
bic_arma <- numeric(8)
fitting_time_arma <- numeric(8)
mse_arma <- numeric(8)
mse_train_arma <- numeric(8)
mse_test_arma <- numeric(8)
parameters_arma <- vector("list", 8)

deviance_arma <- numeric(8)
deviance_train_arma <- numeric(8)
deviance_test_arma <- numeric(8)

for(i in 1:8){
  tic()
  fit <- glmstarma(rota, model = list(past_obs = rep(2, i), past_mean = 1), wlist = list(diag(412), W1, W2),
                   covariates = covariates_log,
                   family = vpoisson("log"), control = list(method = "nloptr", maxit = 10000))
  time <- toc()
  qic_arma[i] <- fit$qic
  aic_arma[i] <- fit$aic
  bic_arma[i] <- fit$bic
  fitting_time_arma[i] <- time$toc - time$tic
  mse_arma[i] <- mean(residuals(fit)^2)
  
  ll1 <- sum(dpois(fit$ts[, -seq(fit$max_time_lag)], fit$fitted.values[, -seq(fit$max_time_lag)], log = TRUE) - dpois(fit$ts[, -seq(fit$max_time_lag)], fit$ts[, -seq(fit$max_time_lag)], log = TRUE))
  ll2 <- sum(dpois(fit$ts[, -seq(fit$max_time_lag)], mean(fit$ts), log = TRUE) - dpois(fit$ts[, -seq(fit$max_time_lag)], fit$ts[, -seq(fit$max_time_lag)], log = TRUE))
  deviance_arma[i] <- 1 - ll1 / ll2
  
  parameters_arma[[i]] <- summary(fit)
  
  fit_train <- glmstarma(rota[,1:835], model = list(past_obs = rep(2, i), past_mean = 2), wlist = list(diag(412), W1, W2),
                         covariates = covariates_log_train, 
                         family = vpoisson("log"), control = list(method = "nloptr", maxit = 10000))
  predictions <- predict(fit_train, newobs = rota[, 836:903], newxreg = covariates_log_test, n.ahead = 68)
  
  mse_train_arma[i] <- mean(residuals(fit_train)^2)
  mse_test_arma[i] <- mean((predictions - rota[, 836:903])^2)
  
  ll1 <- sum(dpois(fit_train$ts[, -seq(fit_train$max_time_lag)], fit_train$fitted.values[, -seq(fit_train$max_time_lag)], log = TRUE) - dpois(fit_train$ts[, -seq(fit_train$max_time_lag)], fit_train$ts[, -seq(fit_train$max_time_lag)], log = TRUE))
  ll2 <- sum(dpois(fit_train$ts[, -seq(fit_train$max_time_lag)], mean(fit_train$ts), log = TRUE) - dpois(fit_train$ts[, -seq(fit_train$max_time_lag)], fit_train$ts[, -seq(fit_train$max_time_lag)], log = TRUE))
  deviance_train_arma[i] <- 1 - ll1 / ll2
  
  ll1 <- sum(dpois(rota[, 836:903], predictions, log = TRUE) - dpois(rota[, 836:903], rota[, 836:903], log = TRUE))
  ll2 <- sum(dpois(rota[, 836:903], mean(rota[, 836:903]), log = TRUE) - dpois(rota[, 836:903], rota[, 836:903], log = TRUE))
  deviance_test_arma[i] <- 1 - ll1 / ll2
  
}

arma_models <- list(
  qic = qic_arma,
  aic = aic_arma,
  bic = bic_arma,
  fitting_time = fitting_time_arma,
  mse = mse_arma,
  mse_train = mse_train_arma,
  mse_test = mse_test_arma,
deviance = deviance_arma,
deviance_train = deviance_train_arma,
deviance_test = deviance_test_arma,
  parameters = parameters_arma
)

saveRDS(arma_models, "../results/pstarma/rota_results_arma.rds")



## Linear model for comparison

# covariates_linear <- list()
# covariates_linear$population <- TimeConstant(population / 100000)
# covariates_linear$is_east <- TimeConstant(ddr_feature)


# Fitting only possible without constraints to parameters
# model_linear_ar <- glmstarma(rota, model = list(past_obs = c(1), past_mean = 1), wlist = list(diag(412), W1),
#                          covariates = covariates_linear, 
#                          family = vpoisson("identity"), control = list(method = "optim", maxit = 1000))
