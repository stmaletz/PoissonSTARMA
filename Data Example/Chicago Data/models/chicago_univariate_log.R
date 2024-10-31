library("PNAR")
library("mgcv")
library("tictoc")
library("glmSTARMA")
data("crime")
data("crime_W")
library("tscount")

population <- read.csv("../data/pop.csv")$x
unemployed <- read.csv("../data/unemp.csv")$x
wealth <- read.csv("../data/wealth.csv")$x
young_males <- read.csv("../data/ym.csv")$x
temperature <- read.csv2("../data/temperature.csv")$Temperatur

covariates <- cbind(temperature, time)
covariates_train <- covariates[1:60,]
covariates_test <- covariates[61:72,]



tic()
fit_univariate <- lapply(seq(ncol(crime)), function(i) tsglm(crime[, i], model = list(past_mean = 1, past_obs = 1), xreg = covariates,
                                                             link = "log"))
time <- toc()
fitting_time <- time$toc - time$tic

mse_full <- mean(sapply(fit_univariate, function(x) mean(residuals(x)^2)))
ll1 <- sum(sapply(fit_univariate, function(x){
  sum(dpois(x$ts[-1], x$fitted.values[-1], log = TRUE) - dpois(x$ts[-1], x$ts[-1], log = TRUE))
}))
ll2 <- sum(dpois(crime[-1, ], mean(crime), log = TRUE) - dpois(crime[-1, ], crime[-1, ], log = TRUE))
deviance_full <- 1 - ll1 / ll2



fit_univariate_train <- vector(mode = "list", ncol(crime))
test_pred <- matrix(NA, ncol = ncol(crime), nrow = 12)
for(i in seq(ncol(crime))){
  cat("Dimension", i, "\n")
  fit_univariate_train[[i]] <- tsglm(crime[1:60, i], model = list(past_mean = 1, past_obs = 1), xreg = covariates_train,
                                                                     link = "log")
  test_pred[, i] <- predict(fit_univariate_train[[i]], newobs = crime[61:72, i], n.ahead = 12, newxreg = covariates_test)$pred
}

mse_train <- mean(sapply(fit_univariate_train, function(x) mean(residuals(x, type = "response")^2)))
mse_test <- mean((test_pred - crime[61:72, ])^2)


ll1 <- sum(sapply(fit_univariate_train, function(x){
  sum(dpois(x$ts[-1], x$fitted.values[-1], log = TRUE) - dpois(x$ts[-1], x$ts[-1], log = TRUE))
}))
ll2 <- sum(dpois(crime[-c(1, 61:72),], mean(crime[-c(61:72),]), log = TRUE) - dpois(crime[-c(1, 61:72),], crime[-c(1, 61:72),], log = TRUE))
deviance_train <- 1 - ll1 / ll2

ll1 <- sum(dpois(crime[61:72,], test_pred, log = TRUE) - dpois(crime[61:72,], crime[61:72,], log = TRUE))
ll2 <- sum(dpois(crime[61:72,], mean(crime[61:72,]), log = TRUE) - dpois(crime[61:72,], crime[61:72,], log = TRUE))
deviance_test <- 1 - ll1 / ll2

result <- list(time = fitting_time,
               mse_full = mse_full, deviance_full = deviance_full,
               mse_train = mse_train, deviance_train = deviance_train,
               mse_test = mse_test, deviance_test = deviance_test)


saveRDS(result, "../results/univariate/result_univariate_log.rds")
