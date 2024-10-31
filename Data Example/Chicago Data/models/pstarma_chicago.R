library("PNAR")
library("mgcv")
# library("xtable")
library("glmSTARMA")
library("tictoc")
data("crime")
data("crime_W")

population <- read.csv("../data/pop.csv")$x
unemployed <- read.csv("../data/unemp.csv")$x
wealth <- read.csv("../data/wealth.csv")$x
young_males <- read.csv("../data/ym.csv")$x
temperature <- read.csv2("../data/temperature.csv")$Temperatur

W2 <- 1 * (crime_W > 0)
W2 <- W2 %*% W2
diag(W2) <- 0
W2[crime_W > 0] <- 0
W2 <- W2 / colSums(W2)


# EDA
table(rowSums(crime_W > 0))
sum(crime)
sum(crime) / 552
sum(crime) / 72
max(crime)
which(crime == max(crime), arr.ind = TRUE)

min(population)
max(population)

min(unemployed)
max(unemployed)

cor(cbind(population, unemployed, wealth, young_males))
cor(cbind(population, unemployed * population, wealth, young_males))
cor(wealth, population)


#pdf("../plots/crime_temperature.pdf", width = 11, height = 9)
postscript("../plots/crime_temperature.ps", width = 11, height = 9)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(apply(crime, 1, sum) ~ temperature, ylab = "Monthly number of burglaries", xlab = "Temperature (°F)",
     cex = 2.5, cex.lab=2.5, cex.axis=2.5)
crime_temp <- lm(apply(crime, 1, sum) ~ temperature)
abline(crime_temp, lwd = 1.5)
dev.off()

crime_sum <- ts(apply(crime, 1, sum), start = c(2010, 1), end = c(2015, 12), frequency = 12)
time <- 1:72


### Model fitting without covariates:

orders <- expand.grid(obs = c(1, 2, 3), means = c(0, 1, 2), link = c("identity", "log"), stringsAsFactors = FALSE)
orders <- subset(orders, obs >= means)
models <- vector(mode = "list", length = nrow(orders))
fitting_times <- numeric(length(models))
for(i in seq(length(models))){
  if(orders$means[i] > 0){
    tic()
    models[[i]] <- glmstarma(t(crime), model = list(past_obs = rep(2, orders$obs[i]), 
                                                    past_mean = rep(1, orders$means[i])),
                             wlist = list(diag(552), crime_W, W2), family = vpoisson(orders$link[i]),
                             control = list(method = "nloptr", maxit = 1000))
    time <- toc()
    fitting_times[i] <- time$toc - time$tic
  } else {
    tic()
    models[[i]] <- glmstarma(t(crime), model = list(past_obs = rep(2, orders$obs[i])),
                             wlist = list(diag(552), crime_W, W2), family = vpoisson(orders$link[i]),
                             control = list(method = "nloptr", maxit = 1000))
    time <- toc()
    fitting_times[i] <- time$toc - time$tic
  }
}

orders$qic <- sapply(models, QIC)
orders$time <- fitting_times
orders$mse <- sapply(models, function(x) mean(residuals(x)^2))

orders$deviance_explained <- sapply(models, function(x) {
  ll1 <- sum(dpois(x$ts[, -seq(x$max_time_lag)], x$fitted.values[, -seq(x$max_time_lag)], log = TRUE) - dpois(x$ts[, -seq(x$max_time_lag)], x$ts[, -seq(x$max_time_lag)], log = TRUE))
  ll2 <- sum(dpois(x$ts[, -seq(x$max_time_lag)], mean(x$ts), log = TRUE) - dpois(x$ts[, -seq(x$max_time_lag)], x$ts[, -seq(x$max_time_lag)], log = TRUE))
  return(1 - ll1 / ll2)
})

orders$variance_explained <- sapply(models, function(x) {
  ll1 <- mean(residuals(x)^2)
  ll2 <- var(c(x$ts))
  return(1 - ll1 / ll2)
})


to_consider <- orders$obs == orders$means
fitting_times <- fitting_times[to_consider]


### Modelfitting with covariates
covariates_log <- list(population = TimeConstant(log(population / 1000)), 
                   unemployed = TimeConstant(unemployed),
                   wealth = TimeConstant(wealth),
                   youngmales = TimeConstant(young_males / population),
                   temperature = SpatialConstant(temperature), 
                   trend = SpatialConstant(72:1)
)

covariates_linear <- list(population = TimeConstant(population / 1000), 
                          unemployed = TimeConstant(unemployed * population / 100),
                          wealth = TimeConstant(log1p(exp(wealth))),
                          youngmales = TimeConstant(young_males / 100),
                          temperature = SpatialConstant(temperature), 
                          trend = SpatialConstant(72:1)
)

covariates_log_train <- list(population = TimeConstant(log(population / 1000)), 
                             unemployed = TimeConstant(unemployed),
                             wealth = TimeConstant(wealth),
                             youngmales = TimeConstant(young_males / population),
                             temperature = SpatialConstant(temperature[1:60]), 
                             trend = SpatialConstant(72:13)
)

covariates_log_test <- list(population = TimeConstant(log(population / 1000)), 
                             unemployed = TimeConstant(unemployed),
                             wealth = TimeConstant(wealth),
                             youngmales = TimeConstant(young_males / population),
                             temperature = SpatialConstant(temperature[61:72]), 
                             trend = SpatialConstant(12:1)
)


covariates_linear_train <- list(population = TimeConstant(population / 1000), 
                          unemployed = TimeConstant(unemployed * population / 100),
                          wealth = TimeConstant(log1p(exp(wealth))),
                          youngmales = TimeConstant(young_males / 100),
                          temperature = SpatialConstant(temperature[1:60]), 
                          trend = SpatialConstant(72:13)
)

covariates_linear_test <- list(population = TimeConstant(population / 1000), 
                          unemployed = TimeConstant(unemployed * population / 100),
                          wealth = TimeConstant(log1p(exp(wealth))),
                          youngmales = TimeConstant(young_males / 100),
                          temperature = SpatialConstant(temperature[61:72]), 
                          trend = SpatialConstant(12:1)
)



# orders_covariates <- expand.grid(obs = c(1, 2, 3), means = c(0, 1, 2), link = c("identity", "log"), stringsAsFactors = FALSE)
models_covariates <- vector(mode = "list", length = nrow(orders))
fitting_times_cov <- numeric(length(models))
for(i in seq(length(models))){
  covariates <- covariates_linear
  meth <- "optim"
  if(orders$link[i] == "log"){
    covariates <- covariates_log
    meth <- "nloptr"
  }
  if(orders$means[i] > 0){
    tic()
    models_covariates[[i]] <- try(glmstarma(t(crime), model = list(past_obs = rep(2, orders$obs[i]), 
                                                    past_mean = rep(1, orders$means[i])),
                             wlist = list(diag(552), crime_W, W2), covariates = covariates,
                             family = vpoisson(orders$link[i]),
                             control = list(method = meth, maxit = 1000)))
    time <- toc()
    fitting_times_cov[i] <- time$toc - time$tic
  } else {
    tic()
    models_covariates[[i]] <- glmstarma(t(crime), model = list(past_obs = rep(2, orders$obs[i])),
                                        covariates = covariates, wlist = list(diag(552), crime_W, W2), family = vpoisson(orders$link[i]),
                                        control = list(method = "nloptr", maxit = 1000))
    time <- toc()
    fitting_times_cov[i] <- time$toc - time$tic
  }
}

models_covariates <- models_covariates[orders$obs == orders$means]
fitting_times_cov <- fitting_times_cov[orders$obs == orders$means]
orders <- orders[orders$obs == orders$means,]

orders$qic_covariates <- sapply(models_covariates, QIC)
orders$time_covariates <- fitting_times_cov
orders$mse_covariates <- sapply(models_covariates, function(x) mean(residuals(x)^2))

orders$deviance_explained_covariates <- sapply(models_covariates, function(x) {
  ll1 <- sum(dpois(x$ts[, -seq(x$max_time_lag)], x$fitted.values[, -seq(x$max_time_lag)], log = TRUE) - dpois(x$ts[, -seq(x$max_time_lag)], x$ts[, -seq(x$max_time_lag)], log = TRUE))
  ll2 <- sum(dpois(x$ts[, -seq(x$max_time_lag)], mean(x$ts), log = TRUE) - dpois(x$ts[, -seq(x$max_time_lag)], x$ts[, -seq(x$max_time_lag)], log = TRUE))
  return(1 - ll1 / ll2)
})

orders$variance_explained_covariates <- sapply(models_covariates, function(x) {
  ll1 <- mean(residuals(x)^2)
  ll2 <- var(c(x$ts))
  return(1 - ll1 / ll2)
})




ft <- c(fitting_times, fitting_times_cov)

models_considered <- models[to_consider]
models_considered <- c(models_considered, models_covariates)
models_considered <- models_considered[c(1, 5, 2, 6, 3, 7, 4, 8)]
ft <- ft[c(1, 5, 2, 6, 3, 7, 4, 8)]
mses <- c(orders$mse, orders$mse_covariates)
mses <- mses[c(1, 5, 2, 6, 3, 7, 4, 8)]



extract_coefs <- function(models){
  params <- matrix(NA, nrow = 20, length(models))
  for(i in seq(length(models))){
    params[1, i] <- models[[i]]$coefficients_list$intercept
    if(!is.null(models[[i]]$coefficients_list$past_mean)){
      params[2:5, i] <- c(models[[i]]$coefficients_list$past_mean, rep(NA, 4 - length(models[[i]]$coefficients_list$past_mean)))
    }
    params[6:11, i] <- c(models[[i]]$coefficients_list$past_obs, rep(NA, 6 - length(models[[i]]$coefficients_list$past_obs)))
    if(!is.null(models[[i]]$coefficients_list$covariates)){
      params[12:17, i] <- c(models[[i]]$coefficients_list$covariates) 
    }
    params[18, i] <- models[[i]]$qic / 1000
    params[19, i] <- ft[i]
    params[20, i] <- mses[i]
  }
  return(params)
}

to_save <- extract_coefs(models_considered)
# saveRDS(to_save, "../results/pstarma/params_est_pstarma.rds")
# xtable(extract_coefs(models_considered), digits = 4)

# summaries <- lapply(models_considered, summary)


####

# pdf("../plots/crime_accu.pdf", width = 11, height = 9)
postscript("../plots/crime_accu.ps", width = 11, height = 9)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(rowSums(crime), ylab = "Accumulated Buglaries", xlab = "Month",
     cex = 2.5, cex.lab=2.5, cex.axis=2.5)
points(c(NA, NA, colSums(fitted.values(models_considered[[4]]))[-1]), type = "l", col = "lightblue", lwd = 1.5)
points(c(NA, NA, colSums(fitted.values(models_considered[[8]]))[-1]), type = "l", col = "orange", lwd = 1.5)
legend("topright", legend = c("linear", "log-linear"), lty = 1, col = c("lightblue", "orange"),
       cex = 2.5, lwd = 1.5)
dev.off()
##



# Train-Test
qic_cov <- numeric(4)
aic_cov <- numeric(4)
bic_cov <- numeric(4)
fitting_time_cov <- numeric(4)
mse_cov <- numeric(4)
mse_train_cov <- numeric(4)
mse_test_cov <- numeric(4)
parameters_cov <- vector("list", 4)

deviance_explained_train <- numeric(4)
deviance_explained_covariates_train <- numeric(4)
deviance_explained_test <- numeric(4)
deviance_explained_covariates_test <- numeric(4)

variance_explained_train <- numeric(4)
variance_explained_covariates_train <- numeric(4)
variance_explained_test <- numeric(4)
variance_explained_covariates_test <- numeric(4)


for(i in seq(4)){
  covariates <- covariates_linear_train
  meth <- "optim"
  if(orders$link[i] == "log"){
    covariates <- covariates_log_train
    meth <- "nloptr"
  }
  if(orders$means[i] > 0){
    tic()
    fit_train <- try(glmstarma(t(crime)[, 1:60], model = list(past_obs = rep(2, orders$obs[i]), 
                                                                   past_mean = rep(1, orders$means[i])),
                                            wlist = list(diag(552), crime_W, W2), covariates = covariates,
                                            family = vpoisson(orders$link[i]),
                                            control = list(method = meth, maxit = 1000)))
    time <- toc()
    fitting_times_cov[i] <- time$toc - time$tic
    predictions <- predict(fit_train, newobs = t(crime)[, 61:72], newxreg = covariates_log_test, n.ahead = 12)
  } else {
    tic()
    fit_train <- try(glmstarma(t(crime)[, 1:60], model = list(past_obs = rep(2, orders$obs[i])),
                               wlist = list(diag(552), crime_W, W2), covariates = covariates,
                               family = vpoisson(orders$link[i]),
                               control = list(method = "nloptr", maxit = 1000)))
    time <- toc()
    fitting_times_cov[i] <- time$toc - time$tic
    predictions <- predict(fit_train, newobs = t(crime)[, 61:72], newxreg = covariates_log_test, n.ahead = 12)
    time <- toc()
    fitting_times_cov[i] <- time$toc - time$tic
  }
  
  ll1 <- sum(dpois(fit_train$ts[, -seq(fit_train$max_time_lag)], fit_train$fitted.values[, -seq(fit_train$max_time_lag)], log = TRUE) - dpois(fit_train$ts[, -seq(fit_train$max_time_lag)], fit_train$ts[, -seq(fit_train$max_time_lag)], log = TRUE))
  ll2 <- sum(dpois(fit_train$ts[, -seq(fit_train$max_time_lag)], mean(fit_train$ts), log = TRUE) - dpois(fit_train$ts[, -seq(fit_train$max_time_lag)], fit_train$ts[, -seq(fit_train$max_time_lag)], log = TRUE))
  deviance_explained_covariates_train[i] <- 1 - ll1 / ll2 
  
  ll1 <- mean(residuals(fit_train)^2)
  ll2 <- var(c(fit_train$ts))
  variance_explained_covariates_train[i] <- 1 - ll1 / ll2
  
  mse_train_cov[i] <- mean(residuals(fit_train)^2)
  mse_test_cov[i] <- mean((predictions - t(crime)[, 61:72])^2)
  
  
  ll1 <- sum(dpois(t(crime)[, 61:72], predictions, log = TRUE) - dpois(t(crime)[, 61:72], t(crime)[, 61:72], log = TRUE))
  ll2 <- sum(dpois(t(crime)[, 61:72], mean(t(crime)[, 61:72]), log = TRUE) - dpois(t(crime)[, 61:72], t(crime)[, 61:72], log = TRUE))
  deviance_explained_covariates_test[i] <- 1 - ll1 / ll2 
  
  ll1 <- mean((predictions - t(crime)[, 61:72])^2)
  ll2 <- var(c(t(crime)[, 61:72]))
  variance_explained_covariates_test[i] <- 1 - ll1 / ll2
}
orders$mse_train_cov <- mse_train_cov
orders$mse_test_cov <- mse_test_cov
orders$deviance_explained_covariates_train <- deviance_explained_covariates_train
orders$deviance_explained_covariates_test <- deviance_explained_covariates_test
orders$variance_explained_covariates_train <- variance_explained_covariates_train
orders$variance_explained_covariates_test <- variance_explained_covariates_test



mse_train_ar <- numeric(4)
mse_test_ar <- numeric(4)




for(i in seq(4)){
  meth <- "optim"
  if(orders$link[i] == "log"){
    meth <- "nloptr"
  }
  if(orders$means[i] > 0){
    tic()
    fit_train <- try(glmstarma(t(crime)[, 1:60], model = list(past_obs = rep(2, orders$obs[i]), 
                                                              past_mean = rep(1, orders$means[i])),
                               wlist = list(diag(552), crime_W, W2),
                               family = vpoisson(orders$link[i]),
                               control = list(method = meth, maxit = 1000)))
    time <- toc()
    fitting_times_cov[i] <- time$toc - time$tic
    predictions <- predict(fit_train, newobs = t(crime)[, 61:72], n.ahead = 12)
  } else {
    tic()
    fit_train <- try(glmstarma(t(crime)[, 1:60], model = list(past_obs = rep(2, orders$obs[i])),
                               wlist = list(diag(552), crime_W, W2),
                               family = vpoisson(orders$link[i]),
                               control = list(method = "nloptr", maxit = 1000)))
    time <- toc()
    fitting_times_cov[i] <- time$toc - time$tic
    predictions <- predict(fit_train, newobs = t(crime)[, 61:72], n.ahead = 12)
    time <- toc()
    fitting_times_cov[i] <- time$toc - time$tic
  }
  
  ll1 <- sum(dpois(fit_train$ts[, -seq(fit_train$max_time_lag)], fit_train$fitted.values[, -seq(fit_train$max_time_lag)], log = TRUE) - dpois(fit_train$ts[, -seq(fit_train$max_time_lag)], fit_train$ts[, -seq(fit_train$max_time_lag)], log = TRUE))
  ll2 <- sum(dpois(fit_train$ts[, -seq(fit_train$max_time_lag)], mean(fit_train$ts), log = TRUE) - dpois(fit_train$ts[, -seq(fit_train$max_time_lag)], fit_train$ts[, -seq(fit_train$max_time_lag)], log = TRUE))
  deviance_explained_train[i] <- 1 - ll1 / ll2 
  
  ll1 <- mean(residuals(fit_train)^2)
  ll2 <- var(c(fit_train$ts))
  variance_explained_train[i] <- 1 - ll1 / ll2
  
  mse_train_ar[i] <- mean(residuals(fit_train)^2)
  mse_test_ar[i] <- mean((predictions - t(crime)[, 61:72])^2)
  
  ll1 <- sum(dpois(t(crime)[, 61:72], predictions, log = TRUE) - dpois(t(crime)[, 61:72], t(crime)[, 61:72], log = TRUE))
  ll2 <- sum(dpois(t(crime)[, 61:72], mean(t(crime)[, 61:72]), log = TRUE) - dpois(t(crime)[, 61:72], t(crime)[, 61:72], log = TRUE))
  deviance_explained_test[i] <- 1 - ll1 / ll2 
  
  ll1 <- mean((predictions - t(crime)[, 61:72])^2)
  ll2 <- var(c(t(crime)[, 61:72]))
  variance_explained_test[i] <- 1 - ll1 / ll2
}


orders$mse_train_ar <- mse_train_ar
orders$mse_test_ar <- mse_test_ar
orders$deviance_explained_train <- deviance_explained_train
orders$deviance_explained_test <- deviance_explained_test
orders$variance_explained_train <- variance_explained_train
orders$variance_explained_test <- variance_explained_test



# saveRDS(orders, "../results/pstarma/results.rds")

