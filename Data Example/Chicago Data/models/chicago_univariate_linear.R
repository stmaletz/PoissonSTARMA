library("PNAR")
library("mgcv")
# library("xtable")
library("glmSTARMA")
data("crime")
data("crime_W")
library("tscount")

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


plot(apply(crime, 1, sum) ~ temperature, ylab = "Monthly number of burglaries", xlab = "Temperature (°F)",
     cex = 2.5, cex.lab=2.5, cex.axis=2.5, pch = paste0(rep(1:6, each = 12)))
crime_temp <- lm(apply(crime, 1, sum) ~ temperature)
abline(crime_temp, lwd = 1.5)




crime_sum <- ts(apply(crime, 1, sum), start = c(2010, 1), end = c(2015, 12), frequency = 12)


time <- 1:72
plot(c(crime_sum) ~ time, type = "o")
mod <- gam(c(crime_sum) ~ s(time) + s(temperature))
points(predict(mod) ~ time, type = "l", col = "red")

### Model fitting without covariates:


covariates <- cbind(temperature, 72 - time)
covariates_train <- covariates[1:60,]
covariates_test <- covariates[61:72,]


tic()
fit_univariate <- lapply(seq(ncol(crime)), function(i) tsglm(crime[, i], model = list(past_mean = 1, past_obs = 1), xreg = covariates,
                                                             link = "identity"))
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
                                     link = "identity")
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


saveRDS(result, "../results/result_univariate_linear.rds")
