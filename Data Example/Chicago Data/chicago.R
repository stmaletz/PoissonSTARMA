library("PNAR")
library("mgcv")
library("xtable")
library("glmSTARMA")
data("crime")
data("crime_W")

population <- read.csv("pop.csv")$x
unemployed <- read.csv("unemp.csv")$x
wealth <- read.csv("wealth.csv")$x
young_males <- read.csv("ym.csv")$x
temperature <- read.csv2("temperature.csv")$Temperatur


# EDA
sum(crime)
sum(crime) / 552
sum(crime) / 72
max(crime)
which(crime == max(crime), arr.ind = TRUE)


pdf("plots/crime_temperature.pdf", width = 11, height = 9)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(apply(crime, 1, sum) ~ temperature, ylab = "Monthly number of burglaries", xlab = "Temperature (°F)",
     cex = 2.5, cex.lab=2.5, cex.axis=2.5)
crime_temp <- lm(apply(crime, 1, sum) ~ temperature)
abline(crime_temp, lwd = 1.5)
dev.off()

crime_sum <- ts(apply(crime, 1, sum), start = c(2010, 1), end = c(2015, 12), frequency = 12)


time <- 1:72
plot(c(crime_sum) ~ time, type = "o")
mod <- gam(c(crime_sum) ~ s(time) + s(temperature))
points(predict(mod) ~ time, type = "l", col = "red")

### Model fitting without covariates:

orders <- expand.grid(obs = c(1, 2), means = c(0, 1, 2), link = c("identity", "log"), stringsAsFactors = FALSE)
models <- vector(mode = "list", length = nrow(orders))
for(i in seq(length(models))){
  if(orders$means[i] > 0){
    models[[i]] <- glmstarma(t(crime), model = list(past_obs = rep(1, orders$obs[i]), 
                                                    past_mean = rep(1, orders$means[i])),
                             wlist = list(diag(552), crime_W), family = vpoisson(orders$link[i]),
                             control = list(method = "nloptr", maxit = 1000))
  } else {
    models[[i]] <- glmstarma(t(crime), model = list(past_obs = rep(1, orders$obs[i])),
                             wlist = list(diag(552), crime_W), family = vpoisson(orders$link[i]),
                             control = list(method = "nloptr", maxit = 1000))
  }
}

orders$qic <- sapply(models, QIC)


### Modelfitting with covariates
covariates_log <- list(population = TimeConstant(log(population / 1000)), 
                   unemployed = TimeConstant(unemployed),
                   wealth = TimeConstant(log1p(exp(wealth))),
                   youngmales = TimeConstant(young_males / population),
                   temperature = SpatialConstant(temperature), 
                   trend = SpatialConstant(72:1)
)

covariates_linear <- list(population = TimeConstant(population / 1000), 
                          unemployed = TimeConstant(unemployed),
                          wealth = TimeConstant(log1p(exp(wealth))),
                          youngmales = TimeConstant(young_males / population),
                          temperature = SpatialConstant(temperature), 
                          trend = SpatialConstant(72:1)
)



orders_covariates <- expand.grid(obs = c(1, 2), means = c(0, 1, 2), link = c("identity", "log"), stringsAsFactors = FALSE)
models_covariates <- vector(mode = "list", length = nrow(orders_covariates))
for(i in seq(length(models))){
  covariates <- covariates_linear
  meth <- "optim"
  if(orders_covariates$link[i] == "log"){
    covariates <- covariates_log
    meth <- "nloptr"
  }
  if(orders_covariates$means[i] > 0){
    models_covariates[[i]] <- glmstarma(t(crime), model = list(past_obs = rep(1, orders_covariates$obs[i]), 
                                                    past_mean = rep(1, orders_covariates$means[i])),
                             wlist = list(diag(552), crime_W), covariates = covariates,
                             family = vpoisson(orders_covariates$link[i]),
                             control = list(method = meth, maxit = 1000))
  } else {
    models_covariates[[i]] <- glmstarma(t(crime), model = list(past_obs = rep(1, orders_covariates$obs[i])),
                                        covariates = covariates, wlist = list(diag(552), crime_W), family = vpoisson(orders_covariates$link[i]),
                                        control = list(method = "nloptr", maxit = 1000))
  }
}

orders$qic_covariates <- sapply(models_covariates, QIC)


sapply(models, function(x) mean(residuals(x)^2))
sapply(models_covariates, function(x) mean(residuals(x)^2))


models_considered <- models[(orders$obs == orders$means)]
models_considered <- c(models_considered, models_covariates[orders$obs == orders$means])
models_considered <- models_considered[c(1, 5, 2, 6, 3, 7, 4, 8)]

subset(orders, obs == means)


extract_coefs <- function(models){
  params <- matrix(0, nrow = 16, length(models))
  for(i in seq(length(models))){
    params[1, i] <- models[[i]]$coefficients_list$intercept
    if(!is.null(models[[i]]$coefficients_list$past_mean)){
      params[2:5, i] <- c(models[[i]]$coefficients_list$past_mean)
    }
    params[6:9, i] <- c(models[[i]]$coefficients_list$past_obs)
    if(!is.null(models[[i]]$coefficients_list$covariates)){
      params[10:15, i] <- c(models[[i]]$coefficients_list$covariates) 
    }
    params[16, i] <- models[[i]]$qic / 1000
  }
  return(params)
}
xtable(extract_coefs(models_considered), digits = 4)

summaries <- lapply(models_considered, summary)

round(extract_coefs(models_considered), digits = 4)


subset(orders, (obs == 2 & means == 0) | (obs == 2 & means == 2))




####

pdf("plots/crime_accu.pdf", width = 11, height = 9)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(rowSums(crime), ylab = "Accumulated Buglaries", xlab = "Month",
     cex = 2.5, cex.lab=2.5, cex.axis=2.5)
points(colSums(fitted.values(models_considered[[1]])), type = "l", col = "lightblue", lwd = 1.5)
points(colSums(fitted.values(models_considered[[8]])), type = "l", col = "orange", lwd = 1.5)
legend("topright", legend = c("linear", "log-linear"), lty = 1, col = c("lightblue", "orange"),
       cex = 2.5, lwd = 1.5)

dev.off()
##
