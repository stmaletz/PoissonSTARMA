library("glmSTARMA")
library("copula")
library("tictoc")

# replicate=1-4590
index <- as.numeric(Sys.getenv("PBS_ARRAYID"))
iterations <- 1000

# parameter_value <- seq(from = -1, to = 3, by = 0.02)
parameter_value <- seq(from = -0.25, to = 0.5, by = 0.01)
settings <- expand.grid(dim = c(5, 7, 9, 50),
                        link = c("log", "identity"),
                        obs = c(5, 10, 20, 50, 100, 250),
                        covariate = c(TRUE, FALSE),
                        test = c("alpha", "beta", "covariate"),
                        param_value = parameter_value,
                        stringsAsFactors = FALSE)
settings <- subset(settings, !(!covariate & test == "covariate"))
settings <- subset(settings, !(link == "identity" & param_value < 0))
settings <- subset(settings, !(link == "log" & param_value > 0.25))
settings <- subset(settings, !(dim < 50 & obs < 100))
settings <- subset(settings, !(dim == 50 & obs > 50))
settings <- settings[order(settings$dim, settings$obs),]

W <- generateW("rectangle", dim = settings$dim[index]^2, 3, width = settings$dim[index])


# Generate covariate process
set.seed(42)
covariates <- t(replicate(settings$dim[index]^2,
                          arima.sim(n = 1000,
                                    list(ar = c(0.89, -0.3), ma = c(-0.1, 0.28)),
                                    rand.gen = runif)))
covariates <- covariates / max(covariates)
covariates[covariates < 0] <- 0
covariates <- covariates[, seq(settings$obs[index])]
if(settings$covariate[index]){
  covariate <- list(covariates)
} else {
  covariate <- list()
}


## Set model and params
if(settings$covariate[index]){
  model <- list(intercept = "homo", past_obs = 1, past_mean = 1, covariates = 0)
  if(settings$link[index] == "log"){
    params <- list(intercept = 0.6, past_mean = c(0.2, 0.1), past_obs = c(0.2, 0.1), covariates = 0.9)
  } else {
    params <- list(intercept = 5, past_mean = c(0.2, 0.1), past_obs = c(0.2, 0.1), covariates = 2)
  }
} else if(!settings$covariate[index]){
  model <- list(intercept = "homo", past_obs = 1, past_mean = 1)
  if(settings$link[index] == "log"){
    params <- list(intercept = 0.6, past_mean = c(0.2, 0.1), past_obs = c(0.2, 0.1))
  } else {
    params <- list(intercept = 5, past_mean = c(0.2, 0.1), past_obs = c(0.2, 0.1))
  }
}


## Adjust parameters for testing

if(settings$test[index] == "alpha"){
  params$past_mean <- c(settings$param_value[index], 0.1)
}
if(settings$test[index] == "beta"){
  params$past_obs <- c(settings$param_value[index], 0.1)
}
if(settings$test[index] == "covariate"){
  params$covariates <- settings$param_value[index]
}


# Set up family:
fam <- vpoisson(link = settings$link[index], copula = "clayton", copula_param = 2)

param_est <- vector(mode = "list", length = iterations)
fitting_times <- numeric(iterations)
wald_result <- vector(mode = "list", length = iterations)
convergence <- logical(iterations)


for(i in seq(iterations)){
  set.seed(i)
  if(i %% 10 == 0){
    cat("Iteration", i, "\n")
  }
  if(i %% 50 == 0){
    gc()
    Sys.sleep(5)
  }
  
  sim <- glmstarma.sim(settings$obs[index], params, model, W, covariate, family = fam, n_start = 100)
  tic()
  fit <- try(glmstarma(sim$observations, model = model, wlist = W, covariates = covariate, family = fam,
                       control = list(parameter_init = "zero", maxit = 10000L,
                                      method = "nloptr", constrained = TRUE)), TRUE)
  time <- toc(quiet = TRUE)
  summa <- try(summary(fit)$coefficients, TRUE)
  
  if(!inherits(summa, "try-error") && settings$test[index] == "alpha"){
    wald_result[[i]] <- list(test_statistic = summa$`z value`[2]^2, p_value = summa$`Pr(>|z|)`[2])
  }
  if(!inherits(summa, "try-error") && settings$test[index] == "beta"){
    wald_result[[i]] <- list(test_statistic = summa$`z value`[4]^2, p_value = summa$`Pr(>|z|)`[4])
  }
  if(!inherits(summa, "try-error") && settings$test[index] == "covariate"){
    wald_result[[i]] <- list(test_statistic = summa$`z value`[6]^2, p_value = summa$`Pr(>|z|)`[6])
  }
  
  if(inherits(fit, "try-error")){
    param_est[[i]] <- NA
    fitting_times[i] <- NA
    convergence[i] <- FALSE
  } else {
    fitting_times[i] <- time$toc - time$tic
    param_est[[i]] <- fit$coefficients_list
    convergence[i] <- fit$convergence$convergence
  }
}

results <- list(fitting_times = fitting_times, 
                param_est = param_est,
                wald_result = wald_result,
                convergence = convergence,
                true_param = params,
                setting = settings[index,])

saveRDS(results, file = paste0("power/power_sim_setting_", index, ".rds"))
