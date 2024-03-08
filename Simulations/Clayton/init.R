library("glmSTARMA")
library("copula")

index <- as.numeric(Sys.getenv("PBS_ARRAYID"))
iterations <- 1000

# 144 settings

settings <- expand.grid(init_link = c("first_obs", "mean", "transformed_mean", "zero", "true"),
                        link = c("log", "identity"),
                        obs = c(50, 100, 250, 500),
                        order = c(1, 2),
                        covariate = c(TRUE, FALSE),
                        stringsAsFactors = FALSE)
settings <- subset(settings, !(init_link == "transformed_mean" & link == "identity"))
W <- generateW("rectangle", dim = 81, 6, width = 9)


# Generate covariate process
set.seed(42)
covariates <- t(replicate(81,
                          arima.sim(n = 1000, 
                                    list(ar = c(0.89, -0.3), ma = c(-0.1, 0.28)), 
                                    rand.gen = runif)))
covariates <- covariates / max(covariates)
covariates <- covariates[, seq(settings$obs[index])]

if(settings$covariate[index]){
  covariate <- list(covariates)
} else {
  covariate <- list()
}


## Set model and params
if(settings$order[index] == 1 && settings$covariate[index]){
  model <- list(intercept = "homo", past_obs = 1, past_mean = 1, covariates = 0)
  if(settings$link[index] == "log"){
    params <- list(intercept = 0.6, past_mean = c(0.2, 0.1), past_obs = c(0.2, 0.1), covariates = 0.9)
  } else {
    params <- list(intercept = 5, past_mean = c(0.2, 0.1), past_obs = c(0.2, 0.1), covariates = 2)
  }
} else if(settings$order[index] == 2 && settings$covariate[index]){
  model <- list(intercept = "homo", past_obs = c(1, 1), past_mean = c(1, 1), covariates = 0)
  if(settings$link[index] == "log"){
    params <- list(intercept = 0.6, past_mean = c(0.2, 0.1, 0.05, 0.05), past_obs = c(0.2, 0.1, 0.05, 0.05), covariates = 0.9)
  } else {
    params <- list(intercept = 5, past_mean = c(0.2, 0.1, 0.05, 0.05), past_obs = c(0.2, 0.1, 0.05, 0.05), covariates = 2)
  }
} else if(settings$order[index] == 1 && !settings$covariate[index]){
  model <- list(intercept = "homo", past_obs = 1, past_mean = 1)
  if(settings$link[index] == "log"){
    params <- list(intercept = 0.6, past_mean = c(0.2, 0.1), past_obs = c(0.2, 0.1))
  } else {
    params <- list(intercept = 5, past_mean = c(0.2, 0.1), past_obs = c(0.2, 0.1))
  }
} else {
  model <- list(intercept = "homo", past_obs = c(1, 1), past_mean = c(1, 1))
  if(settings$link[index] == "log"){
    params <- list(intercept = 0.6, past_mean = c(0.2, 0.1, 0.05, 0.05), past_obs = c(0.2, 0.1, 0.05, 0.05))
  } else {
    params <- list(intercept = 5, past_mean = c(0.2, 0.1, 0.05, 0.05), past_obs = c(0.2, 0.1, 0.05, 0.05))
  }
}


# Set up family:
fam <- vpoisson(link = settings$link[index], copula = "clayton", copula_param = 2)

param_est <- vector(mode = "list", length = iterations)
param_est_constrained <- vector(mode = "list", length = iterations)
fitting_times <- numeric(iterations)
fitting_times_constrained <- numeric(iterations)
converged <- logical(iterations)
converged_constrained <- logical(iterations)

for(i in seq(iterations)){
  set.seed(i)
  if(i %% 10 == 0){
    cat("Iteration", i, "\n")
  }
  if(i %% 50){
    gc()
    Sys.sleep(10)
  }
  sim <- glmstarma.sim(settings$obs[index], params, model, W, covariate, family = fam, n_start = 100)
  if(settings$init_link[index] == "true"){
    fit <- try(glmstarma(sim$observations, model = model, wlist = W, covariates = covariate, family = fam,
                         control = list(parameter_init = "zero", maxit = 1000L,
                                        init_link = sim$link_values[, seq(settings$order[index]), drop = FALSE],
                                        method = "nloptr", constrained = TRUE)), TRUE)
  } else {
    fit <- try(glmstarma(sim$observations, model = model, wlist = W, covariates = covariate, family = fam,
                                     control = list(parameter_init = "zero", maxit = 1000L,
                                                    init_link = settings$init_link[index],
                                                    method = "nloptr", constrained = TRUE)), TRUE)
  }
  if(inherits(fit, "try-error")){
    param_est[[i]] <- NA
    fitting_times[i] <- NA
    converged[i] <- FALSE
  } else {
    fitting_times[i] <- fit$convergence$fitting_time
    param_est[[i]] <- fit$coefficients_list
    converged[i] <- fit$convergence$convergence
  }
}

results <- list(fitting_times = fitting_times, 
                param_est = param_est,
                converged = converged,
                true_param = params,
                setting = settings[index,])

saveRDS(results, file = paste0("init/init_sim_setting_", index, ".rds"))
