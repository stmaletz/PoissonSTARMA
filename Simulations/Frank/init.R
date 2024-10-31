library("glmSTARMA")
library("copula")
library("tictoc")

index <- as.numeric(Sys.getenv("PBS_ARRAYID"))
iterations <- 1000

# replicate=1-288
settings <- expand.grid(init_link = c("first_obs", "mean", "transformed_mean", "zero", "true"),
                        link = c("log", "identity"),
                        obs = c(5, 10, 20, 50, 100, 250, 500),
                        dim = c(9, 50),
                        order = c(1, 2),
                        covariate = c(TRUE, FALSE),
                        stringsAsFactors = FALSE)
settings <- subset(settings, !(init_link == "transformed_mean" & link == "identity"))
settings <- subset(settings, !(dim == 9 & obs < 50))
settings <- subset(settings, !(dim == 50 & obs > 50))
settings <- settings[order(settings$dim),]

W <- generateW("rectangle", dim = settings$dim[index]^2, 3, width = settings$dim[index])


# Generate covariate process
set.seed(42)
covariates <- t(replicate(settings$dim[index]^2,
                          arima.sim(n = 1000,
                                    list(ar = c(0.89, -0.3), ma = c(-0.1, 0.28)),
                                    rand.gen = runif)))
covariates <- covariates / max(covariates)
covariates[covariates < 0] <- 0
covariates <- covariates[seq(settings$dim[index]),]
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
fam <- vpoisson(link = settings$link[index], copula = "frank", copula_param = 1)

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
    tic()
    fit <- try(glmstarma(sim$observations, model = model, wlist = W, covariates = covariate, family = fam,
                         control = list(parameter_init = "zero", maxit = 10000L,
                                        init_link = sim$link_values[, seq(settings$order[index]), drop = FALSE],
                                        method = "nloptr", constrained = TRUE)), TRUE)
    time <- toc(quiet = TRUE)
  } else {
    tic()
    fit <- try(glmstarma(sim$observations, model = model, wlist = W, covariates = covariate, family = fam,
                                     control = list(parameter_init = "zero", maxit = 10000L,
                                                    init_link = settings$init_link[index],
                                                    method = "nloptr", constrained = TRUE)), TRUE)
    time <- toc(quiet = TRUE)
  }
  if(inherits(fit, "try-error")){
    param_est[[i]] <- NA
    fitting_times[i] <- NA
    converged[i] <- FALSE
  } else {
    fitting_times[i] <- time$toc - time$tic
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
