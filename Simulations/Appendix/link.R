library("glmSTARMA")
library("copula")

index <- as.numeric(Sys.getenv("PBS_ARRAYID"))
iterations <- 1000

settings <- expand.grid(dim = 9,
                        link = c("log", "identity"),
                        param = c("positive", "negative"),
                        obs = c(50, 100, 250, 500),
                        stringsAsFactors = FALSE)
settings <- subset(settings, !(link == "identity" & param == "negative"))

if(settings$link[index] == "log"){
  if(settings$param[index] == "negative"){
    params <- list(intercept = 0.6, past_mean = c(-0.2, 0.1), past_obs = c(0.2, 0.1))
  } else {
    params <- list(intercept = 0.6, past_mean = c(0.2, 0.1), past_obs = c(0.2, 0.1))
  }
} else {
  params <- list(intercept = 5, past_mean = c(0.2, 0.1), past_obs = c(0.2, 0.1))
}
W <- generateW("rectangle", dim = settings$dim[index]^2, 1, width = settings$dim[index])


model <- list(intercept = "homo", past_obs = 1, past_mean = 1)

# Set up family:
fam <- vpoisson(link = settings$link[index], copula = "clayton", copula_param = 2)

param_est_log <- vector(mode = "list", length = iterations)
param_est_linear <- vector(mode = "list", length = iterations)

fitting_times_log <- numeric(iterations)
fitting_times_linear <- numeric(iterations)

convergence_log <- logical(iterations)
convergence_linear <- logical(iterations)

mse_log <- numeric(iterations)
mse_linear <- numeric(iterations)
mae_log <- numeric(iterations)
mae_linear <- numeric(iterations)

aic_log <- numeric(iterations)
bic_log <- numeric(iterations)
qic_log <- numeric(iterations)
aic_linear <- numeric(iterations)
bic_linear <- numeric(iterations)
qic_linear <- numeric(iterations)


for(i in seq(iterations)){
  set.seed(i)
  if(i %% 10 == 0){
    cat("Iteration", i, "\n")
  }
  if(i %% 50 == 0){
    gc()
    Sys.sleep(5)
  }
  sim <- glmstarma.sim(settings$obs[index], params, model, W, family = fam, n_start = 100)
  fit_log <- try(glmstarma(sim$observations, model = model, wlist = W, family = vpoisson("log"),
                              control = list(maxit = 10000L,
                                             method = "nloptr", constrained = TRUE)), TRUE)
  fit_linear <- try(glmstarma(sim$observations, model = model, wlist = W, family = vpoisson("identity"),
                            control = list(maxit = 10000L,
                                           method = "nloptr", constrained = TRUE)), TRUE)
  
  
  if(inherits(fit_log, "try-error")){
    param_est_log[[i]] <- NA
    fitting_times_log[i] <- NA
    aic_log[i] <- NA
    bic_log[i] <- NA
    qic_log[i] <- NA
    mse_log[i] <- NA
    mae_log[i] <- NA
    convergence_log[i] <- FALSE
  } else {
    fitting_times_log[i] <- fit_log$convergence$fitting_time
    param_est_log[[i]] <- fit_log$coefficients_list
    convergence_log[i] <- fit_log$convergence$convergence
    aic_log[i] <- fit_log$aic
    bic_log[i] <- fit_log$bic
    qic_log[i] <- fit_log$qic
    mse_log[i] <- mean(residuals(fit_log)^2)
    mae_log[i] <- mean(abs(residuals(fit_log)))
  }
  
  if(inherits(fit_linear, "try-error")){
    param_est_linear[[i]] <- NA
    fitting_times_linear[i] <- NA
    aic_linear[i] <- NA
    bic_linear[i] <- NA
    qic_log[i] <- NA
    mae_linear[i] <- NA
    mse_linear[i] <- NA
    convergence_linear[i] <- FALSE
  } else {
    fitting_times_linear[i] <- fit_linear$convergence$fitting_time
    param_est_linear[[i]] <- fit_linear$coefficients_list
    convergence_linear[i] <- fit_linear$convergence$convergence
    aic_linear[i] <- fit_linear$aic
    bic_linear[i] <- fit_linear$bic
    qic_linear[i] <- fit_linear$qic
    mse_linear[i] <- mean(residuals(fit_linear)^2)
    mae_linear[i] <- mean(abs(residuals(fit_linear)))
  }
}


results <- list(params_log = param_est_log,
                params_linear = param_est_linear,
                fitting_times_log = fitting_times_log,
                fitting_times_linear = fitting_times_linear,
                convergence_log = convergence_log,
                convergence_linear = convergence_linear,
                mse_log = mse_log,
                mae_log = mae_log, 
                mse_linear = mse_linear,
                mae_linear = mae_linear,
                aic_log = aic_log,
                bic_log = bic_log,
                qic_log = qic_log,
                aic_linear = aic_linear,
                bic_linear = bic_linear,
                qic_linear = qic_linear,
                setting = settings[index,])

saveRDS(results, file = paste0("link/link_sim_setting_", index, ".rds"))

