library("glmSTARMA")
library("copula")

index <- as.numeric(Sys.getenv("PBS_ARRAYID"))
iterations <- 1000


settings <- expand.grid(dim = 9,
                        link = c("log", "identity"),
                        obs = c(50, 100, 250, 500),
                        stringsAsFactors = FALSE)


if(settings$link[index] == "log"){
  params <- list(intercept = 0.6 * 1:81 / 81, past_mean = c(0.2, 0.1), past_obs = c(0.2, 0.1))
} else {
  params <- list(intercept = 2 + 3 * 1:81 / 81, past_mean = c(0.2, 0.1), past_obs = c(0.2, 0.1))
}
W <- generateW("rectangle", dim = settings$dim[index]^2, 4, width = settings$dim[index])


model <- list(intercept = "inhomo", past_obs = 1, past_mean = 1)
model_homo <- list(intercept = "homo", past_obs = 1, past_mean = 1)

# Set up family:
fam <- vpoisson(link = settings$link[index], copula = "clayton", copula_param = 2)

param_est <- vector(mode = "list", length = iterations)
param_est_homo <- vector(mode = "list", length = iterations)

fitting_times <- numeric(iterations)
fitting_times_homo <- numeric(iterations)

convergence <- logical(iterations)
convergence_homo <- logical(iterations)

mae_data_inhomo <- numeric(iterations)
mae_data_homo <- numeric(iterations)
mse_data_inhomo <- numeric(iterations)
mse_data_homo <- numeric(iterations)

mae_param_inhomo <- numeric(iterations)
mae_param_homo <- numeric(iterations)
mse_param_inhomo <- numeric(iterations)
mse_param_homo <- numeric(iterations)




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
  fit_inhomo <- try(glmstarma(sim$observations, model = model, wlist = W, family = fam,
                       control = list(parameter_init = "zero", maxit = 10000L,
                                      method = "nloptr", constrained = TRUE)), TRUE)
  fit_homo <- try(glmstarma(sim$observations, model = model_homo, wlist = W, family = fam,
                              control = list(parameter_init = "zero", maxit = 10000L,
                                             method = "nloptr", constrained = TRUE)), TRUE)
  
  if(inherits(fit_inhomo, "try-error")){
    param_est[[i]] <- NA
    fitting_times[i] <- NA
    convergence[i] <- FALSE
    
    mae_data_inhomo[i] <- NA
    mse_data_inhomo[i] <- NA
    mae_param_inhomo[i] <- NA
    mse_param_inhomo[i] <- NA
  } else {
    fitting_times[i] <- fit_inhomo$convergence$fitting_time
    param_est[[i]] <- fit_inhomo$coefficients_list
    convergence[i] <- fit_inhomo$convergence$convergence
    
    mae_data_inhomo[i] <- mean(abs(residuals(fit_inhomo)))
    mse_data_inhomo[i] <- mean(residuals(fit_inhomo)^2)
    bias_est <- c(drop(fit_inhomo$coefficients_list$past_obs) - params$past_obs,
                  drop(fit_inhomo$coefficients_list$past_mean) - params$past_mean)
    mae_param_inhomo[i] <- mean(abs(bias_est))
    mse_param_inhomo[i] <- mean(bias_est^2)
  }
  
  if(inherits(fit_homo, "try-error")){
    param_est_homo[[i]] <- NA
    fitting_times_homo[i] <- NA
    convergence_homo[i] <- FALSE
    
    mae_data_homo[i] <- NA
    mse_data_homo[i] <- NA
    mae_param_homo[i] <- NA
    mse_param_homo[i] <- NA
    
  } else {
    fitting_times_homo[i] <- fit_homo$convergence$fitting_time
    param_est_homo[[i]] <- fit_homo$coefficients_list
    convergence_homo[i] <- fit_homo$convergence$convergence
    
    mae_data_homo[i] <- mean(abs(residuals(fit_homo)))
    mse_data_homo[i] <- mean(residuals(fit_homo)^2)
    bias_est <- c(drop(fit_homo$coefficients_list$past_obs) - params$past_obs,
                  drop(fit_homo$coefficients_list$past_mean) - params$past_mean)
    mae_param_homo[i] <- mean(abs(bias_est))
    mse_param_homo[i] <- mean(bias_est^2)
  }
}


results <- list(params_inhomo = param_est,
                params_homo = param_est_homo,
                fitting_times_inhomo = fitting_times,
                fitting_times_homo = fitting_times_homo,
                convergence_inhomo = convergence,
                convergence_homo = convergence_homo,
                mae_data_inhomo = mae_data_inhomo, 
                mse_data_inhomo = mse_data_inhomo, 
                mae_param_inhomo = mae_param_inhomo, 
                mse_param_inhomo = mse_param_inhomo, 
                mae_data_homo = mae_data_homo, 
                mse_data_homo = mse_data_homo, 
                mae_param_homo = mae_param_homo, 
                mse_param_homo = mse_param_homo, 
                setting = settings[index,])

saveRDS(results, file = paste0("intercept/intercept_sim_setting_", index, ".rds"))



