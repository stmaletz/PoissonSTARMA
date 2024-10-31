library("glmSTARMA")
library("copula")
library("tictoc")

# replicate=1-272
index <- as.numeric(Sys.getenv("PBS_ARRAYID"))
iterations <- 1000

# joe > 1
# frank > 0
# clayton > 0

settings <- expand.grid(dim = c(9, 50),
                        copula = c("frank", "clayton", "joe", "independent"),
                        copula_param = seq(from = 0.5, to = 3, by = 0.5),
                        link = c("log", "identity"),
                        obs = c(5, 10, 20, 50, 100, 250, 500), stringsAsFactors = FALSE)
settings <- subset(settings, !(copula == "joe" & copula_param <= 1))
settings <- subset(settings, !(copula == "independent" & copula_param != 0.5))
settings <- subset(settings, !(dim == 9 & obs < 50))
settings <- subset(settings, !(dim == 50 & obs > 50))
settings <- settings[order(settings$dim),]

W <- generateW("rectangle", dim = settings$dim[index]^2, 2, width = settings$dim[index])

## Data generation
model <- list(intercept = "homo", past_obs = 1, past_mean = 1)

if(settings$link[index] == "log"){
  params <- list(intercept = 0.6, past_mean = c(0.2, 0.1), past_obs = c(0.2, 0.1))
} else {
  params <- list(intercept = 5, past_mean = c(0.2, 0.1), past_obs = c(0.2, 0.1))
}

if(settings$copula[index] == "independent"){
  fam <- vpoisson(link = settings$link[index])
} else {
  fam <- vpoisson(link = settings$link[index], 
                  copula = settings$copula[index],
                  copula_param = settings$copula_param[index])
}





param_est <- vector(mode = "list", length = iterations)
mse_param <- numeric(iterations)
mae_param <- numeric(iterations)
mse_data <- numeric(iterations)
mae_data <- numeric(iterations)
fitting_times <- numeric(iterations)


for(i in seq(iterations)){
  set.seed(i)
  cat("Iteration", i, "\n")
  sim <- glmstarma.sim(settings$obs[index], params, model, W, family = fam)
  tic()
  fit <- try(glmstarma(sim$observations, model = model, wlist = W, family = fam,
            control = list(parameter_init = "zero", maxit = 10000L, method = "nloptr")), TRUE)
  time <- toc(quiet = TRUE)
  if(inherits(fit, "try-error")){
      param_est[[i]] <- NA
      fitting_times[i] <- NA
      mse_data[i] <- NA
      mse_param[i] <- NA
      mae_data[i] <- NA
      mae_param[i] <- NA
  } else {
      fitting_times[i] <- fit$convergence$fitting_time
      param_est[[i]] <- fit$coefficients_list
      
      mse_data[i] <- mean(residuals(fit)^2)
      mae_data[i] <- mean(abs(residuals(fit)))
      
      bias_est <- unlist(params) - coef(fit)
      mse_param[i] <- mean(bias_est^2)
      mae_param[i] <- mean(abs(bias_est))
  }
}


results <- list(fitting_times = fitting_times,
                param_est = param_est,
                mse_data = mse_data,
                mae_data = mae_data,
                mse_param = mse_param,
                mae_param = mae_param,
                setting = settings[index,])

saveRDS(results, file = paste0("copula/copula_sim_setting_", index, ".rds"))

