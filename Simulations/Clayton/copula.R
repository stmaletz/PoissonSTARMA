library("glmSTARMA")
library("copula")

index <- as.numeric(Sys.getenv("PBS_ARRAYID"))
iterations <- 1000

# joe > 1
# frank > 0
# clayton > 0

settings <- expand.grid(copula = c("frank", "clayton", "joe", "independent"), 
                        copula_param = seq(from = 0.5, to = 3, by = 0.5),
                        link = c("log", "identity"),
                        obs = c(50, 100, 250, 500), stringsAsFactors = FALSE)
settings <- subset(settings, !(copula == "joe" & copula_param <= 1))
settings <- subset(settings, !(copula == "independent" & copula_param != 0.5))


W <- generateW("rectangle", dim = 81, 6, width = 9)

## Data generation
model <- list(intercept = "homo", past_obs = 3, past_mean = 1)

if(settings$link[index] == "log"){
  params <- list(intercept = 0.7, past_mean = c(0.2, 0.1), past_obs = rep(0.1, 4))
} else {
  params <- list(intercept = 7, past_mean = c(0.2, 0.1), past_obs = rep(0.1, 4))
}

if(settings$copula[index] == "independent"){
  fam <- vpoisson(link = settings$link[index])
} else {
  fam <- vpoisson(link = settings$link[index], 
                  copula = settings$copula[index],
                  copula_param = settings$copula_param[index])
}





param_est <- vector(mode = "list", length = iterations)
fitting_times <- numeric(iterations)

for(i in seq(iterations)){
  set.seed(i)
  cat("Iteration", i, "\n")
  sim <- glmstarma.sim(settings$obs[index], params, model, W, family = fam)
  fit <- try(glmstarma(sim$observations, model = model, wlist = W, family = fam,
            control = list(parameter_init = "zero", maxit = 1000L)), TRUE)
  if(inherits(fit, "try-error")){
    param_est[[i]] <- NA
    fitting_times[i] <- NA
  } else {
    fitting_times[i] <- fit$fitting_time
    param_est[[i]] <- fit$coefficients_list
  }
}


results <- list(fitting_times = fitting_times, 
                param_est = param_est,
                setting = settings[index,])

saveRDS(results, file = paste0("copula/copula_sim_setting_", index, ".rds"))

