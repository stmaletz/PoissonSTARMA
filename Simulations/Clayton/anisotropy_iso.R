library("glmSTARMA")
library("copula")
library("tictoc")

# replicate=1-16
index <- as.numeric(Sys.getenv("PBS_ARRAYID"))
iterations <- 1000
# Function for Wald-Test:

settings <- expand.grid(dim = c(9, 50),
                        link = c("log", "identity"),
                        obs = c(5, 10, 20, 50, 100, 250, 500),
                        stringsAsFactors = FALSE)
settings <- subset(settings, !(dim == 50 & obs > 50))
settings <- subset(settings, !(dim == 9 & obs < 50))
settings <- settings[order(settings$dim),]

if(settings$link[index] == "log"){
  params <- list(intercept = 0.6, past_mean = c(0.2, 0.1, 0.05), past_obs = c(0.2, 0.1, 0.05))
} else {
  params <- list(intercept = 5, past_mean = c(0.2, 0.1, 0.05), past_obs = c(0.2, 0.1, 0.05))
}

if(settings$dim[index] == 50){
  params <- list(intercept = 0.6, past_obs = c(0.2, 0.1))
} else {
  params <- list(intercept = 5, past_obs = c(0.2, 0.1))
}


## Generate anisotropic neighborhood:

if(settings$dim[index] == 9){
    coordinates <- expand.grid(x = 1:9, y = 1:9)

    W1 <- matrix(0, 81, 81)
    W2 <- matrix(0, 81, 81)

    for(i in 1:81){
      # Find neighbors:
      W1[i, which((coordinates$y == coordinates$y[i]) & (coordinates$x == (coordinates$x[i] + 1) | coordinates$x == (coordinates$x[i] - 1)))] <- 1
      W2[i, which((coordinates$x == coordinates$x[i]) & (coordinates$y == (coordinates$y[i] + 1) | coordinates$y == (coordinates$y[i] - 1)))] <- 1
    }

    W1 <- W1 / colSums(W1)
    W2 <- W2 / colSums(W2)

    W_anisotropy <- list(diag(81), W1, W2)
    W_isotropy <- generateW("rectangle", dim = settings$dim[index]^2, 4, width = settings$dim[index])
} else {
    coordinates <- expand.grid(x = 1:50, y = 1:50)

    W1 <- matrix(0, 2500, 2500)
    W2 <- matrix(0, 2500, 2500)

    for(i in 1:2500){
      # Find neighbors:
      W1[i, which((coordinates$y == coordinates$y[i]) & (coordinates$x == (coordinates$x[i] + 1) | coordinates$x == (coordinates$x[i] - 1)))] <- 1
      W2[i, which((coordinates$x == coordinates$x[i]) & (coordinates$y == (coordinates$y[i] + 1) | coordinates$y == (coordinates$y[i] - 1)))] <- 1
    }

    W1 <- W1 / colSums(W1)
    W2 <- W2 / colSums(W2)

    W_anisotropy <- list(diag(2500), W1, W2)
    W_isotropy <- generateW("rectangle", dim = settings$dim[index]^2, 2, width = settings$dim[index])
}

if(settings$dim[index] == 9)
{
    model <- list(intercept = "homo", past_obs = 2, past_mean = 2)
    model_isotrop <- list(intercept = "homo", past_obs = 1, past_mean = 1)
} else {
    model <- list(intercept = "homo", past_obs = 2)
    model_isotrop <- list(intercept = "homo", past_obs = 1)
}


# Set up family:
fam <- vpoisson(link = settings$link[index], copula = "clayton", copula_param = 2)

param_est_isotropy <- vector(mode = "list", length = iterations)
fitting_times_isotropy <- numeric(iterations)
convergence_isotropy <- logical(iterations)
residuals_isotropy <- numeric(iterations)# vector("list", iterations)
aic_isotropy <- numeric(iterations)
bic_isotropy <- numeric(iterations)
qic_isotropy <- numeric(iterations)
mse_isotropy <- numeric(iterations)
mae_isotropy <- numeric(iterations)


param_est_anisotropy <- vector(mode = "list", length = iterations)
fitting_times_anisotropy <- numeric(iterations)
convergence_anisotropy <- logical(iterations)
residuals_anisotropy <- numeric(iterations)# vector("list", iterations)
residuals_anisotropy <- numeric(iterations)# vector("list", iterations)
aic_anisotropy <- numeric(iterations)
bic_anisotropy <- numeric(iterations)
qic_anisotropy <- numeric(iterations)
mse_anisotropy <- numeric(iterations)
mae_anisotropy <- numeric(iterations)
# Test for anisotropy:

wald_results <- vector("list", iterations)
wald_results_isotropy <- vector("list", iterations)


# Do Simulation

for(i in seq(iterations)){
  set.seed(i)
  if(i %% 10 == 0){
    cat("Iteration", i, "\n")
  }
  if(i %% 50 == 0){
    gc()
    Sys.sleep(5)
  }
  sim <- glmstarma.sim(settings$obs[index], params, model_isotrop, W_isotropy, family = fam, n_start = 100)
  tic()
  fit_anisotropy <- try(glmstarma(sim$observations, model = model, wlist = W_anisotropy, family = fam,
                           control = list(parameter_init = "zero", maxit = 10000L,
                                          method = "nloptr", constrained = TRUE)), TRUE)
  time_anisotropy <- toc(quiet = TRUE)
  tic()
  fit_isotropy <- try(glmstarma(sim$observations, model = model_isotrop, wlist = W_isotropy, family = fam,
                              control = list(parameter_init = "zero", maxit = 10000L,
                                             method = "nloptr", constrained = TRUE)), TRUE)
  time_isotropy <- toc(quiet = TRUE)
  
  if(inherits(fit_isotropy, "try-error")){
    fitting_times_isotropy[i] <- NA
    convergence_isotropy[i] <- FALSE
    aic_isotropy[i] <- NA
    bic_isotropy[i] <- NA
    qic_isotropy[i] <- NA
    mse_isotropy[i] <- NA
    mae_isotropy[i] <- NA
  } else {
    param_est_isotropy[[i]] <- fit_isotropy$coefficients_list
    fitting_times_isotropy[i] <- time_isotropy$toc - time_isotropy$tic
    convergence_isotropy[i] <- fit_isotropy$convergence$convergence
    residuals_isotropy[i] <- mean(residuals(fit_isotropy)^2)
    aic_isotropy[i] <- fit_isotropy$aic
    bic_isotropy[i] <- fit_isotropy$bic
    qic_isotropy[i] <- fit_isotropy$qic
    mse_isotropy[i] <- mean(residuals(fit_isotropy)^2)
    mae_isotropy[i] <- mean(abs(residuals(fit_isotropy)))
  }
  
  if(inherits(fit_anisotropy, "try-error")){
    fitting_times_anisotropy[i] <- NA
    convergence_anisotropy[i] <- FALSE
    aic_anisotropy[i] <- NA
    bic_anisotropy[i] <- NA
    qic_anisotropy[i] <- NA
    mse_anisotropy[i] <- NA
    mae_anisotropy[i] <- NA
  } else {
    param_est_anisotropy[[i]] <- fit_anisotropy$coefficients_list
    fitting_times_anisotropy[i] <- time_anisotropy$toc - time_anisotropy$tic
    convergence_anisotropy[i] <- fit_anisotropy$convergence$convergence
    residuals_anisotropy[i] <- mean(residuals(fit_anisotropy)^2)
    aic_anisotropy[i] <- fit_anisotropy$aic
    bic_anisotropy[i] <- fit_anisotropy$bic
    qic_anisotropy[i] <- fit_anisotropy$qic
    mse_anisotropy[i] <- mean(residuals(fit_anisotropy)^2)
    mae_anisotropy[i] <- mean(abs(residuals(fit_anisotropy)))
  }
}


results <- list(param_est_isotropy = param_est_isotropy,
                fitting_times_isotropy = fitting_times_isotropy,
                convergence_isotropy = convergence_isotropy,
                residuals_isotropy = residuals_isotropy,
                aic_isotropy = aic_isotropy,
                bic_isotropy = bic_isotropy,
                qic_isotropy = qic_isotropy,
                param_est_anisotropy = param_est_anisotropy,
                fitting_times_anisotropy = fitting_times_anisotropy,
                convergence_anisotropy = convergence_anisotropy,
                residuals_anisotropy = residuals_anisotropy,
                residuals_anisotropy = residuals_anisotropy,
                aic_anisotropy = aic_anisotropy,
                bic_anisotropy = bic_anisotropy,
                qic_anisotropy = qic_anisotropy,
                mse_anisotropy = mse_anisotropy,
                mse_isotropy = mse_isotropy,
                mae_anisotropy = mae_anisotropy,
                mae_isotropy = mae_isotropy,
                setting = settings[index,])

saveRDS(results, file = paste0("anisotropy_ar/anisotropy_ar_sim_setting_", index, ".rds"))




