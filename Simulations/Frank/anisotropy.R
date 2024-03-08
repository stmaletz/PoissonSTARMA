library("glmSTARMA")
library("copula")

index <- as.numeric(Sys.getenv("PBS_ARRAYID"))
iterations <- 1000
# Function for Wald-Test:

isotropy_test <- function(object){
  variance <- variance_estimation(object$model_pointer)
  
  # Test1:
  C0 <- c(0, 0, 1, -1, 0, 0, 0)
  var1 <- drop(t(C0) %*% variance %*% C0)
  statistic1 <- drop(t(C0) %*% coef(object))^2 / var1
  p_val1 <- pchisq(statistic1, df = 1, lower.tail = FALSE)
  
  # Test2:
  C0 <- c(0, 0, 0, 0, 0, 1, -1)
  var2 <- drop(t(C0) %*% variance %*% C0)
  statistic2 <- drop(t(C0) %*% coef(object))^2 / var2
  p_val2 <- pchisq(statistic2, df = 1, lower.tail = FALSE)
  
  return(list(past_mean = list(p.value = p_val1, statistics = statistic1), past_obs = list(p.value = p_val2, statistics = statistic2)))
}


settings <- expand.grid(dim = 9,
                        link = c("log", "identity"),
                        obs = c(50, 100, 250, 500),
                        stringsAsFactors = FALSE)


if(settings$link[index] == "log"){
  params <- list(intercept = 0.6, past_mean = c(0.2, 0.1, 0.05), past_obs = c(0.2, 0.1, 0.05))
} else {
  params <- list(intercept = 5, past_mean = c(0.2, 0.1, 0.05), past_obs = c(0.2, 0.1, 0.05))
}

## Generate anisotropic neighborhood:

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
#

W_isotropy <- generateW("rectangle", dim = settings$dim[index]^2, 4, width = settings$dim[index])



model <- list(intercept = "homo", past_obs = 2, past_mean = 2)
model_isotrop <- list(intercept = "homo", past_obs = 1, past_mean = 1)

# Set up family:
fam <- vpoisson(link = settings$link[index], copula = "frank", copula_param = 1)

param_est_isotropy <- vector(mode = "list", length = iterations)
fitting_times_isotropy <- numeric(iterations)
convergence_isotropy <- logical(iterations)
residuals_isotropy <- numeric(iterations)# vector("list", iterations)
aic_isotropy <- numeric(iterations)
bic_isotropy <- numeric(iterations)
qic_isotropy <- numeric(iterations)


param_est_anisotropy <- vector(mode = "list", length = iterations)
fitting_times_anisotropy <- numeric(iterations)
convergence_anisotropy <- logical(iterations)
residuals_anisotropy <- numeric(iterations)# vector("list", iterations)
residuals_anisotropy <- numeric(iterations)# vector("list", iterations)
aic_anisotropy <- numeric(iterations)
bic_anisotropy <- numeric(iterations)
qic_anisotropy <- numeric(iterations)
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
  sim <- glmstarma.sim(settings$obs[index], params, model, W_anisotropy, family = fam, n_start = 100)
  fit_anisotropy <- try(glmstarma(sim$observations, model = model, wlist = W_anisotropy, family = fam,
                           control = list(parameter_init = "zero", maxit = 1000L,
                                          method = "nloptr", constrained = TRUE)), TRUE)
  fit_isotropy <- try(glmstarma(sim$observations, model = model_isotrop, wlist = W_isotropy, family = fam,
                              control = list(parameter_init = "zero", maxit = 1000L,
                                             method = "nloptr", constrained = TRUE)), TRUE)
  
  
  if(inherits(fit_isotropy, "try-error")){
    fitting_times_isotropy[i] <- NA
    convergence_isotropy[i] <- FALSE
    aic_isotropy[i] <- NA
    bic_isotropy[i] <- NA
    qic_isotropy[i] <- NA
  } else {
    param_est_isotropy[[i]] <- fit_isotropy$coefficients_list
    fitting_times_isotropy[i] <- fit_isotropy$convergence$fitting_time
    convergence_isotropy[i] <- fit_isotropy$convergence$convergence
    residuals_isotropy[i] <- mean(residuals(fit_isotropy)^2)
    aic_isotropy[i] <- fit_isotropy$aic
    bic_isotropy[i] <- fit_isotropy$bic
    qic_isotropy[i] <- fit_isotropy$qic
    summa <- summary(fit_isotropy)$coefficients$`Pr(>|z|)`
    wald_results_isotropy[[i]] <- list(past_mean = summa[3],
                                       past_obs = summa[5])
  }
  
  if(inherits(fit_anisotropy, "try-error")){
    fitting_times_anisotropy[i] <- NA
    convergence_anisotropy[i] <- FALSE
    aic_anisotropy[i] <- NA
    bic_anisotropy[i] <- NA
    qic_anisotropy[i] <- NA
  } else {
    param_est_anisotropy[[i]] <- fit_anisotropy$coefficients_list
    fitting_times_anisotropy[i] <- fit_anisotropy$convergence$fitting_time
    convergence_anisotropy[i] <- fit_anisotropy$convergence$convergence
    residuals_anisotropy[i] <- mean(residuals(fit_anisotropy)^2)
    aic_anisotropy[i] <- fit_anisotropy$aic
    bic_anisotropy[i] <- fit_anisotropy$bic
    qic_anisotropy[i] <- fit_anisotropy$qic
    wald_results[[i]] <- isotropy_test(fit_anisotropy)
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
                wald_results = wald_results,
                wald_results_isotropy = wald_results_isotropy,
                setting = settings[index,])

saveRDS(results, file = paste0("anisotropy/anisotropy_sim_setting_", index, ".rds"))




