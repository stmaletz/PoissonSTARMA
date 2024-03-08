library("data.table")
library("ggplot2")
library("xtable")
library("scales")

# copula <- "Clayton"
copula <- "Frank"



#### Initialisierung

results <- lapply(list.files(paste0(copula, "/init/"), full.names = TRUE), readRDS)
settings <- lapply(results, function(x) x$setting)
settings <- rbindlist(settings)

fitting_times <- sapply(results, function(x) x$fitting_times)
sum(is.na(fitting_times))
which(is.na(fitting_times), arr.ind = TRUE)
# Insgesamt 12 Fehlschlaege im linearen Modell mit 50 Beobachtungen ohne Kovariablen
# bei initialisierung ueber 0
# Fehlschlaege im linearen Modell

calc_mse <- function(elem){
  # browser()
  elem$param_est <- elem$param_est[!sapply(elem$param_est, function(x) all(is.na(x)))]
  runs <- length(elem$param_est)
  # estimations <- rbindlist(estimations)
  
  mse_intercept <- numeric(runs)
  mse_alpha <- numeric(runs)
  mse_beta <- numeric(runs)
  mse_covariates <- numeric(runs)
  mse_total <- numeric(runs)
  
  for(i in 1:runs){
    mse_intercept[i] <- sum(elem$param_est[[i]]$intercept - elem$true_param$intercept)^2
    mse_alpha[i] <- sum(elem$param_est[[i]]$past_mean - elem$true_param$past_mean)^2
    mse_beta[i] <-  sum(elem$param_est[[i]]$past_obs - elem$true_param$past_obs)^2#
    if(!is.null(elem$true_param$covariates)){
      mse_covariates[i] <- sum(elem$param_est[[i]]$covariates - elem$true_param$covariates)^2#
    }
    mse_total[i] <- mse_intercept[i] + mse_alpha[i] + mse_beta[i] + mse_covariates[i]
  }
  mse_alpha <- mse_alpha / length(elem$true_param$past_mean)
  mse_beta <- mse_beta / length(elem$true_param$past_obs)
  df <- 1 + length(elem$true_param$past_mean) + length(elem$true_param$past_obs)
  if(!is.null(elem$true_param$covariates)){
    mse_covariates <- mse_covariates / length(elem$true_param$covariates)
    df <- df + length(elem$true_param$covariates)
  } else {
    mse_covariates <- rep(NA, runs)
  }
  mse_total <- mse_total / df
  return(data.frame(mse_total = mse_total, mse_intercept = mse_intercept, mse_alpha = mse_alpha, mse_beta = mse_beta, mse_covariates = mse_covariates))
}


# Get MSEs

res <- lapply(results, calc_mse)
for(i in seq(length(res))){
  res[[i]] <- cbind(res[[i]], settings[i, ])
}
res <- rbindlist(res)
res <- subset(res, init_link != "transformed_mean")

# New facet label names for order and covariates
order.labs <- c("q = 1", "q = 2")
names(order.labs) <- c("1", "2")
covariate.labs <- c("No covariates", "Covariates")
names(covariate.labs) <- c("FALSE", "TRUE")


init_log <- ggplot(subset(res, obs == 250 & link == "log")) + 
  geom_boxplot(aes(x = init_link, y = log(mse_total), fill = init_link)) + 
  facet_wrap(order~covariate, labeller = labeller(order = order.labs, covariate = covariate.labs)) +
  theme_bw() + xlab("Initialization Method") + ylab("log(MSE)") +
  labs(fill = "Init method") + 
  scale_x_discrete(labels = c("first_obs" = "First Obs.", "mean" = "Mean", "zero" = "Zero", "true" = "True")) + 
  scale_fill_manual(labels = c("first_obs" = "First Obs.", "mean" = "Mean", "zero" = "Zero", "true" = "True"),
                    values = hue_pal()(4)) + 
  theme(legend.position = "none", text = element_text(size = 30))

init_linear <- ggplot(subset(res, obs == 250 & link == "identity")) + 
  geom_boxplot(aes(x = init_link, y = log(mse_total), fill = init_link)) + 
  facet_wrap(order~covariate, labeller = labeller(order = order.labs, covariate = covariate.labs)) +
  theme_bw() + xlab("Initialization Method") + ylab("log(MSE)") +
  labs(fill = "Init method") + 
  scale_x_discrete(labels = c("first_obs" = "First Obs.", "mean" = "Mean", "zero" = "Zero", "true" = "True")) +
  scale_fill_manual(labels = c("first_obs" = "First Obs.", "mean" = "Mean", "zero" = "Zero", "true" = "True"),
                    values = hue_pal()(4)) + 
  theme(legend.position = "none",text = element_text(size = 30))


ggsave(paste0("plots/", copula, "/init_log.pdf"), init_log, width = 12, height = 8) 
ggsave(paste0("plots/", copula, "/init_linear.pdf"), init_linear, width = 12, height = 8)

#







