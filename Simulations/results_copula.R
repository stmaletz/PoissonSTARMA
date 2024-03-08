library("data.table")
library("ggplot2")
library("xtable")
library("scales")


#### Copula

results <- lapply(list.files("Appendix/copula/", full.names = TRUE), readRDS)
settings <- lapply(results, function(x) x$setting)
settings <- rbindlist(settings)

fitting_times <- sapply(results, function(x) x$fitting_times)
sum(is.na(fitting_times))
which(is.na(fitting_times), arr.ind = TRUE)
# Insgesamt 2 Fehlschlaege im linearen Modell mit 50 Beobachtungen


# Make data.frame
res_data <- lapply(results, function(x){
  data.frame(copula = x$setting$copula, param = x$setting$copula_param,
             link = x$setting$link, obs = x$setting$obs,
             mae_data = x$mae_data, mse_data = x$mse_data,
             mae_param = x$mae_param, mse_param = x$mse_param)
})


res_data <- rbindlist(res_data)


## Plot linear
res_data$copula <- gsub("clayton", "Clayton", res_data$copula)
res_data$copula <- gsub("frank", "Frank", res_data$copula)
res_data$copula <- gsub("joe", "Joe", res_data$copula)

for(t in c(50, 100, 250, 500)){
  plot_mae_data <- ggplot(subset(res_data, obs == t & link == "identity" & copula != "independent")) + 
    geom_boxplot(aes(x = param, group = param, y = mae_data)) + 
    facet_wrap(~copula) + 
    labs(x = "Copula Parameter", y = "MAE in Data", fill = "Copula") + 
    theme_bw() + theme(legend.position = "none",text = element_text(size = 30))
  
  plot_mse_data <- ggplot(subset(res_data, obs == t & link == "identity" & copula != "independent")) + 
    geom_boxplot(aes(x = param, group = param, y = mse_data)) + 
    facet_wrap(~copula) + 
    labs(x = "Copula Parameter", y = "MSE in Data", fill = "Copula") + 
    theme_bw() + theme(legend.position = "none",text = element_text(size = 30))
  
  plot_mae_param <- ggplot(subset(res_data, obs == t & link == "identity" & copula != "independent")) + 
    geom_boxplot(aes(x = param, group = param, y = mae_param)) + 
    facet_wrap(~copula) + 
    labs(x = "Copula Parameter", y = "MAE of parameter estimation", fill = "Copula") + 
    theme_bw() + theme(legend.position = "none",text = element_text(size = 30))
  
  plot_mse_param <- ggplot(subset(res_data, obs == t & link == "identity" & copula != "independent")) + 
    geom_boxplot(aes(x = param, group = param, y = mse_param)) + 
    facet_wrap(~copula) + 
    labs(x = "Copula Parameter", y = "MSE of parameter estimation", fill = "Copula") + 
    theme_bw() + theme(legend.position = "none",text = element_text(size = 30))
  
  ggsave(paste0("plots/Appendix/copula/mae_data_linear_obs_", t, ".pdf"), plot_mae_data, width = 12, height = 6) 
  ggsave(paste0("plots/Appendix/copula/mse_data_linear_obs_", t, ".pdf"), plot_mse_data, width = 12, height = 6) 
  ggsave(paste0("plots/Appendix/copula/mae_param_linear_obs_", t, ".pdf"), plot_mae_param, width = 12, height = 6) 
  ggsave(paste0("plots/Appendix/copula/mse_param_linear_obs_", t, ".pdf"), plot_mse_param, width = 12, height = 6) 
}





## Plot log-linear


for(t in c(50, 100, 250, 500)){
  plot_mae_data <- ggplot(subset(res_data, obs == t & link == "log" & copula != "independent")) + 
    geom_boxplot(aes(x = param, group = param, y = mae_data)) + 
    facet_wrap(~copula) + 
    labs(x = "Copula Parameter", y = "MAE in Data", fill = "Copula") + 
    theme_bw() + theme(legend.position = "none",text = element_text(size = 30))
  
  plot_mse_data <- ggplot(subset(res_data, obs == t & link == "log" & copula != "independent")) + 
    geom_boxplot(aes(x = param, group = param, y = mse_data)) + 
    facet_wrap(~copula) + 
    labs(x = "Copula Parameter", y = "MSE in Data", fill = "Copula") + 
    theme_bw() + theme(legend.position = "none",text = element_text(size = 30))
  
  plot_mae_param <- ggplot(subset(res_data, obs == t & link == "log" & copula != "independent")) + 
    geom_boxplot(aes(x = param, group = param, y = mae_param)) + 
    facet_wrap(~copula) + 
    labs(x = "Copula Parameter", y = "MAE of parameter estimation", fill = "Copula") + 
    theme_bw() + theme(legend.position = "none",text = element_text(size = 30))
  
  plot_mse_param <- ggplot(subset(res_data, obs == t & link == "log" & copula != "independent")) + 
    geom_boxplot(aes(x = param, group = param, y = mse_param)) + 
    facet_wrap(~copula) + 
    labs(x = "Copula Parameter", y = "MSE of parameter estimation", fill = "Copula") + 
    theme_bw() + theme(legend.position = "none",text = element_text(size = 30))
  
  ggsave(paste0("plots/Appendix/copula/mae_data_log_obs_", t, ".pdf"), plot_mae_data, width = 12, height = 6) 
  ggsave(paste0("plots/Appendix/copula/mse_data_log_obs_", t, ".pdf"), plot_mse_data, width = 12, height = 6) 
  ggsave(paste0("plots/Appendix/copula/mae_param_log_obs_", t, ".pdf"), plot_mae_param, width = 12, height = 6) 
  ggsave(paste0("plots/Appendix/copula/mse_param_log_obs_", t, ".pdf"), plot_mse_param, width = 12, height = 6) 
}


bias <- lapply(results, function(x){
  y <- sapply(x$param_est, unlist)
  if(is.list(y)){
    y <- y[sapply(y, length) == 5]
    y <- t(do.call("rbind", y))
  }
  rowMeans(y)
})
bias <- do.call("rbind", bias)


variance <- lapply(results, function(x){
  y <- sapply(x$param_est, unlist)
  if(is.list(y)){
    y <- y[sapply(y, length) == 5]
    y <- t(do.call("rbind", y))
  }
  apply(y, 1, var)
})
variance <- do.call("rbind", variance)
variance <- rowMeans(variance)

settings$bias_intercept <- bias[, 1] - ifelse(settings$link == "identity", 5, 0.6)
settings$bias_alpha1 <- bias[, 4] - 0.2
settings$bias_alpha2 <- bias[, 5] - 0.1
settings$bias_beta1 <- bias[, 2] - 0.2
settings$bias_beta2 <- bias[, 3] - 0.1
settings$bias_squared <- (settings$bias_intercept^2 + settings$bias_alpha1^2 + settings$bias_alpha2^2 + settings$bias_beta1^2 + settings$bias_beta2^2) / 5
settings$mse <- sapply(results, function(x) mean(x$mse_param, na.rm = TRUE))
settings$variance <- settings$mse - settings$bias_squared
settings$variance2 <- variance


tab <- subset(settings, copula == "clayton" & obs == 250 & link == "identity")
tab <- tab[order(tab$copula_param)]
tab$copula <- NULL
tab$link <- NULL
tab$obs <- NULL
tab

subset(settings, copula == "independent" & obs == 250 & link == "identity")
