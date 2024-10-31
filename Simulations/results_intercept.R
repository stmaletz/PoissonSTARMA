library("data.table")
library("ggplot2")
library("xtable")
library("scales")


#### Intercept

results <- lapply(list.files("Appendix/intercept/", full.names = TRUE), readRDS)
settings <- lapply(results, function(x) x$setting)
settings <- rbindlist(settings)

fitting_times <- sapply(results, function(x) x$fitting_times)
sum(is.na(fitting_times))
which(is.na(fitting_times), arr.ind = TRUE)



# Make data.frame
res_data <- lapply(results, function(x){
  data.frame(link = rep(x$setting$link, 1000), obs = rep(x$setting$obs, 1000),
             mae_data_inhomo = x$mae_data_inhomo, mae_data_homo = x$mae_data_homo, 
             mse_param_inhomo = x$mse_param_inhomo, mse_param_homo = x$mse_param_homo)
})

res_data <- rbindlist(res_data)


# Get bias

bias_inhomo <- lapply(results, function(x){
  y <- sapply(x$params_inhomo, unlist)
  if(is.list(y)){
    y <- y[sapply(y, length) == 5]
    y <- t(do.call("rbind", y))
  }
  rowMeans(y)
})
bias_inhomo <- do.call("rbind", bias_inhomo)
bias_inhomo <- bias_inhomo[, -(1:81)]

bias_homo <- lapply(results, function(x){
  y <- sapply(x$params_homo, unlist)
  if(is.list(y)){
    y <- y[sapply(y, length) == 5]
    y <- t(do.call("rbind", y))
  }
  rowMeans(y)
})
bias_homo <- do.call("rbind", bias_homo)
bias_homo <- bias_homo[, -1]


settings$inhomo_bias_alpha1 <- bias_inhomo[, 3] - 0.2
settings$inhomo_bias_alpha2 <- bias_inhomo[, 4] - 0.1
settings$inhomo_bias_beta1 <- bias_inhomo[, 1] - 0.2
settings$inhomo_bias_beta2 <- bias_inhomo[, 2] - 0.1

settings$homo_bias_alpha1 <- bias_homo[, 3] - 0.2
settings$homo_bias_alpha2 <- bias_homo[, 4] - 0.1
settings$homo_bias_beta1 <- bias_homo[, 1] - 0.2
settings$homo_bias_beta2 <- bias_homo[, 2] - 0.1

settings$dim <- NULL

xtable(settings[order(settings$link, settings$obs),], digits = 3)





###

res_data2 <- res_data

res_data$mae <- res_data$mae_data_homo
res_data$model <- "homogenous"
res_data2$mae <- res_data2$mae_data_inhomo
res_data2$model <- "inhomogenous"
res <- rbind(res_data, res_data2)


plot_mae_data_linear <- ggplot(subset(res, link == "identity")) + 
    geom_boxplot(aes(x = obs, group = interaction(obs, model), y = mae, fill = model)) + 
    labs(x = "Observation times", y = "MAE", fill = "Intercept") + 
    theme_bw() + theme(legend.position = "bottom", text = element_text(size = 30))

plot_mae_data_log <- ggplot(subset(res, link == "log")) + 
  geom_boxplot(aes(x = obs, group = interaction(obs, model), y = mae, fill = model)) + 
  labs(x = "Observation times", y = "MAE", fill = "Intercept") + 
  theme_bw() + theme(legend.position = "bottom", text = element_text(size = 30))  

ggsave(paste0("plots/Appendix/intercept/intercept_mae_data_log.pdf"), plot_mae_data_log, width = 12, height = 10)  
ggsave(paste0("plots/Appendix/intercept/intercept_mae_data_linear.pdf"), plot_mae_data_linear, width = 12, height = 10)
















