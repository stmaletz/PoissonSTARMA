library("data.table")
library("ggplot2")
library("xtable")
library("scales")


#### Link

results <- lapply(list.files("Appendix/link/", full.names = TRUE), readRDS)
settings <- lapply(results, function(x) x$setting)
settings <- rbindlist(settings)

fitting_times <- sapply(results, function(x) x$fitting_times)
sum(is.na(fitting_times))
which(is.na(fitting_times), arr.ind = TRUE)


# Make data.frame
res_data <- lapply(results, function(x){
  data.frame(link = rep(x$setting$link, 1000), obs = rep(x$setting$obs, 1000),
             param = rep(x$setting$param, 1000),
             mae_data_linear = x$mae_linear, mae_data_log = x$mae_log, 
             qic_linear = x$qic_linear, qic_log = x$qic_log)
})

res_data <- rbindlist(res_data)
res_data$log_preferred <- res_data$qic_linear > res_data$qic_log

aggregate(cbind(mae_data_linear, mae_data_log, log_preferred, qic_linear, qic_log) ~ param + link + obs, data = res_data, FUN = mean)

 


## plotting
plot_data <- res_data
plot_data$mae <- res_data$mae_data_log
plot_data$model_fit <- "log"

pl2 <- res_data
pl2$mae <- res_data$mae_data_linear
pl2$model_fit <- "linear"
plot_data <- rbind(plot_data, pl2)


ggplot(subset(plot_data)) + 
  geom_boxplot(aes(x = obs, group = interaction(obs, model_fit), y = mae, fill = model_fit)) + 
  labs(x = "Observation times", y = "MAE in Data", fill = "Intercept") + 
  facet_wrap(~interaction(link, param))
  theme_bw() + theme(legend.position = "bottom", text = element_text(size = 30)) 






# Get bias

bias_log <- lapply(results, function(x){
  y <- sapply(x$params_log, unlist)
  if(is.list(y)){
    y <- y[sapply(y, length) == 5]
    y <- t(do.call("rbind", y))
  }
  rowMeans(y)
})
bias_log <- do.call("rbind", bias_log)
bias_log <- bias_log[,-1]

bias_linear <- lapply(results, function(x){
  y <- sapply(x$params_linear, unlist)
  if(is.list(y)){
    y <- y[sapply(y, length) == 5]
    y <- t(do.call("rbind", y))
  }
  rowMeans(y)
})
bias_linear <- do.call("rbind", bias_linear)
bias_linear <- bias_linear[, -1]



settings$log_bias_alpha1 <- bias_log[, 3] - ifelse(settings$param == "positive", 0.2, -0.2)
settings$log_bias_alpha2 <- bias_log[, 4] - 0.1
settings$log_bias_beta1 <- bias_log[, 1] - 0.2
settings$log_bias_beta2 <- bias_log[, 2] - 0.1

settings$linear_bias_alpha1 <- bias_linear[, 3] - ifelse(settings$param == "positive", 0.2, -0.2)
settings$linear_bias_alpha2 <- bias_linear[, 4] - 0.1
settings$linear_bias_beta1 <- bias_linear[, 1] - 0.2
settings$linear_bias_beta2 <- bias_linear[, 2] - 0.1

settings[order(settings$link, settings$param, settings$obs),]




res_data2 <- res_data

res_data$mae <- res_data$mae_data_homo
res_data$model <- "homogenous"
res_data2$mae <- res_data2$mae_data_inhomo
res_data2$model <- "inhomogenous"
res <- rbind(res_data, res_data2)


plot_mae_data_linear <- ggplot(subset(res, link == "identity")) + 
    geom_boxplot(aes(x = obs, group = interaction(obs, model), y = mae, fill = model)) + 
    labs(x = "Observation times", y = "MAE in Data", fill = "Intercept") + 
    theme_bw() + theme(legend.position = "bottom", text = element_text(size = 30))

plot_mae_data_log <- ggplot(subset(res, link == "log")) + 
  geom_boxplot(aes(x = obs, group = interaction(obs, model), y = mae, fill = model)) + 
  labs(x = "Observation times", y = "MAE in Data", fill = "Intercept") + 
  theme_bw() + theme(legend.position = "bottom", text = element_text(size = 30))  

ggsave(paste0("plots/Appendix/intercept/intercept_mae_data_log.pdf"), plot_mae_data_log, width = 12, height = 10)  
ggsave(paste0("plots/Appendix/intercept/intercept_mae_data_linear.pdf"), plot_mae_data_linear, width = 12, height = 10)
















