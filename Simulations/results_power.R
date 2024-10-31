library("data.table")
library("ggplot2")
library("xtable")
library("scales")

copula <- "Clayton"
# copula <- "Frank"

## Function for power calculation
power_calc <- function(elem){
  mean(unlist(sapply(elem$wald_result, "[[", "p_value")) < 0.05, na.rm = TRUE)
}

power_calc2 <- function(elem){
  x <- sapply(elem$wald_result, "[[", "p_value")
  y <- !sapply(x, is.null)
  x <- x[y]
  conv <- elem$convergence[y]
  mean(x < 0.05, na.rm = TRUE)
}


## read in data
results <- lapply(list.files(paste0(copula, "/power"), full.names = TRUE), readRDS)
settings <- rbindlist(lapply(results, "[[", "setting"))

params <- lapply(results, "[[", "param_est")
params <- lapply(params, function(x) lapply(x, unlist))
for(i in seq_along(params)){
  if(results[[i]]$setting$link == "log"){
    params[[i]] <- which(sapply(params[[i]], function(y) length(y) == 1 | y[6] == 0))
  } else {
    params[[i]] <- which(sapply(params[[i]], function(y) length(y) == 1 | y[6] == 1))
  }
}

settings$errors <- sapply(params, length) / 1000

for(i in seq_along(results)){
  if(length(params[[i]]) > 0){
    results[[i]]$fitting_times <- results[[i]]$fitting_times[-params[[i]]]
    results[[i]]$param_est <- results[[i]]$param_est[-params[[i]]]
    results[[i]]$wald_result <- results[[i]]$wald_result[-params[[i]]]
  }
}


settings$power <- sapply(results, power_calc)
# settings$power <- sapply(results, power_calc2)

set2 <- settings
set2$power <- NULL
set2$param_value <- NULL
set2$test <- NULL
set2 <- unique(set2)
set_log <- subset(set2, link == "log")
set_linear <- subset(set2, link == "identity")


fitting_times <- sapply(results, function(x) x$fitting_times)
sum(is.na(fitting_times))
which(is.na(fitting_times), arr.ind = TRUE)

## generate plots

# log
mapply(function(n, di, li, cov){
  plotI <- ggplot(subset(settings, obs == n & dim == di & link == li & covariate == cov)) +
    geom_line(aes(x = param_value, y = power, col = test), lwd = 1.2) + xlim(c(-0.25, 0.25)) + ylim(c(0, 1)) + 
    geom_hline(yintercept = 0.05) + theme_bw() + 
    labs(x = "Parameter Value", y = "Empirical Power", col = "Tested Parameter") +
    scale_color_manual(labels = c("alpha" = expression(alpha), "beta" = expression(beta), "covariate" = expression(gamma)),
                       values = c("blue", "red", "yellow")) + 
    theme(legend.position = "bottom",
          text = element_text(size = 30))
  ggsave(paste0("plots/", copula, "/power/power_", li, "_", n, "_dim_", di,"_cov_", cov, ".pdf"), plotI, width = 12, height = 9)
}, n = set_log$obs, li = set_log$link, di = set_log$dim, cov = TRUE)

# linear
mapply(function(n, di, li, cov){
  plotI <- ggplot(subset(settings, obs == n & dim == di & link == li & covariate == cov)) +
    geom_line(aes(x = param_value, y = power, col = test), lwd = 1.2) + xlim(c(0, 0.5)) + ylim(c(0, 1)) + 
    geom_hline(yintercept = 0.05) + theme_bw() + 
    labs(x = "Parameter Value", y = "Empirical Power", col = "Tested Parameter") +
    scale_color_manual(labels = c("alpha" = expression(alpha), "beta" = expression(beta), "covariate" = expression(gamma)),
                       values = c("blue", "red", "yellow")) + 
    theme(legend.position = "bottom",
          text = element_text(size = 30))
  ggsave(paste0("plots/", copula, "/power/power_", li, "_", n, "_dim_", di,"_cov_", cov, ".pdf"), plotI, width = 12, height = 9)
}, n = set_linear$obs, li = set_linear$link, di = set_linear$dim, cov = TRUE)









