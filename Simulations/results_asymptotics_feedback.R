library("data.table")
library("ggplot2")
library("xtable")
library("scales")

copula <- "Clayton"
#copula <- "Frank"


### Asymptotics:

results <- lapply(list.files(paste0(copula, "/asymptotics_feedback/"), full.names = TRUE), readRDS)
settings <- lapply(results, "[[", "setting")
settings <- do.call("rbind", settings)


# Remove non-converged results

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



subset(settings, errors > 0)



settings$time <- sapply(sapply(results, "[[", "fitting_times"), mean)




# Funktion um empirical size zu extrahieren:

empirical_size <- function(constraint){
  wald_result <- lapply(results, "[[", "wald_result")
  
  sub_settings <- settings[constraint,]
  sub_results <- wald_result[constraint]
  sub_results <- lapply(sub_results, function(x) sapply(x, "[[", "p_value"))
  sub_results <- lapply(sub_results, unlist)
  
  size <- sapply(sub_results, function(x) mean(x < 0.05))
  
  size <- data.frame(size = size, 
                     obs = sub_settings$obs,
                     dim = sub_settings$dim,
                     link = sub_settings$link)
  size <- size[order(size$obs, size$link, size$dim),]
  x <- reshape(size, direction = "wide",
               idvar = c("dim", "link"), timevar = "obs", # mandatory
               v.names = c("size"), varying = c("5", "10", "20", "50", "100"))
  subset(x, select = names(x) != "400")
}





covariate_results <- empirical_size(settings$test == "covariate")
alpha_results_covariate <- empirical_size(settings$test == "alpha" & settings$covariate)
beta_results_covariate <- empirical_size(settings$test == "beta" & settings$covariate)







settings <- settings[order(settings$link, settings$test, settings$dim, settings$obs),]
settings$weighted_time <- settings$time * (1 - settings$errors)
fitting <- aggregate(weighted_time ~ link + dim + obs, data = settings, FUN = sum)
weight <- aggregate(errors ~ link + dim + obs, data = settings, FUN = function(x) sum(1 - x))
fitting$weighted_time <- fitting$weighted_time / weight$errors

fitting <- fitting[order(fitting$link, fitting$dim, fitting$obs),]

fitting <- reshape(fitting, direction = "wide",
              idvar = c("dim", "link"), timevar = "obs", # mandatory
              v.names = c("weighted_time"), varying = c("5", "10", "20", "50", "100"))







