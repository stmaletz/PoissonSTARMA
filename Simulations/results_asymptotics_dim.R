library("data.table")
library("ggplot2")
library("xtable")
library("scales")

copula <- "Clayton"
#copula <- "Frank"


### Asymptotics:

results <- lapply(list.files(paste0(copula, "/asymptotics_dim/"), full.names = TRUE), readRDS)
settings <- lapply(results, "[[", "setting")
settings <- do.call("rbind", settings)

fitting_times <- sapply(results, "[[", "fitting_times")
settings[unique(which(is.na(fitting_times), arr.ind = TRUE)[, 2]),]
settings[which(is.na(colMeans(fitting_times))),]


params <- lapply(results, "[[", "param_est")
params <- lapply(params, function(x) lapply(x, unlist))
for(i in seq_along(params)){
  if(results[[i]]$setting$link == "log"){
    params[[i]] <- which(sapply(params[[i]], function(y) length(y) == 1 | y[4] == 0))
  } else {
    params[[i]] <- which(sapply(params[[i]], function(y) length(y) == 1 | y[4] == 1))
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





# Get Parameter estimates

extract_params <- function(elem){
  params <- sapply(elem$param_est, unlist)
  rowMeans(params)
}


x <- cbind(settings, t(sapply(results, extract_params)))
x[order(x$link, x$dim, x$test, x$obs),]



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
beta_results_covariate <- empirical_size(settings$test == "beta" & settings$covariate)


settings$fitting_times <- colMeans(fitting_times)

settings <- settings[order(settings$link, settings$test, settings$dim, settings$obs),]


ft <- aggregate(fitting_times ~ link + dim + obs, data = settings, FUN = mean)
ft <- ft[order(ft$link, ft$dim, ft$obs),]

ft <- reshape(ft, direction = "wide",
             idvar = c("dim", "link"), timevar = "obs", # mandatory
             v.names = c("fitting_times"), varying = c("5", "10", "20", "50", "100"))












