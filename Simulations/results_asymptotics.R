library("data.table")
library("ggplot2")
library("xtable")
library("scales")

# copula <- "Clayton"
copula <- "Frank"


### Asymptotics:

results <- lapply(list.files(paste0(copula, "/asymptotics/"), full.names = TRUE), readRDS)
settings <- lapply(results, "[[", "setting")
settings <- do.call("rbind", settings)

fitting_times <- sapply(results, "[[", "fitting_times")
settings[unique(which(is.na(fitting_times), arr.ind = TRUE)[, 2]),]
settings[which(is.na(colMeans(fitting_times))),]

# Probleme bei Settings mit 1000 Beobachtungen:

results <- results[settings$obs < 1000]
settings <- settings[settings$obs < 1000,]

# Allgemeine Asymptotik (Normalverteilung) über QQ-Plots:

asym_settings <- settings[settings$test == "none" & settings$covariate,]
asym_results <- results[settings$test == "none" & settings$covariate]
asym_results <- lapply(asym_results, function(x) do.call("rbind", lapply(x$param_est, unlist)))

plotIt <- function(index){
  mat <- matrix(1:6, 3)
  labels <- list(expression(delta[0]), expression(beta["0,0"]), expression(beta["0,1"]), expression(alpha["0,0"]), expression(alpha["0,1"]), expression(gamma))
  layout(mat)
  for(i in 1:6){
    qqnorm(asym_results[[index]][, i], main = labels[[i]])
    qqline(asym_results[[index]][, i])
  }
}

for(i in 1:36){
  pdf(file = paste0("plots/", copula, "/asymptotics/qq_dim_", asym_settings$dim[i], 
                    "_link_", asym_settings$link[i], 
                    "_obs_", asym_settings$obs[i], ".pdf"),
      width = 8, height = 12)
  plotIt(i)
  dev.off()
  Sys.sleep(0.5)
  
}

#
mean(asym_results[[34]][, 5] < 1e-5)

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
  x <- reshape(size, direction = "wide",
          idvar = c("dim", "link"), timevar = "obs", # mandatory
          v.names = c("size"), varying = c("50", "100", "250", "400", "500", "750"))
  subset(x, select = names(x) != "400")
}


covariate_results <- empirical_size(settings$test == "covariate")
alpha_results_covariate <- empirical_size(settings$test == "alpha" & settings$covariate)
beta_results_covariate <- empirical_size(settings$test == "beta" & settings$covariate)

alpha_results_no_covariate <- empirical_size(settings$test == "alpha" & !settings$covariate)
beta_results_no_covariate <- empirical_size(settings$test == "beta" & !settings$covariate)
