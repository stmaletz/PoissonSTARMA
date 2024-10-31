library("mgcv")
library("PNAR")
library("glmSTARMA")
library("spdep")
library("tictoc")

data("crime")
data("crime_W")

index <- as.numeric(Sys.getenv("PBS_ARRAYID"))

population <- read.csv("../data/pop.csv")$x
unemployed <- read.csv("../data/unemp.csv")$x
wealth <- read.csv("../data/wealth.csv")$x
young_males <- read.csv("../data/ym.csv")$x
temperature <- read.csv2("../data/temperature.csv")$Temperatur

crime <- t(crime)


crime_df <- data.frame(crimes_count = c(crime),
           location = factor(rep(paste0("X", 1:552), ncol(crime)), levels = paste0("X", 1:552)),
           population = rep(log(population / 1000), ncol(crime)),
           unemployed = rep(unemployed, ncol(crime)),
           wealth = rep(wealth, ncol(crime)),
           young_males = rep(young_males / population, ncol(crime)),
           temperature = rep(temperature, each = 552),
           trend = rep(1:72, each = 552))


colnames(crime_W) <- paste0("X", 1:552)
row.names(crime_W) <- paste0("X", 1:552)

wlist <- mat2listw(crime_W, style = "W")
x <- lapply(wlist$neighbours, function(y){
  z <- paste0("X", y)
  factor(z, levels = paste0("X",1:552))
})
names(x) <- paste0("X", 1:552)


dfs <- expand.grid(test = c(FALSE, TRUE), trend_covariate = c(TRUE, FALSE),
time = c(5, 10, 15, 20, 25, 30), space = c(5, 10, 15, 20, 25, 30))

df_time <- dfs$time[index]
df_space <- dfs$space[index]
test_train <- dfs$test[index]
trend_cov <- dfs$trend_covariate[index]

# Fit Models without autoregression
if(!test_train){
  if(trend_cov){
      tic()
      gam_model <- gam(crimes_count ~ population + unemployed + wealth + trend +
                         young_males + temperature +
                         te(trend, location,
                            bs = c("cr", "mrf"),
                            k = c(df_time, df_space),
                            xt = list(list(), list(nb = x))),
                       data = crime_df, family = poisson(),
                       method = "REML")
      time <- toc()
  } else {
      tic()
      gam_model <- gam(crimes_count ~ population + unemployed + wealth +
                         young_males + temperature +
                         te(trend, location,
                            bs = c("cr", "mrf"),
                            k = c(df_time, df_space),
                            xt = list(list(), list(nb = x))),
                       data = crime_df, family = poisson(),
                       method = "REML")
      time <- toc()
  }
  time <- time$toc - time$tic
  file_name <- paste0("../results/mgcv/mgcv_crime_trend_", trend_cov, "_space_", df_space, "_time_", df_time, ".rds")
  res <- list(fitting_time = time, model = gam_model)
  saveRDS(res, file_name)
}

## Use 5 years for training and then make predictions for year 6:

crime_df_train <- crime_df[seq(60*552),]
crime_df_train <- crime_df_train[complete.cases(crime_df),]
crime_df_test <- crime_df[-seq(60*552),]

# Fit Models without autoregression
if(test_train){
    if(trend_cov){
        tic()
        gam_model <- gam(crimes_count ~ population + unemployed + wealth + trend + young_males + temperature +
                           te(trend, location,
                              bs = c("cr", "mrf"),
                              k = c(df_time, df_space),
                              xt = list(list(), list(nb = x))),
                         data = crime_df_train, family = poisson(),
                         method = "REML")
        time <- toc()
    } else {
        tic()
        gam_model <- gam(crimes_count ~ population + unemployed + wealth +
                           young_males + temperature +
                           te(trend, location,
                              bs = c("cr", "mrf"),
                              k = c(df_time, df_space),
                              xt = list(list(), list(nb = x))),
                         data = crime_df_train, family = poisson(),
                         method = "REML")
        time <- toc()
    }
    time <- time$toc - time$tic
    file_name <- paste0("../results/mgcv/predict_mgcv_crime_trend_", trend_cov, "_space_", df_space, "_time_", df_time, ".rds")
  predictions <- predict(gam_model, newdata = crime_df_test, type = "response")
  truth <- crime_df_test$crimes_count
  res <- list(fitting_time = time, model = gam_model, prediction = predictions,
              truth = truth)
  saveRDS(res, file_name)
}

