# Rota
library("mgcv")
# library("PNAR")
library("glmSTARMA")
library("spdep")
library("tictoc")

## First two numbers give the state
## 05 is NRW

index <- as.numeric(Sys.getenv("PBS_ARRAYID"))

shape <- read_sf("../data/germany_county_shapes.json")
nb <- poly2nb(shape, row.names = shape$RKI_ID)
W1 <- nb2mat(nb, style="W", zero.policy = TRUE)
row.names(W1) <- colnames(W1) <- shape$RKI_ID

W2 <- 1 * (W1 > 0)
W2 <- W2 %*% W2
diag(W2) <- 0
W2[W1 > 0] <- 0
W2 <- W2 / colSums(W2)


## Read in features:
ddr_feature <- read.csv("../data/region_features.csv")
ddr_feature$RKI_ID <- ifelse(nchar(as.character(ddr_feature$RKI_ID)) == 4, paste0("0", ddr_feature$RKI_ID), as.character(ddr_feature$RKI_ID))
ddr_feature <- ddr_feature[match(shape$RKI_ID, ddr_feature$RKI_ID),]
rownames(ddr_feature) <- ddr_feature$RKI_ID
ddr_feature <- ddr_feature$is_east

population <- read.csv("../data/germany_population_data.csv")
population <- aggregate(population ~ county + year, data = population, FUN = sum)
population <- reshape(population, timevar = "year", idvar = "county", direction = "wide")
population <- population[match(shape$RKI_NameDE, population$county),]
population <- population$population.2001

## Rota Fit
rota <- read.csv2("../data/rotavirus.csv")
rota <- t(as.matrix(rota[-1]))
rota <- rota[-413,]
rownames(rota) <- gsub("X", "", rownames(rota))
rota <- rota[match(shape$RKI_ID, row.names(rota)),]

rota_lag_1 <- rota[, -ncol(rota)]
rota_lag_1 <- cbind(NA, rota_lag_1)
rota_neigh_lag_1 <- W1 %*% rota_lag_1
rota_neigh_lag_2 <- W2 %*% rota_lag_1

rota_df <- data.frame(rota_count = c(rota),
                      time = c(rep(1:903, each = 412)),
                      location = factor(rep(row.names(rota), ncol(rota)), levels = row.names(rota)),
                      population = rep(log(population / 100000), ncol(rota)),
                      is_east = rep(ddr_feature, ncol(rota)),
                      season_cos = rep(cos(2 * pi / 52 * seq(ncol(rota))), each = nrow(rota)),
                      season_sin = rep(sin(2 * pi / 52 * seq(ncol(rota))), each = nrow(rota)),
                      vaccine = rep(1 * (seq(ncol(rota)) >= 654), each = nrow(rota)),
                      rota_lag = c(rota_lag_1),
                      rota_neigh_lag = c(rota_neigh_lag_1),
                      rota_neigh_lag_2 = c(rota_neigh_lag_2))


wlist <- mat2listw(W1, style = "W")
x <- lapply(wlist$neighbours, function(y){
  z <- names(y)
  factor(z, levels = row.names(rota))
})
names(x) <- row.names(rota)


dfs <- expand.grid(test = c(FALSE, TRUE), time = c(5, 10, 15, 20, 25, 30, 40), space = c(5, 10, 15, 20, 25, 30, 40))
df_time <- dfs$time[index]
df_space <- dfs$space[index]
test_train <- dfs$test[index]



# Fit Models without autoregression
# for(df in dfs){
if(!test_train){
  tic()
  gam_model <- gam(rota_count ~ population + is_east + season_cos + season_sin + vaccine +
                     te(time, location, 
                        bs = c("cr", "mrf"),
                        k = c(df_time, df_space),
                        xt = list(list(), list(nb = x))),
                   data = rota_df, family = poisson(),
                   method = "REML")
  time <- toc()
  time <- time$toc - time$tic
  file_name <- paste0("../results/mgcv/mgcv_rota_no_ar_df_", df_time, "_", df_space, ".rds")
  res <- list(fitting_time = time, model = gam_model)
  saveRDS(res, file_name)
}

## Use 5 years for training and then make predictions for year 6:
# Predict year of 01-2017 until end of data (2018-16)


rota_df_train <- rota_df[seq(835*412),]
rota_df_train <- rota_df_train[complete.cases(rota_df_train),]
rota_df_test <- rota_df[-seq(835*412),]



# Fit Models without autoregression
# for(df in dfs){
if(test_train){
  tic()
  gam_model <- gam(rota_count ~ population + is_east + season_cos + season_sin + vaccine +
                                  te(time, location, 
                                     bs = c("cr", "mrf"),
                                     k = c(df_time, df_space),
                                     xt = list(list(), list(nb = x))),
                                data = rota_df_train, family = poisson(),
                                method = "REML")
  time <- toc()
  time <- time$toc - time$tic
  predictions <- predict(gam_model, newdata = rota_df_test, type = "response")
  truth <- rota_df_test$rota_count
  file_name <- paste0("../results/mgcv/predict_mgcv_rota_no_ar_df_", df_time, "_", df_space, ".rds")
  res <- list(fitting_time = time, model = gam_model, prediction = predictions,
              truth = truth)
  saveRDS(res, file_name)
}
