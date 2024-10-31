library("surveillance")
library("glmSTARMA")
# library("rjson")
library("spdep")
library("tictoc")

## First two numbers give the state
## 05 is NRW

shape <- read_sf("../data/germany_county_shapes.json")
nb <- poly2nb(shape, row.names = shape$RKI_ID)
W1 <- nb2mat(nb, style="W")
row.names(W1) <- colnames(W1) <- shape$RKI_ID

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

covariates_log <- list()
covariates_log$is_east <- TimeConstant(ddr_feature)
covariates_log$population <- TimeConstant(log(population / 100000))
covariates_log$season_cos <- SpatialConstant(cos(2 * pi / 52 * seq(ncol(rota))))
covariates_log$season_sin <- SpatialConstant(sin(2 * pi / 52 * seq(ncol(rota))))


rota_sts <- sts(t(rota), population = population, neighbourhood = 1 * (W1 > 0))

# With covariates in 1 component

model_components_endemic <- list(
  end = list(f = addSeason2formula(~1 + population(rota_sts) + ddr + vaccine, period = 52)),
  ar = list(f = ~1),
  ne = list(f = ~1),
  data = list(ddr = t(replicate(ncol(rota), ddr_feature)), 
              vaccine = replicate(412, 1 * ((rota_sts@epoch - min(rota_sts@epoch)) > 654)),
              t = rota_sts@epoch - min(rota_sts@epoch))
)

model_components_auto <- list(
  end = list(f = ~ 1),
  ar = list(f = addSeason2formula(~1 + population(rota_sts) + ddr + vaccine, period=52)),
  ne = list(f = ~ 1),
  data = list(ddr = t(replicate(ncol(rota), ddr_feature)),
              vaccine = replicate(412, 1 * ((rota_sts@epoch - min(rota_sts@epoch)) > 654)),
              t = rota_sts@epoch - min(rota_sts@epoch))
)

model_components_neighbor <- list(
  end = list(f = ~1),
  ar = list(f = ~1),
  ne = list(f = addSeason2formula(~1 + population(rota_sts) + ddr + vaccine, period=52)),
  data = list(ddr = t(replicate(ncol(rota), ddr_feature)),
              vaccine = replicate(412, 1 * ((rota_sts@epoch - min(rota_sts@epoch)) > 654)),
              t = rota_sts@epoch - min(rota_sts@epoch))
)

# Covariates in 2 Components:
model_components_endemic_auto_neighbor <- list(
  end = list(f = addSeason2formula(~1 + population(rota_sts) + ddr + vaccine, period = 52)),
  ar = list(f = addSeason2formula(~1 + population(rota_sts) + ddr + vaccine, period = 52)),
  ne = list(f = ~1),
  data = list(ddr = t(replicate(ncol(rota), ddr_feature)),
              vaccine = replicate(412, 1 * ((rota_sts@epoch - min(rota_sts@epoch)) > 654)),
              t = rota_sts@epoch - min(rota_sts@epoch))
)

model_components_auto_neighbor <- list(
  end = list(f = ~ 1),
  ar = list(f = addSeason2formula(~1 + population(rota_sts) + ddr + vaccine, period=52)),
  ne = list(f = addSeason2formula(~1 + population(rota_sts) + ddr + vaccine, period = 52)),
  data = list(ddr = t(replicate(ncol(rota), ddr_feature)), 
              vaccine = replicate(412, 1 * ((rota_sts@epoch - min(rota_sts@epoch)) > 654)),
              t = rota_sts@epoch - min(rota_sts@epoch))
)

model_components_endemic_neighbor <- list(
  end = list(f = addSeason2formula(~1 + population(rota_sts) + ddr + vaccine, period = 52)),
  ar = list(f = ~1),
  ne = list(f = addSeason2formula(~1 + population(rota_sts) + ddr + vaccine, period=52)),
  data = list(ddr = t(replicate(ncol(rota), ddr_feature)), 
              vaccine = replicate(412, 1 * ((rota_sts@epoch - min(rota_sts@epoch)) > 654)),
              t = rota_sts@epoch - min(rota_sts@epoch))
)


# In 3 Components
model_components_all <- list(
  end = list(f = addSeason2formula(~1 + population(rota_sts) + ddr + vaccine, period = 52)),
  ar = list(f = addSeason2formula(~1 + population(rota_sts) + ddr + vaccine, period=52)),
  ne = list(f = addSeason2formula(~1 + population(rota_sts) + ddr + vaccine, period=52)),
  data = list(ddr = t(replicate(ncol(rota), ddr_feature)), 
              vaccine = replicate(412, 1 * ((rota_sts@epoch - min(rota_sts@epoch)) > 654)),
              t = rota_sts@epoch - min(rota_sts@epoch))
)


# With time_lag 2
# With covariates in 1 component

model_components_endemic_2 <- list(
  end = list(f = addSeason2formula(~1 + population(rota_sts) + ddr + vaccine, period = 52)),
  ar = list(f = ~1, lag = 2),
  ne = list(f = ~1, lag = 2),
  subset = 3:903,
  data = list(ddr = t(replicate(ncol(rota), ddr_feature)), 
              vaccine = replicate(412, 1 * ((rota_sts@epoch - min(rota_sts@epoch)) > 654)),
              t = rota_sts@epoch - min(rota_sts@epoch))
)

model_components_auto_2 <- list(
  end = list(f = ~ 1),
  ar = list(f = addSeason2formula(~1 + population(rota_sts) + ddr + vaccine, period=52), lag = 2),
  ne = list(f = ~ 1, lag = 2),
  subset = 3:903,
  data = list(ddr = t(replicate(ncol(rota), ddr_feature)), 
              vaccine = replicate(412, 1 * ((rota_sts@epoch - min(rota_sts@epoch)) > 654)),
              t = rota_sts@epoch - min(rota_sts@epoch))
)

model_components_neighbor_2 <- list(
  end = list(f = ~1),
  ar = list(f = ~1, lag = 2),
  ne = list(f = addSeason2formula(~1 + population(rota_sts) + ddr + vaccine, period=52), lag = 2),
  subset = 3:903,
  data = list(ddr = t(replicate(ncol(rota), ddr_feature)), 
              vaccine = replicate(412, 1 * ((rota_sts@epoch - min(rota_sts@epoch)) > 654)),
              t = rota_sts@epoch - min(rota_sts@epoch))
)

# Covariates in 2 Components:
model_components_endemic_auto_2 <- list(
  end = list(f = addSeason2formula(~1 + population(rota_sts) + ddr + vaccine, period = 52)),
  ar = list(f = addSeason2formula(~1 + population(rota_sts) + ddr + vaccine, period = 52), lag = 2),
  ne = list(f = ~1, lag = 2),
  subset = 3:903,
  data = list(ddr = t(replicate(ncol(rota), ddr_feature)), 
              vaccine = replicate(412, 1 * ((rota_sts@epoch - min(rota_sts@epoch)) > 654)),
              t = rota_sts@epoch - min(rota_sts@epoch))
)

model_components_auto_neighbor_2 <- list(
  end = list(f = ~ 1),
  ar = list(f = addSeason2formula(~1 + population(rota_sts) + ddr + vaccine, period=52), lag = 2),
  ne = list(f = addSeason2formula(~1 + population(rota_sts) + ddr + vaccine, period = 52), lag = 2),
  subset = 3:903,
  data = list(ddr = t(replicate(ncol(rota), ddr_feature)), 
              vaccine = replicate(412, 1 * ((rota_sts@epoch - min(rota_sts@epoch)) > 654)),
              t = rota_sts@epoch - min(rota_sts@epoch))
)

model_components_endemic_neighbor_2 <- list(
  end = list(f = addSeason2formula(~1 + population(rota_sts) + ddr + vaccine, period = 52)),
  ar = list(f = ~1, lag = 2),
  ne = list(f = addSeason2formula(~1 + population(rota_sts) + ddr + vaccine, period=52), lag = 2),
  subset = 3:903,
  data = list(ddr = t(replicate(ncol(rota), ddr_feature)), 
              vaccine = replicate(412, 1 * ((rota_sts@epoch - min(rota_sts@epoch)) > 654)),
              t = rota_sts@epoch - min(rota_sts@epoch))
)


# In 3 Components
model_components_all_2 <- list(
  end = list(f = addSeason2formula(~1 + population(rota_sts) + ddr + vaccine, period = 52)),
  ar = list(f = addSeason2formula(~1 + population(rota_sts) + ddr + vaccine, period=52), lag = 2),
  ne = list(f = addSeason2formula(~1 + population(rota_sts) + ddr + vaccine, period=52), lag = 2),
  subset = 3:903,
  data = list(ddr = t(replicate(ncol(rota), ddr_feature)), 
              vaccine = replicate(412, 1 * ((rota_sts@epoch - min(rota_sts@epoch)) > 654)),
              t = rota_sts@epoch - min(rota_sts@epoch))
)



to_fit <- ls()
to_fit <- to_fit[grepl("model", to_fit)]

for(mod in to_fit){
  mod_to_fit <- get(mod)
  tic()
  fit <- hhh4(rota_sts, control = mod_to_fit)
  time <- toc()
  file_name <- paste0("../results/hhh4/hhh4_rota_", mod, ".rds")
  res <- list(fitting_time = time, model = fit)
  saveRDS(res, file_name)
}


# Fit on subset and predict last year
for(mod in to_fit){
  mod_to_fit <- get(mod)
  if(is.null(mod_to_fit$subset)){
    mod_to_fit$subset <- 2:835
  } else {
    mod_to_fit$subset <- 3:835
  }
  tic()
  fit <- hhh4(rota_sts, control = mod_to_fit)
  time <- toc()
  predictions <- predict(fit, newSubset = 836:903, type = "response")
  truth <- t(rota[,836:903])
  file_name <- paste0("../results/hhh4/predict_hhh4_rota", mod, ".rds")
  res <- list(fitting_time = time, model = fit, prediction = predictions,
              truth = truth)
  saveRDS(res, file_name)
}

