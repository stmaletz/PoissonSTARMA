library("mgcv")
library("PNAR")
library("glmSTARMA")
library("spdep")
library("tictoc")
library("surveillance")

data("crime")
data("crime_W")

population <- read.csv("../data/pop.csv")$x
unemployed <- read.csv("../data/unemp.csv")$x
wealth <- read.csv("../data/wealth.csv")$x
young_males <- read.csv("../data/ym.csv")$x / population
temperature <- read.csv2("../data/temperature.csv")$Temperatur
population <- log(population)

crime <- t(crime)

crime_sts <- sts(t(crime), start = c(2010, 1), frequency = 12, 
                 population = population,
                 neighbourhood = 1 * (crime_W > 0))


# Models with time lag 1:

# Mit Kovariablen in lediglich einer Komponente
model_components_endemic <- list(
  end = list(f = addSeason2formula(~ 1 + t + population(crime_sts) + unemployed + wealth + young_males + temperature, period = 12)),
  ar = list(f = ~ 1),
  ne = list(f = ~ 1),
  data = list(unemployed = t(replicate(ncol(crime), unemployed)), 
              wealth = t(replicate(ncol(crime), wealth)), 
              young_males = t(replicate(ncol(crime), young_males)), 
              temperature = replicate(nrow(crime), temperature),
              t = crime_sts@epoch - min(crime_sts@epoch))
)

model_components_auto <- list(
  end = list(f = ~ 1),
  ar = list(f = addSeason2formula(~ 1 + t + population(crime_sts) + unemployed + wealth + young_males + temperature, period = 12)),
  ne = list(f = ~ 1),
  data = list(unemployed = t(replicate(ncol(crime), unemployed)), 
              wealth = t(replicate(ncol(crime), wealth)), 
              young_males = t(replicate(ncol(crime), young_males)), 
              temperature = replicate(nrow(crime), temperature),
              t = crime_sts@epoch - min(crime_sts@epoch))
)

model_components_neighbor <- list(
  end = list(f = ~ 1),
  ar = list(f = ~ 1),
  ne = list(f = addSeason2formula(~ 1 + t + population(crime_sts) + unemployed + wealth + young_males + temperature, period = 12)),
  data = list(unemployed = t(replicate(ncol(crime), unemployed)), 
              wealth = t(replicate(ncol(crime), wealth)), 
              young_males = t(replicate(ncol(crime), young_males)), 
              temperature = replicate(nrow(crime), temperature),
              t = crime_sts@epoch - min(crime_sts@epoch))
)

# Covariates in 2 components:
model_components_endemic_auto <- list(
  end = list(f = addSeason2formula(~ 1 + t + population(crime_sts) + unemployed + wealth + young_males + temperature, period = 12)),
  ar = list(f = addSeason2formula(~ 1 + t + population(crime_sts) + unemployed + wealth + young_males + temperature, period = 12)),
  ne = list(f = ~ 1),
  data = list(unemployed = t(replicate(ncol(crime), unemployed)), 
              wealth = t(replicate(ncol(crime), wealth)), 
              young_males = t(replicate(ncol(crime), young_males)), 
              temperature = replicate(nrow(crime), temperature),
              t = crime_sts@epoch - min(crime_sts@epoch))
)

model_components_endemic_neighbor <- list(
  end = list(f = addSeason2formula(~ 1 + t + population(crime_sts) + unemployed + wealth + young_males + temperature, period = 12)),
  ar = list(f = ~ 1),
  ne = list(f = addSeason2formula(~ 1 + t + population(crime_sts) + unemployed + wealth + young_males + temperature, period = 12)),
  data = list(unemployed = t(replicate(ncol(crime), unemployed)), 
              wealth = t(replicate(ncol(crime), wealth)), 
              young_males = t(replicate(ncol(crime), young_males)), 
              temperature = replicate(nrow(crime), temperature),
              t = crime_sts@epoch - min(crime_sts@epoch))
)

model_components_auto_neighbor <- list(
  end = list(f = ~ 1),
  ar = list(f = addSeason2formula(~ 1 + t + population(crime_sts) + unemployed + wealth + young_males + temperature, period = 12)),
  ne = list(f = addSeason2formula(~ 1 + t + population(crime_sts) + unemployed + wealth + young_males + temperature, period = 12)),
  data = list(unemployed = t(replicate(ncol(crime), unemployed)), 
              wealth = t(replicate(ncol(crime), wealth)), 
              young_males = t(replicate(ncol(crime), young_males)), 
              temperature = replicate(nrow(crime), temperature),
              t = crime_sts@epoch - min(crime_sts@epoch))
)


# Covariates everywhere:
model_components_all <- list(
  end = list(f = addSeason2formula(~ 1 + t + population(crime_sts) + unemployed + wealth + young_males + temperature, period = 12)),
  ar = list(f = addSeason2formula(~ 1 + t + population(crime_sts) + unemployed + wealth + young_males + temperature, period = 12)),
  ne = list(f = addSeason2formula(~ 1 + t + population(crime_sts) + unemployed + wealth + young_males + temperature, period = 12)),
  data = list(unemployed = t(replicate(ncol(crime), unemployed)), 
              wealth = t(replicate(ncol(crime), wealth)), 
              young_males = t(replicate(ncol(crime), young_males)), 
              temperature = replicate(nrow(crime), temperature),
              t = crime_sts@epoch - min(crime_sts@epoch))
)


# Timelag 2:

# 1 component
model_components_endemic_2 <- list(
  end = list(f = addSeason2formula(~ 1 + t + population(crime_sts) + unemployed + wealth + young_males + temperature, period = 12)),
  ar = list(f = ~ 1, lag = 2),
  ne = list(f = ~ 1, lag = 2),
  subset = 3:72,
  data = list(unemployed = t(replicate(ncol(crime), unemployed)), 
              wealth = t(replicate(ncol(crime), wealth)), 
              young_males = t(replicate(ncol(crime), young_males)), 
              temperature = replicate(nrow(crime), temperature),
              t = crime_sts@epoch - min(crime_sts@epoch))
)

model_components_auto_2 <- list(
  end = list(f = ~ 1),
  ar = list(f = addSeason2formula(~ 1 + t + population(crime_sts) + unemployed + wealth + young_males + temperature, period = 12), lag = 2),
  ne = list(f = ~ 1, lag = 2),
  subset = 3:72,
  data = list(unemployed = t(replicate(ncol(crime), unemployed)), 
              wealth = t(replicate(ncol(crime), wealth)), 
              young_males = t(replicate(ncol(crime), young_males)), 
              temperature = replicate(nrow(crime), temperature),
              t = crime_sts@epoch - min(crime_sts@epoch))
)

model_components_neighbor_2 <- list(
  end = list(f = ~ 1),
  ar = list(f = ~ 1, lag = 2),
  ne = list(f = addSeason2formula(~ 1 + t + population(crime_sts) + unemployed + wealth + young_males + temperature, period = 12), lag = 2),
  subset = 3:72,
  data = list(unemployed = t(replicate(ncol(crime), unemployed)), 
              wealth = t(replicate(ncol(crime), wealth)), 
              young_males = t(replicate(ncol(crime), young_males)), 
              temperature = replicate(nrow(crime), temperature),
              t = crime_sts@epoch - min(crime_sts@epoch))
)

# 2:
model_components_endemic_auto_2 <- list(
  end = list(f = addSeason2formula(~ 1 + t + population(crime_sts) + unemployed + wealth + young_males + temperature, period = 12)),
  ar = list(f = addSeason2formula(~ 1 + t + population(crime_sts) + unemployed + wealth + young_males + temperature, period = 12), lag = 2),
  ne = list(f = ~ 1, lag = 2),
  subset = 3:72,
  data = list(unemployed = t(replicate(ncol(crime), unemployed)), 
              wealth = t(replicate(ncol(crime), wealth)), 
              young_males = t(replicate(ncol(crime), young_males)), 
              temperature = replicate(nrow(crime), temperature),
              t = crime_sts@epoch - min(crime_sts@epoch))
)

model_components_endemic_neighbor_2 <- list(
  end = list(f = addSeason2formula(~ 1 + t + population(crime_sts) + unemployed + wealth + young_males + temperature, period = 12)),
  ar = list(f = ~ 1, lag = 2),
  ne = list(f = addSeason2formula(~ 1 + t + population(crime_sts) + unemployed + wealth + young_males + temperature, period = 12), lag = 2),
  subset = 3:72,
  data = list(unemployed = t(replicate(ncol(crime), unemployed)), 
              wealth = t(replicate(ncol(crime), wealth)), 
              young_males = t(replicate(ncol(crime), young_males)), 
              temperature = replicate(nrow(crime), temperature),
              t = crime_sts@epoch - min(crime_sts@epoch))
)

model_components_auto_neighbor_2 <- list(
  end = list(f = ~ 1),
  ar = list(f = addSeason2formula(~ 1 + t + population(crime_sts) + unemployed + wealth + young_males + temperature, period = 12), lag = 2),
  ne = list(f = addSeason2formula(~ 1 + t + population(crime_sts) + unemployed + wealth + young_males + temperature, period = 12), lag = 2),
  subset = 3:72,
  data = list(unemployed = t(replicate(ncol(crime), unemployed)), 
              wealth = t(replicate(ncol(crime), wealth)), 
              young_males = t(replicate(ncol(crime), young_males)), 
              temperature = replicate(nrow(crime), temperature),
              t = crime_sts@epoch - min(crime_sts@epoch))
)

# everywhere:
model_components_all_2 <- list(
  end = list(f = addSeason2formula(~ 1 + t + population(crime_sts) + unemployed + wealth + young_males + temperature, period = 12)),
  ar = list(f = addSeason2formula(~ 1 + t + population(crime_sts) + unemployed + wealth + young_males + temperature, period = 12),
            lag = 2),
  ne = list(f = addSeason2formula(~ 1 + t + population(crime_sts) + unemployed + wealth + young_males + temperature, period = 12), 
            lag = 2),
  subset = 3:72,
  data = list(unemployed = t(replicate(ncol(crime), unemployed)), 
              wealth = t(replicate(ncol(crime), wealth)), 
              young_males = t(replicate(ncol(crime), young_males)), 
              temperature = replicate(nrow(crime), temperature),
              t = crime_sts@epoch - min(crime_sts@epoch))
)



## Model fitting

to_fit <- ls()
to_fit <- to_fit[grepl("model", to_fit)]

for(mod in to_fit){
  mod_to_fit <- get(mod)
  tic()
  fit <- hhh4(crime_sts, control = mod_to_fit)
  time <- toc()
  file_name <- paste0("../results/hhh4/hhh4_crime_", mod, ".rds")
  res <- list(fitting_time = time, model = fit)
  saveRDS(res, file_name)
}


# Fit on subset and predict last year

for(mod in to_fit){
  mod_to_fit <- get(mod)
  if(is.null(mod_to_fit$subset)){
    mod_to_fit$subset <- 2:60
  } else {
    mod_to_fit$subset <- 3:60
  }
  tic()
  fit <- hhh4(crime_sts, control = mod_to_fit)
  time <- toc()
  predictions <- predict(fit, newSubset = 61:72, type = "response")
  truth <- t(crime[,61:72])
  file_name <- paste0("../results/hhh4/predict_hhh4_crime_", mod, ".rds")
  res <- list(fitting_time = time, model = fit, prediction = predictions,
              truth = truth)
  saveRDS(res, file_name)
}

