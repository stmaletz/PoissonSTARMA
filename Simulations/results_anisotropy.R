library("data.table")
library("ggplot2")
library("xtable")
library("scales")

copula <- "Clayton"
# copula <- "Frank"

# Comments are regarding clayton copula


## Anisotropy
results <- lapply(list.files(paste0(copula, "/anisotropy/"), full.names = TRUE), readRDS)

fitting_times_isotropy <- sapply(results, "[[", "fitting_times_isotropy")
fitting_times_anisotropy <- sapply(results, "[[", "fitting_times_anisotropy")
sum(is.na(fitting_times_isotropy))
sum(is.na(fitting_times_anisotropy))
# No Problems


settings <- rbindlist(lapply(results, "[[", "setting"))


# How often is the anisotropy detected as anisotropy?

# Correct results for missing constant
wald <- lapply(results, function(x) data.frame(past_mean = sapply(x$wald_results, function(y) (x$setting$obs - 1) * y$past_mean$statistics),
                                               past_obs = sapply(x$wald_results, function(y) (x$setting$obs - 1) * y$past_obs$statistics)))

settings$past_obs_detect <- sapply(wald, function(x) mean(x$past_obs > qchisq(0.95, 1)))
settings$past_mean_detect <- sapply(wald, function(x) mean(x$past_mean > qchisq(0.95, 1)))

xtable(
  reshape(subset(settings, select = c(link, obs, past_obs_detect)), idvar = "link", timevar = "obs", v.names = "past_obs_detect", direction = "wide"),
digits = 3)

xtable(
  reshape(subset(settings, select = c(link, obs, past_mean_detect)), idvar = "link", timevar = "obs", v.names = "past_mean_detect", direction = "wide"),
digits = 3)


# How often is spatial dependence detected by the isotropic model?

wald_iso <- lapply(results, "[[", "wald_results_isotropy")
wald_iso <- lapply(wald_iso, function(x) rbindlist(lapply(x, function(y) data.frame(past_mean = y$past_mean, past_obs = y$past_obs))))

settings$past_obs_iso_test <- sapply(wald_iso, function(x) mean(x$past_obs < 0.05))
settings$past_mean_iso_test <- sapply(wald_iso, function(x) mean(x$past_mean < 0.05))


xtable(
  reshape(subset(settings, select = c(link, obs, past_obs_iso_test), dim == 9), idvar = "link", timevar = "obs", v.names = "past_obs_iso_test", direction = "wide"),
  digits = 3)

xtable(
  reshape(subset(settings, select = c(link, obs, past_mean_iso_test)), idvar = "link", timevar = "obs", v.names = "past_mean_iso_test", direction = "wide"),
  digits = 3)



## How often is the isotropic model preferred over the anisotropic by qic?

qic_isotropy <- lapply(results, "[[", "qic_isotropy")
qic_anisotropy <- lapply(results, "[[", "qic_anisotropy")
settings$iso_preferred <- mapply(function(x, y) mean(x < y), x = qic_isotropy, y = qic_anisotropy)

settings$mse_iso <- sapply(results, function(x) mean(x$residuals_isotropy))
settings$mse_aniso <- sapply(results, function(x) mean(x$residuals_anisotropy))

xtable(
  reshape(subset(settings, select = c(link, obs, iso_preferred)), idvar = "link", timevar = "obs", v.names = "iso_preferred", direction = "wide"),
  digits = 3)



