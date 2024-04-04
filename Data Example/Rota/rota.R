library("glmSTARMA")
library("rjson")
library("spdep")

## Erste beiden Ziffern geben Bundesland an.
## 05 ist NRW


shape <- read_sf("germany_county_shapes.json")
nb <- poly2nb(shape, row.names = shape$RKI_ID)
W1 <- nb2mat(nb, style="W")
row.names(W1) <- colnames(W1) <- shape$RKI_ID

bl <- substr(shape$RKI_ID, 1, 2)
W2 <- matrix(0, nrow = 412, ncol = 412)
for(i in 1:412){
  W2[i, bl == bl[i]] <- 1
}
diag(W2) <- 0
W2 <- W2 / ifelse(colSums(W2) == 0, 1, colSums(W2))
row.names(W2) <- colnames(W2) <- shape$RKI_ID


## Read in features:
ddr_feature <- read.csv("region_features.csv")
ddr_feature$RKI_ID <- ifelse(nchar(as.character(ddr_feature$RKI_ID)) == 4, paste0("0", ddr_feature$RKI_ID), as.character(ddr_feature$RKI_ID))
ddr_feature <- ddr_feature[match(shape$RKI_ID, ddr_feature$RKI_ID),]
rownames(ddr_feature) <- ddr_feature$RKI_ID
ddr_feature <- ddr_feature$is_east

population <- read.csv("germany_population_data.csv")
population <- aggregate(population ~ county + year, data = population, FUN = sum)
population <- reshape(population, timevar = "year", idvar = "county", direction = "wide")
population <- population[match(shape$RKI_NameDE, population$county),]
population <- population$population.2001



## Rota Fit
rota <- read.csv2("rotavirus.csv")
rota <- t(as.matrix(rota[-1]))
rota <- rota[-413,]
rownames(rota) <- gsub("X", "", rownames(rota))
rota <- rota[match(shape$RKI_ID, row.names(rota)),]

means <- log(rowMeans(rota))

pdf("incidence_rota.pdf", width = 8, height = 12)
x <- shape
x$LogIncidence <- means - log(population / 100000)
x$RKI_NameDE <- NULL
x$RKI_NameEN <- NULL
x$RKI_ID <- NULL
plot(x, main = NULL)
dev.off()

pdf("rota_cases.pdf", width = 10, height = 12)
par(mar = c(5, 5, 4, 2) + 0.1, cex.axis=1.5,cex.lab=1.5,cex=1.5)
colsum_rota <- ts(colSums(rota), start = c(2001, 1), frequency = 365.25/7)
plot(colsum_rota, type = "o", 
     ylab = "Weekly number of rota cases", xlab = "Date", 
     pch = 16)
dev.off()

covariates_log <- list()
covariates_log$is_east <- TimeConstant(ddr_feature)
covariates_log$population <- TimeConstant(log(population / 100000))
covariates_log$season_cos <- SpatialConstant(cos(2 * pi / 52 * seq(ncol(rota))))
covariates_log$season_sin <- SpatialConstant(sin(2 * pi / 52 * seq(ncol(rota))))


qics <- sapply(1:8, function(i) glmstarma(rota, model = list(past_obs = rep(1, i)), wlist = list(diag(412), W1),
                                          covariates = covariates_log, 
                                          family = vpoisson("log"), control = list(method = "nloptr", maxit = 1000))$qic)

qics_pm <- sapply(1:8, function(i) glmstarma(rota, model = list(past_obs = rep(1, i), past_mean = 1), wlist = list(diag(412), W1),
                                          covariates = covariates_log, 
                                          family = vpoisson("log"), control = list(method = "nloptr", maxit = 1000))$qic)


model_log_ar <- glmstarma(rota, model = list(past_obs = rep(1, 4)), wlist = list(diag(412), W1),
                       covariates = covariates_log, 
                       family = vpoisson("log"), control = list(method = "nloptr", maxit = 1000))


summary(model_log_ar)


## Linear model for comparison

covariates_linear <- list()
covariates_linear$population <- TimeConstant(population / 100000)
covariates_linear$is_east <- TimeConstant(ddr_feature)


# Fitting only possible without constraints to parameters
model_linear_ar <- glmstarma(rota, model = list(past_obs = c(1), past_mean = 1), wlist = list(diag(412), W1),
                          covariates = covariates_linear, 
                          family = vpoisson("identity"), control = list(method = "optim", maxit = 1000))

