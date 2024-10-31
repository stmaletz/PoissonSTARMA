library("glmSTARMA")
library("rjson")
library("spdep")
library("tictoc")

## First two numbers give the state
## 05 is NRW

shape <- read_sf("data/germany_county_shapes.json")
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
ddr_feature <- read.csv("data/region_features.csv")
ddr_feature$RKI_ID <- ifelse(nchar(as.character(ddr_feature$RKI_ID)) == 4, paste0("0", ddr_feature$RKI_ID), as.character(ddr_feature$RKI_ID))
ddr_feature <- ddr_feature[match(shape$RKI_ID, ddr_feature$RKI_ID),]
rownames(ddr_feature) <- ddr_feature$RKI_ID
ddr_feature <- ddr_feature$is_east

population <- read.csv("data/germany_population_data.csv")
population <- aggregate(population ~ county + year, data = population, FUN = sum)
population <- reshape(population, timevar = "year", idvar = "county", direction = "wide")
population <- population[match(shape$RKI_NameDE, population$county),]
population <- population$population.2001


## Rota Fit
rota <- read.csv2("data/rotavirus.csv")
rota <- t(as.matrix(rota[-1]))
rota <- rota[-413,]
rownames(rota) <- gsub("X", "", rownames(rota))
rota <- rota[match(shape$RKI_ID, row.names(rota)),]

means <- log(rowMeans(rota))

pdf("plots/incidence_rota.pdf", width = 8, height = 12)
# postscript("plots/incidence_rota.ps", width = 8, height = 12)
x <- shape
x$LogIncidence <- means - log(population / 100000)
x$RKI_NameDE <- NULL
x$RKI_NameEN <- NULL
x$RKI_ID <- NULL
plot(x, main = NULL)
dev.off()

pdf("plots/rota_cases.pdf", width = 10, height = 12)
# postscript("plots/rota_cases.ps", width = 10, height = 12)
par(mar = c(5, 5, 4, 2) + 0.1, cex.axis=1.5,cex.lab=1.5,cex=1.5)
colsum_rota <- ts(colSums(rota), start = c(2001, 1), frequency = 365.25/7)
plot(colsum_rota, type = "o", 
     ylab = "Weekly number of rota cases", xlab = "Date", 
     pch = 16)
#rect(xleft=2013.5, ybottom=0, xright=2018.75, ytop=max(colsum_rota) + 100, col=rgb(0.5, 0.5, 0.5, 0.5), border=NA)
rect(xleft = 2013.5, ybottom = 0, xright = 2018.75, ytop = max(colsum_rota) + 1, col = rgb(0.5, 0.5, 0.5), density = 30, angle = 45)
points(colsum_rota, type = "o", pch = 16)
dev.off()


# Descriptive statistics
mean(rota)
max(rota)

# Number of neihbors
summary(rowSums(W1 > 0))
# Maximum at 12 Neighbors


sum(rowSums(W1 > 0) == 1)
## In total 27 enclaves (only one neighbor)


## extract state
state <- substr(rownames(rota), 1, 2)
levels <- unique(state)
sapply(levels, function(lvl) mean(rota[state == lvl,] / (population[state == lvl] / 100000)) )



