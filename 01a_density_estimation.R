rm(list=ls())

library(Distance)

#### DENSITY ESTIMATION - FRANKLIN BAY 2018 ####

# Loading data
data <- read.csv("data/caracaras_2018.csv", sep = ";")
head(data)
hist(data$distance) #relative frequency

# Fitting detection function

# halfnormal, cosine adjustments
halfnorm.pau <- ds(data, key="hn", adjustment="cos", truncation = 140)
plot(halfnorm.pau)

#uniform, cosine adjustments
uniform.pau <- ds(data, key="unif", adjustment="cos", mono="strict", truncation = 140)
plot(uniform.pau)

#hazard-rate, simple polynomial adjustments
hazard.pau <- ds(data, key="hr",  adjustment="poly", truncation = 140)
plot(hazard.pau)

#models fit
#Distance sampling Cramer-von Mises test (unweighted)
#If non-significant good of fit, plausible model

fit.test <- ddf.gof(uniform.pau$ddf)
caracara_qq <- gof_ds(uniform.pau)
caracara_qq

fit.test <- ddf.gof(halfnorm.pau$ddf)
caracara_qq <- gof_ds(halfnorm.pau)
caracara_qq

fit.test <- ddf.gof(hazard.pau$ddf)
caracara_qq <- gof_ds(hazard.pau)
caracara_qq

# All models are statistically plausible

summarize_ds_models(halfnorm.pau, hazard.pau, uniform.pau)
#Hazard-rate discarded, as delta AIC is higher than 2
#only valid if models have the same truncation level!

#Population inference ####

826000/4035150 #area covered 

summary(uniform.pau)

summary(halfnorm.pau)


# % of region covered by the sampling (see 01b_area_density_estimation)
784000/2900000

# para estimar la abundancia de cada grupo uso las abundancias relativas

table(data$age)

51/59 # australis
28/59 # adults
6/59 #immatures
17/59 #juveniles
8/59 # Caracara plancus




