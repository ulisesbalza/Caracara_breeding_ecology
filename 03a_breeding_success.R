ls()
rm(list=ls())

# required libraries
library(gplots)
library(nlme)
library(MuMIn)
library(Amelia)
library(mlbench)
library(corrplot)
library(lme4)
library(ggplot2)
library(MASS)
library(mgcv)
library(DHARMa)
library(mgcViz)
library(ggeffects)
library(gridExtra)
library(ggpubr)

par(mfrow = c(1, 1))

# correlation among early and late breeding performance ####
data_correlation <- read.csv("data/data_brood_correlation.csv",header=T,sep=";")
head(data_correlation)
cor.test(data_correlation$early, data_correlation$late)

#Loading breeding success data
data <- read.csv("data/base_franklin_breed.csv",header=T,sep=";")
names(data)
head(data)

data$temp <- as.factor(data$temp)

# Missing data in the database
missmap(data, col=c("red", "green"), legend=T)

# pairwise correlations among variables
correlations <- cor(data[,6:17], use = "pairwise")
corrplot(correlations, method="circle")

# Generate subsets for each response variable (occupation, success and productivity)

data_occup <-  data[complete.cases(cbind(data$id, data$occup, data$temp,
                                                       data$dist.PPA, data$elevacion, data$slope)),]

data_exito <-  data[complete.cases(cbind(data$id, data$exito, data$temp,
                                         data$area.ppa.voronoi, data$elevacion, data$nndist, data$slope)),]


data_productividad <-  data[complete.cases(cbind(data$id, data$prod, data$temp,
                                                 data$area.ppa.voronoi, data$elevacion, data$area.ppa.circulo, data$nndist, data$slope)),]


# correlated pair of variables
cor.test(data$nndist, data$area.ppa.voronoi)
cor.test(data$dist.PPA, data$area.ppa.voronoi)
cor.test(data$elevacion, data$area.ppa.voronoi)
cor.test(data$slope, data$area.ppa.voronoi)
cor.test(data$elevacion, data$dist.PPA)
cor.test(data$orientation, data$dist.PPA)
cor.test(data$area.ppa.circulo, data$dist.PPA)


#----------------------------------------#
#### Occupation modeling ####
#----------------------------------------#

# Set of candidate models

# Distance to PPA model
glm.occ_ppa <- glmer(occup ~ dist.PPA + (1|id), 
                family=binomial(link = "logit"), data=data_occup)

# Topographic model
# check that variables are not correlated
cor.test(data$elevacion, data$slope)

glm.occ_topo <- glmer(occup ~ elevacion + slope + (1|id), 
                      family=binomial(link = "logit"), data=data_occup, na.action = "na.fail")


# Local area model
glm.occ_area <- glmer(occup ~ area.circulo.miles + (1|id), 
                       family=binomial(link = "logit"), data=data_occup, na.action = "na.fail")

# Null model 
glm.occ_null <- glmer(occup ~ 1 + (1|id), 
                      family=binomial(link = "logit"), data=data_occup, na.action = "na.fail")


# ranking of models
occ_models <-  data.frame(model = c("Distance PPA", "Topography", "Local area", "Null"),
                      AICc = c(AICc(glm.occ_ppa), AICc(glm.occ_topo), AICc(glm.occ_area), AICc(glm.occ_null)))

occ_models$deltaAICc <- occ_models$AICc-72.73063 #AICc of best model
head(occ_models)

summary_occ_models <- occ_models[order(occ_models$AICc),]
summary_occ_models

# summary of models with significant effects
summary(glm.occ_area)


#Variance explained by the best model 
#R2m (marginal): marginal R squared value associated with fixed effects
#R2c (conditional) conditional R2 value associated with fixed effects plus the random effects.
r.squaredGLMM(glm.occ_area) 


# Model diagnostics with DHARMa ####
testDispersion(glm.occ_area)
plotQQunif(glm.occ_area) # left plot in plot.DHARMa()
plotResiduals(glm.occ_area) # right plot in plot.DHARMa()


# Best model predicted values
predicho.ocupacion <-  ggpredict(glm.occ_area, terms=c("area.circulo.miles[0:12.75]"),  type="random")
predicho.ocupacion

# Plot of best model
plot_occ <- ggplot(predicho.ocupacion, aes(x=x, y=predicted)) +
  geom_line(color="blue", lwd=1.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2) +
  xlab("Local colony area (Thousands of m2)") +
  ylab("Reoccupation probability") +
  theme_classic() +
  theme(axis.text = element_text(face = "italic", size=14),
        axis.title.x = element_text(face = "bold", size=16, vjust = -2),
        axis.title.y = element_text(face = "bold", size=16, vjust = +5),
        plot.margin = margin(1,1,1.5,1.5, "cm")) 

plot_occ





#----------------------------------------#
#### Modelado éxito reproductivo ####
#----------------------------------------#

# Set of candidate models
# o. Null

glm.exito.null <- glmer(exito ~ 1 + (1|id), 
                           family=binomial(link = "logit"), data=data_exito, na.action = na.fail)

# 1. Exclusive area

glm.exito.voronoi <- glmer(exito ~ area.voronoi.miles + (1|id), 
                 family=binomial(link = "logit"), data=data_exito, na.action = na.fail)


# 2. Local area
glm.exito.ppa <- glmer(exito ~ area.circulo.miles + (1|id), 
                           family=binomial(link = "logit"), data=data_exito, na.action = na.fail)

# 3. Topographic
glm.exito.topo <- glmer(exito ~ elevacion + slope + (1|id), 
                           family=binomial(link = "logit"), data=data_exito, na.action = na.fail)

# Tabla de modelos
AICc_exito <-  data.frame(model = c("null", "Exclusive area", "Local area", "Topography"),
                        AICc = c(AICc(glm.exito.null), AICc(glm.exito.voronoi), AICc(glm.exito.ppa), AICc(glm.exito.topo)))

AIC_rank <-  AICc_exito[order(AICc_exito$AICc),]
AIC_rank
AIC_rank$deltaAICc <- AIC_rank$AICc - 55.42926
AIC_rank


#Variance explained by the best model 
#R2m (marginal): marginal R squared value associated with fixed effects
#R2c (conditional) conditional R2 value associated with fixed effects plus the random effects.
r.squaredGLMM(glm.exito.ppa) 


# Model diagnostics with DHARMa ####
testDispersion(glm.exito.ppa)
plotQQunif(glm.exito.ppa) # left plot in plot.DHARMa()
plotResiduals(glm.exito.ppa) # right plot in plot.DHARMa()


# Best model predicted values
predicho.exito <-  ggpredict(glm.exito.ppa, terms=c("area.circulo.miles[0:12.75]"),  type="random")
predicho.exito

# Plot of best model
plot_exito <- ggplot(predicho.exito, aes(x=x, y=predicted)) +
  geom_line(color="blue", lwd=1.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2) +
  xlab("Local colony area (Thousands of m2)") +
  ylab("Breeding success probability") +
  theme_classic() +
  theme(axis.text = element_text(face = "italic", size=14),
        axis.title.x = element_text(face = "bold", size=16, vjust = -2),
        axis.title.y = element_text(face = "bold", size=16, vjust = +5),
        plot.margin = margin(1,1,1.5,1.5, "cm")) 

plot_exito



#----------------------------------------#
#### Modelado productividad ####
#----------------------------------------#

# 0. Null model
glm.prod_null  <-lmer(prod ~ 1 + (1|id), 
                             data=data_productividad, na.action=na.fail)


# 1. Exclusive area model
glm.prod_voronoi <- lmer(prod ~ area.voronoi.miles + (1|id), 
                     data=data_productividad, na.action=na.fail)


# 2. Local area model
glm.prod_ppa  <-lmer(prod ~ area.circulo.miles + (1|id), 
                             data=data_productividad,na.action=na.fail)


# 3. Topographic model
glm.prod_topo  <-lmer(prod ~ elevacion + slope + (1|id), 
                         data=data_productividad,na.action=na.fail)

# Tabla de modelos
AICc_prod <-  data.frame(model = c("null", "Exclusive area", "Local area", "Topography"),
                          AICc = c(AICc(glm.prod_null), AICc(glm.prod_voronoi), AICc(glm.prod_ppa), AICc(glm.prod_topo)))

AIC_rank <-  AICc_prod[order(AICc_prod$AICc),]
AIC_rank
AIC_rank$deltaAICc <- AIC_rank$AICc - 136.96
AIC_rank



#Variance explained by the best model 
#R2m (marginal): marginal R squared value associated with fixed effects
#R2c (conditional) conditional R2 value associated with fixed effects plus the random effects.
r.squaredGLMM(glm.prod_null) 

# prediction by season
predicho.prod <-  ggpredict(glm.prod_null, terms=c("temp"),  type="fixed")
predicho.prod

# Local area model
predicho.prod <-  ggpredict(glm.prod_ppa, terms=c("area.circulo.miles[0:12.75]"),  type="random")
predicho.prod

# Plot of best model
plot_prod <- ggplot(predicho.prod, aes(x=x, y=predicted)) +
  geom_line(color="blue", lwd=1.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2) +
  xlab("Local colony area (Thousands of m2)") +
  ylab("Productivity (chicks/nest)") +
  theme_classic() +
  theme(axis.text = element_text(face = "italic", size=14),
        axis.title.x = element_text(face = "bold", size=16, vjust = -2),
        axis.title.y = element_text(face = "bold", size=16, vjust = +5),
        plot.margin = margin(1,1,1.5,1.5, "cm")) 

plot_prod


# Plot both models

figure <-  ggarrange(plot_occ, plot_exito, 
                         labels = c("A", "B"),
                         ncol = 1, nrow = 2)

figure






