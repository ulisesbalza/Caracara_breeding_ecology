rm(list=ls())

# Required libraries
library(ggplot2)
library(RColorBrewer)
library(mgcv)
library(DHARMa)
library(ggpubr)
library(dplyr)
library(nlme)
library(lme4)
library(MuMIn)
library(MASS)
library(boot)
library(ggeffects)
library(insight)

library(gplots)
library(Amelia)
library(mlbench)
library(corrplot)
library(mgcViz)
library(gridExtra)


###
#### REOCCUPATION OF SITES #####
###

occup <- read.csv("data/base_occup.csv",header=T,sep=";")
occup$temp <- as.factor(occup$temp)

# set of models

glm.occ.null <- glmer(occup ~ 1 + (1|id) + (1|temp), 
                      family=binomial(link = "logit"), data=occup, na.action = na.fail)

summary(glm.occ.null) #temp does not have explanation power, so I remove it 

glm.occ.null <- glmer(occup ~ 1 + (1|id), 
                      family=binomial(link = "logit"), data=occup, na.action = na.fail)

glm.occ.distPPA <- glmer(occup ~ dist.PPA + (1|id), 
                      family=binomial(link = "logit"), data=occup)

glm.occ.topo <- glmer(occup ~  cliff + slope + (1|id), 
                         family=binomial(link = "logit"), data=occup)


# ranking of models
occ_models <-  data.frame(model = c("Null", "Distance PPA", "Topographic"),
                           AICc = c(AICc(glm.occ.null), AICc(glm.occ.distPPA), AICc(glm.occ.topo)))

occ_models$deltaAICc <- (occ_models$AICc)-min(occ_models$AICc) 
summary_occ_models <- occ_models[order(occ_models$AICc),]
summary_occ_models


# Global estimates
summary(glm.occ.null)
inv.logit(1.1833)
confint.merMod(glm.occ.null)
inv.logit(0.3717642)
inv.logit(2.420877)
#Variance explained by null model
r.squaredGLMM(glm.occ.null) 

# Diagnosis of the best model
testDispersion(glm.occ.distPPA)
plotQQunif(glm.occ.distPPA) 

# Variance explained by best model 
r.squaredGLMM(glm.occ.distPPA) 
summary(glm.occ.distPPA)

# predicted values
occup_predicted <-  ggpredict(glm.occ.distPPA, terms=c("dist.PPA"),  type="fixed")
occup_predicted

# Plot of best model
plot_occup_ppa <- ggplot(occup_predicted, aes(x=x, y=predicted)) +
  geom_point(data = occup, aes(x = dist.PPA, y = occup)) +
  geom_line(color="black", lwd=1.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2) +
  xlab("Distance to penguin colony (m)") +
  ylab("Reoccupation 
       probability") +
  ylim(0,1) +
  theme_classic() +
  theme(axis.text = element_text(face = "italic", size=14),
        axis.title.x = element_text(face = "bold", size=14, vjust = -2),
        axis.title.y = element_text(face = "bold", size=14, vjust = +5),
        plot.margin = margin(1,1,1.5,1.5, "cm")) 

plot_occup_ppa


###
#### CLUTCH SIZE #####
###

clutch <- read.csv("data/base_clutch.csv",header=T,sep=";")
clutch$temp <- as.factor(clutch$temp)

# Set of models
glm.clutch.null  <-lmer(clutch ~ 1 + (1|id) + (1|temp), 
                      data=clutch, na.action=na.fail)

# Global estimates
summary(glm.clutch.null)
confint(glm.clutch.null)

#Variance explained 
r.squaredGLMM(glm.clutch.null) 
clutch_total_var <- 0.15796 + 0.03009 + 0.19405
clutch_explained_by_temp <- 0.03009/clutch_total_var
clutch_explained_by_id <- 0.19405/clutch_total_var

# predicted values for each nest
clutch_predicted <-  ggpredict(glm.clutch.null, terms=c("id"),  type="random")
clutch_predicted

# Plot by id
plot_clutch_id <- ggplot(clutch_predicted, aes(x = reorder(x, predicted), y=predicted)) +
  geom_point(size=5) +
 # geom_pointrange(aes(ymin = conf.low, ymax = conf.high), size=1,color="id") +
  xlab("Nest id") +
  ylab("Clutch size") +
  ylim(1,3) +
  geom_hline(yintercept = 1.97, linetype='dotted', lwd=2) +
  theme_classic() +
  theme(axis.text = element_text(face = "italic", size=14),
        axis.text.x=element_blank(),
        axis.title.x = element_text(face = "bold", size=14, vjust = -2),
        axis.title.y = element_text(face = "bold", size=14, vjust = +5),
        legend.position="none",
        plot.margin = margin(1,1,1.5,1.5, "cm")) 

plot_clutch_id

# plot by year
glm.clutch.temp  <-lmer(clutch ~ temp + (1|id), 
                        data=clutch, na.action=na.fail)

clutch_predicted_temp <-  ggpredict(glm.clutch.temp, terms=c("temp"),  type="fixed")
clutch_predicted_temp

plot_clutch <- ggplot(clutch_predicted_temp, aes(x = x, y= predicted)) +
  geom_point(size=3) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), size=1,color="black") +
  xlab("Year") +
  ylab("Clutch size") +
  ylim(1,3) +
  geom_hline(yintercept = 1.96, linetype='dotted', lwd=2) +
  theme_classic() +
  theme(axis.text = element_text(face = "italic", size=12),
        axis.title.x = element_text(face = "bold", size=14, vjust = -2),
        axis.title.y = element_text(face = "bold", size=14, vjust = +5),
        plot.margin = margin(1,1,1.5,1.5, "cm")) 

plot_clutch

###
#### HATCHING SUCCESS ####
###

hatch <- read.csv("data/base_hatching.csv",header=T,sep=";")
hatch$temp <- as.factor(hatch$temp)

# set of models

glm.hatch.null <- glmer(hatch ~ 1 + (1|id) + (1|temp), 
                      family=binomial(link = "logit"), data=hatch, na.action = na.fail)

summary(glm.hatch.null) #temp does not have explanation power, so I remove it 

glm.hatch.null <- glmer(hatch ~ 1 + (1|id), 
                        family=binomial(link = "logit"), data=hatch, na.action = na.fail)

glm.hatch.local <- glmer(hatch ~ area.circulo.miles + (1|id), 
                          family=binomial(link = "logit"), data=hatch, na.action = na.fail)

glm.hatch.exclusive <- glmer(hatch ~ area.voronoi.miles + (1|id), 
                       family=binomial(link = "logit"), data=hatch, na.action = na.fail)

glm.hatch.topo <- glmer(hatch ~ cliff + slope + (1|id), 
                           family=binomial(link = "logit"), data=hatch, na.action = na.fail)


# ranking of models
hatch_models <-  data.frame(model = c("Null", "Local area", "Exclusive area", "Topographic"),
                          AICc = c(AICc(glm.hatch.null), AICc(glm.hatch.local), AICc(glm.hatch.exclusive), AICc(glm.hatch.topo)))

hatch_models$deltaAICc <- (hatch_models$AICc)-min(hatch_models$AICc) 
summary_hatc_models <- hatch_models[order(hatch_models$AICc),]
summary_hatc_models


# Global estimates
summary(glm.hatch.null)
inv.logit(1.2490)
confint(glm.hatch.null)
inv.logit(0.7048606)
inv.logit(1.924476)
#Variance explained by null model
r.squaredGLMM(glm.hatch.null) 

# Diagnosis of the best model
testDispersion(glm.hatch.local)
plotQQunif(glm.hatch.local) 

#Variance explained by the best model 
r.squaredGLMM(glm.hatch.local) 

# predicted values
hatch_ppa_predicted <-  ggpredict(glm.hatch.local, terms=c("area.circulo.miles"),  type="fixed")
hatch_ppa_predicted

# Plot of best model
plot_hatch_ppa <- ggplot(hatch_ppa_predicted, aes(x=x, y=predicted)) +
  geom_point(data = hatch, aes(x = area.circulo.miles, y = hatch)) +
  geom_line(color="black", lwd=1.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2) +
  xlab("Local colony area (Thousands of m2)") +
  ylab("Hatching probability") +
  ylim(0,1) +
  theme_classic() +
  theme(axis.text = element_text(face = "italic", size=14),
        axis.title.x = element_text(face = "bold", size=14, vjust = -2),
        axis.title.y = element_text(face = "bold", size=14, vjust = +5),
        plot.margin = margin(1,1,1.5,1.5, "cm")) 

plot_hatch_ppa


###
#### FLEDGLING PROBABILITY ####
###

fled <- read.csv("data/base_fledgling.csv",header=T,sep=";")
fled$temp <- as.factor(fled$temp)

# set of models

glm.fled.null <- glmer(fled ~ 1 + (1|id) + (1|temp), 
                        family=binomial(link = "logit"), data=fled, na.action = na.fail)

summary(glm.fled.null) #temp does not have explanation power, so I remove it 

glm.fled.null <- glmer(fled ~ 1 + (1|id), 
                       family=binomial(link = "logit"), data=fled, na.action = na.fail)

glm.fled.local <- glmer(fled ~ area.circulo.miles + (1|id), 
                         family=binomial(link = "logit"), data=fled, na.action = na.fail)

glm.fled.exclusive <- glmer(fled ~ area.voronoi.miles + (1|id), 
                             family=binomial(link = "logit"), data=fled, na.action = na.fail)

glm.fled.topo <- glmer(fled ~ cliff + slope + (1|id), 
                        family=binomial(link = "logit"), data=fled, na.action = na.fail)


# ranking of models
fled_models <-  data.frame(model = c("Null", "Local area", "Exclusive area", "Topographic"),
                            AICc = c(AICc(glm.fled.null), AICc(glm.fled.local), AICc(glm.fled.exclusive), 
                                     AICc(glm.fled.topo)))

fled_models$deltaAICc <- (fled_models$AICc)-min(fled_models$AICc) 
summary_fled_models <- fled_models[order(fled_models$AICc),]
summary_fled_models

# Global estimates
summary(glm.fled.null)
inv.logit(3.098)
confint.merMod(glm.fled.null, method = "boot")
inv.logit(1.71515)
inv.logit(9.431313)

#Variance explained by null model
r.squaredGLMM(glm.fled.null) 

# Diagnosis of the null model
testDispersion(glm.fled.null)
plotQQunif(glm.fled.null) 

# predicted values for each territory
fled_predicted <-  ggpredict(glm.fled.null, terms=c("id"),  type="random")
fled_predicted

# Plot by id
plot_fled_id <- ggplot(fled_predicted, aes(x = reorder(x, predicted), y=predicted, col=fled_predicted$x)) +
  geom_point(size=5) +
  # geom_pointrange(aes(ymin = conf.low, ymax = conf.high), size=1,color="id") +
  xlab("Nest id") +
  ylab("Fledling probability") +
  ylim(0,1) +
  geom_hline(yintercept = 0.96, linetype='dotted', lwd=2) +
  theme_classic() +
  theme(axis.text = element_text(face = "italic", size=14),
        axis.text.x=element_blank(),
        axis.title.x = element_text(face = "bold", size=14, vjust = -2),
        axis.title.y = element_text(face = "bold", size=14, vjust = +5),
        legend.position="none",
        plot.margin = margin(1,1,1.5,1.5, "cm")) 

plot_fled_id

###
#### PRODUCTIVITY ####
###

prod <- read.csv("data/base_productivity.csv",header=T,sep=";")
prod$temp <- as.factor(prod$temp)

# models
glm.prod.null  <-lmer(prod ~ 1 + (1|id) + (1|temp), 
                      data=prod, na.action=na.fail)
glm.prod_exclusive <- lmer(prod ~ area.voronoi.miles + (1|id) + (1|temp), 
                      data=prod)
glm.prod_local <- lmer(prod ~ area.ppa.miles + (1|id) + (1|temp), 
                         data=prod)
glm.prod_topo <- lmer(prod ~ cliff + slope + (1|id) + (1|temp), 
                       data=prod)

# ranking of models
prod_models <-  data.frame(model = c("Null", "Exclusive area", "Local area", "Topographic"),
                          AICc = c(AICc(glm.prod.null), AICc(glm.prod_exclusive), AICc(glm.prod_local), AICc(glm.prod_topo)))

prod_models$deltaAICc <- (prod_models$AICc)-min(prod_models$AICc) 
summary_prod_models <- prod_models[order(prod_models$AICc),]
summary_prod_models

# Global estimates
summary(glm.prod.null)
confint(glm.prod.null)

# % variance explained by year and id
prod_total_var <- 0.24795 + 0.02611 + 0.53018
prod_id_var <- 0.24795/prod_total_var
prod_temp_var <- 0.02611/prod_total_var  

# predicted values by id
prod_predicted <-  ggpredict(glm.prod_exclusive, terms=c("area.voronoi.miles"),  type="fixed")
prod_predicted

# Plot of exclusive area
plot_prod <- ggplot(prod_predicted, aes(x=x, y=predicted)) +
  geom_point(data = prod, aes(x = area.voronoi.miles, y = prod)) +
  geom_line(color="black", lwd=1.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2) +
  xlab("Exclusive colony area (Thousands of m2)") +
  ylab("Productivity
       (chicks/nest)") +
  ylim(0,3) +
  theme_classic() +
  theme(axis.text = element_text(face = "italic", size=14),
        axis.title.x = element_text(face = "bold", size=13, vjust = -2),
        axis.title.y = element_text(face = "bold", size=14, vjust = +5),
        plot.margin = margin(1,1,1.5,1.5, "cm")) 

plot_prod

###
#### CHICK SEX RATIO ####
###

sex <- read.csv("data/base_sex_ratio.csv",header=T,sep=";")
sex$temp <- as.factor(sex$temp)
sex$brood.size <- as.ordered(sex$brood.size)


# set of models

glm.sex.null <- glmer(male ~ (1|nest) + (1|temp), 
                      family=binomial(link = "logit"), data=sex, na.action = na.fail)

summary(glm.sex.null) #temp does not have explanation power, so I remove it 

glm.sex.null <- glmer(male ~ (1|nest), 
                      family=binomial(link = "logit"), data=sex, na.action = na.fail)

glm.sex.brood <- glmer(male ~ brood.size + (1|nest), 
                      family=binomial(link = "logit"), data=sex, na.action = na.fail)

glm.sex.reduction <- glmer(male ~ brood.reduction + (1|nest), 
                      family=binomial(link = "logit"), data=sex, na.action = na.fail)


# ranking of models
sex_models <-  data.frame(model = c("Null", "Brood size", "Brood reduction"),
                           AICc = c(AICc(glm.sex.null), AICc(glm.sex.brood), AICc(glm.sex.reduction)))

sex_models$deltaAICc <- (sex_models$AICc)-min(sex_models$AICc) 
summary_sex_models <- sex_models[order(sex_models$AICc),]
summary_sex_models

# Global estimates
summary(glm.sex.null)
inv.logit(-0.3719)
confint.merMod(glm.sex.null)
inv.logit(-1.2963)
inv.logit(0.2907489)

#Variance explained by null model
r.squaredGLMM(glm.sex.null) 

#Variance explained by best model
r.squaredGLMM(glm.sex.brood) 


# Diagnosis of the null model
testDispersion(glm.sex.null)
plotQQunif(glm.sex.null) 

#predicted values
sex_predicted <-  ggpredict(glm.sex.brood, terms=c("brood.size"),  type="fixed")
sex_predicted

# Plot of best model
plot_sex <- ggplot(sex_predicted, aes(x=x, y=predicted)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), size=1,color="black") +
  xlab("Brood size") +
  ylab("Proportion of males") +
  ylim(0,1) +
  theme_classic() +
  theme(axis.text = element_text(face = "italic", size=14),
        axis.title.x = element_text(face = "bold", size=14, vjust = -2),
        axis.title.y = element_text(face = "bold", size=14, vjust = +5),
        plot.margin = margin(1,1,1.5,1.5, "cm")) 

plot_sex


# all plots together

ggarrange(plot_clutch_id, plot_clutch, plot_occup_ppa, plot_hatch_ppa, plot_prod, plot_sex,
          labels = c("A", "B", "C", "D", "E", "F"), align = "v",
          ncol = 2, nrow = 3)

