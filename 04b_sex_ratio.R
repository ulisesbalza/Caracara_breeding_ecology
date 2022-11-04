rm(list=ls())


#### packages needed ####

library(gplots)
library(nlme)
library(MuMIn)
library(Amelia)
library(mlbench)
library(corrplot)
library(lme4)
library(ggeffects)
library(ggplot2)
library(MASS)

# loading database

data_sex <- read.csv("data/base_sex_ratio.csv",header=T,sep=";")
head(data_sex)

# exploring correlations
correlations <- cor(data_sex[,4:5], use = "pairwise")
corrplot(correlations, method="circle")

cor.test(data_sex$brood.size, data_sex$brood.reduction)


# set of models
glmm.brood <- glmer(male ~ brood.size + (1|temp) + (1|nest), #the only one with significant marginal effects
                  family=binomial(link = "logit"), data=data_sex)
summary(glmm.brood)


glmm.reduction <- glmer(male ~ brood.reduction + (1|temp) + (1|nest), 
                                family=binomial(link = "logit"), data=data_sex)
summary(glmm.reduction)


glmm.sex.null <- glmer(male ~ 1 + (1|temp) + (1|nest),
                    family=binomial(link = "logit"), data=data_sex)

# Tabla de modelos
AICc_sex <-  data.frame(model = c("brood size", "brood reduction", "null"),
                         AICc = c(AICc(glmm.brood), AICc(glmm.reduction), AICc(glmm.sex.null)))

AIC_rank <-  AICc_sex[order(AICc_sex$AICc),]
AIC_rank
AIC_rank$deltaAICc <- AIC_rank$AICc - 91.46849
AIC_rank


#predicted values
predicho.sex_brood.glm <-  ggpredict(glmm.brood, terms=c("brood.size [1,2,3]"),  type="fixed")
predicho.sex_brood.glm


#Varianza explicada por los modelos
#R2m: marginal R squared value associated with fixed effects
#R2c conditional R2 value associated with fixed effects plus the random effects.
#R2c, delta

r.squaredGLMM(glmm.brood)


#Rbase
boxplot(predicted ~  x, data=predicho.sex_brood.glm, ylim=c(0,1), pch=19, cex=1.2, xlab="Brood size (chicks)", ylab="Estimated proportion of males")
points(conf.low ~  x, data=predicho.sex_brood.glm, pch=19, cex=1)
points(conf.high ~  x, data=predicho.sex_brood.glm, pch=19, cex=1)

#ggplot
ggplot(predicho.sex_brood.glm, aes(x, predicted)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), alpha = 1,position = position_dodge())+
  labs(
    x="Brood size (chicks/sucessful nest)",
    y="Predicted proportion of males") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold"),
        panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        axis.line = element_line(color = "black"))




