#Code & notes by Katie Flowers & Beth Babcock
#Last updated 04 November 2024

##### PREP #####
library(MASS) #negative binomial glm
library(glmmTMB) #tweedie
library(lmtest) #likelihood ratio test
library (MuMIn) #model averaging
library(DHARMa) #residual diagnostics
library(performance) #model check, multicollinearity check
library(patchwork) #combine ggplots
library(ggplot2) #data visualization
library(sjPlot) #forest plots

#citation("DHARMa") #check how to cite packages in pubs

rays<-read.csv("Mean ray covariates_23 Sept 24.csv") #upload data
head(rays) #check header

#make sure numerical variables are numeric 
rays$vis<-as.numeric(rays$vis)
rays$depth<-as.numeric(rays$depth)
rays$temp<-as.numeric(rays$temp)
rays$relief<-as.numeric(rays$relief)
rays$gravity<-as.numeric(rays$gravity)
rays$n_spp<-as.numeric(rays$n_rayspp)

#make categorical variables factors
rays$bait<-factor(rays$bait)
rays$geomorphology<-factor(rays$geomorphology)
rays$reef_protection_status<-factor(rays$reef_protection_status)
rays$ray_fishing<-factor(rays$ray_fishing)
rays$rays_landed<-factor(rays$rays_landed)
rays$bottom_gear<-factor(rays$bottom_gear)
rays$country<-factor(rays$country)

#check levels of factors
levels(rays$bait)
table(rays$bait) #this probably won't work as covariate - note barracuda, gar, herring, tuna locations as different 

#check distributions
hist(rays$ray_maxn) 
hist(rays$shark_maxn)
hist(rays$vis)
hist(rays$gravity)
hist(rays$depth)
hist(rays$temp)
hist(rays$relief)
hist(rays$n_rayspp)

##### RAY MAXN MODELS - REEF LEVEL #####

#Chose mean reef level due to several "nas" in set level numerical covariate data
#Problems when bait included as predicted from above

#Negative binomial GLM
test1<-glm.nb(ray_maxn ~ vis + geomorphology + shark_maxn + temp + depth + 
              relief + gravity + reef_protection_status + rays_landed + 
              ray_fishing + bottom_gear,
              data = rays)

#DHARMa residuals
par(mfrow=c(2,2))
simtest1<- simulateResiduals(fittedModel = test1, n = 250) 
plot(simtest1) #problematic

#Poisson GLM
test2<-glm(ray_maxn ~ vis + geomorphology + shark_maxn + temp + depth + 
                relief + gravity + reef_protection_status + rays_landed + 
                ray_fishing + bottom_gear,
                family = "poisson",
                data = rays)

#DHARMa residuals
simtest2<- simulateResiduals(fittedModel = test2, n = 250) 
plot(simtest2) #problematic

#Tweedie GLM
test3<-glmmTMB(ray_maxn ~ vis + geomorphology + shark_maxn + temp + depth + 
                relief + gravity + reef_protection_status + rays_landed + 
                ray_fishing + bottom_gear,
                family = "tweedie",
                data = rays)

#DHARMa residuals
simtest3<- simulateResiduals(fittedModel = test3, n = 250) 
plot(simtest3) #no significant issues

#Tweedie model from above tests with random effect
mod1<-glmmTMB(ray_maxn ~ vis + geomorphology + shark_maxn + temp + depth + 
              relief + gravity + reef_protection_status + rays_landed + 
              ray_fishing + bottom_gear + (1|country),
              family="tweedie", 
              data = rays)

#DHARMa residuals
simtest4<- simulateResiduals(fittedModel = mod1, n = 250) 
plot(simtest4) #quantile deviations

#check model, multicollinearity 
check_model(mod1)
c<-check_collinearity(mod1)
plot(c) #take out ray_fishing - collinear w/protection status, makes sense

#new model without ray_fishing and without random effect
mod2<-glmmTMB(ray_maxn ~ vis + geomorphology + shark_maxn + temp + depth + 
                relief + gravity + reef_protection_status + rays_landed + 
                bottom_gear,
                family = "tweedie",
                data = rays)

#check model, multicollinearity 
check_model(mod2)
c2<-check_collinearity(mod2)
plot(c2) #looks good now ; devs shared in 2023 error is a misleading msg (https://github.com/easystats/performance/issues/545)

#check residuals with DHARMa
par(mfrow=c(2,2))
simmod2 <- simulateResiduals(fittedModel = mod2, n = 250)
plot(simmod2) #okay
summary(mod2)

#new model without ray_fishing and with random effect
mod3<-glmmTMB(ray_maxn ~ vis + geomorphology + shark_maxn + temp + depth + 
              relief + gravity + reef_protection_status + rays_landed + 
              bottom_gear + (1|country),
              family="tweedie", 
              data = rays) 
#Brooks et al. vignette okay to ignore the Na/NaN function eval warnings

#check model, multicollinearity 
check_model(mod3)
c3<-check_collinearity(mod3)
plot(c3) #looks good

#check residuals with DHARMa
par(mfrow=c(2,2))
simmod3 <- simulateResiduals(fittedModel = mod3, n = 250)
plot(simmod3) #okay
summary(mod3)

#Is model with random effect better?
AIC(mod2, mod3) #mod3
BIC(mod2, mod3) #mod2
#AICs and BICs very close for both models, should do further checks

#Likelihood ratio test for goodness of fit
lrtest(mod2, mod3) #fail to reject null at 0.05 level, therefore use simple model w/o random effect
#but close to 0.05, try more

#Leave one out cross validation
errors <- numeric(nrow(rays)) #for prediction error

#Loop through each observation, leaving one out at a time
##mod2
for (i in 1:nrow(rays)) {
  train_data <- rays[-i, ] 
  test_data <- rays[i, , drop = FALSE] #exclude the ith observation
  mod2_loo <- glmmTMB(ray_maxn ~ vis + geomorphology + shark_maxn + temp + depth + 
                        relief + gravity + reef_protection_status + rays_landed + 
                        bottom_gear,
                        family = "tweedie",
                        data = train_data) #fit model training data
  pred <- predict(mod2_loo, newdata = test_data, type = "response") #predict left out obs
  errors[i] <- (test_data$ray_maxn - pred)^2  #store squared prediction error
}

#Calculate mean squared error (MSE)
mse_loomod2 <- mean(errors)
print(mse_loomod2)

##mod3
for (i in 1:nrow(rays)) {
  train_data <- rays[-i, ] 
  test_data <- rays[i, , drop = FALSE] #exclude the ith observation
  mod3_loo <- glmmTMB(ray_maxn ~ vis + geomorphology + shark_maxn + temp + depth + 
                        relief + gravity + reef_protection_status + rays_landed + 
                        bottom_gear + (1|country),
                        family = "tweedie",
                        data = train_data) #fit model training data
  pred <- predict(mod3_loo, newdata = test_data, type = "response") #predict left out obs
  errors[i] <- (test_data$ray_maxn - pred)^2  #store squared prediction error
}

#Calculate mean squared error (MSE)
mse_loomod3 <- mean(errors)
print(mse_loomod3)

#Calculate MSE for other models
#get predictions
pred_mod2 <- predict(mod2, newdata = rays, type = "response")
pred_mod3 <- predict(mod3, newdata = rays, type = "response")
#get squared errors
squared_errors_mod2 <- (rays$ray_maxn - pred_mod2)^2 
squared_errors_mod3 <- (rays$ray_maxn - pred_mod3)^2
#MSE
mse_mod2 <- mean(squared_errors_mod2) 
mse_mod3 <- mean(squared_errors_mod3) 

#table w/results
mse_table <- data.frame (Model = c("mod2", "mod3", "mod2_loo", "mod3_loo"),
                         MSE = c(mse_mod2,mse_mod3, mse_loomod2, mse_loomod3))
mse_table #mod3 better based on MSE, mod2 better based on leave one out CV but MSEs very close

#model performance similar across multiple criteria, but I think that there is likely 
#variation between countries, selecting mod3 w/random effect

#model averaging
options(na.action = "na.fail")

#model selection table
ms<-dredge(mod3, rank = "AIC")

#supported models (delta AIC < or = 2)
sub<-subset(ms, delta <= 2, recalc.weights = T) 

#model avg coefficients from supported models only 
mod3_avg<-model.avg(sub)
summary(mod3_avg) #not a lot, makes sense b/c pooled maxn similar across countries

#calculate weights
checkw<-model.avg(sub, rank = "AIC") #getting weights only for supported models
Weights(checkw) #weights should sum to 1 (cross-check w/mod3 summary results)

#best model based on delta AIC & weights = relief
bestmod<-glmmTMB(ray_maxn ~ relief + (1|country),
                 family = "tweedie",
                 data = rays)
summary(bestmod)
car::Anova(bestmod) #wald X2

#DATA VIS - EFFECTS
require(ggeffects)
require(ggplot2)
require(patchwork)

#Predicted effects for each covariate with 95% CIs
relief.effect <- ggpredict(bestmod, terms = "relief")

#Plot w/95% CIs
relief.plot <- ggplot(relief.effect, aes (x=x, y=predicted)) +
               ggtitle("") +
               geom_line(color = "black", size = 1.1) +  #hmm color not changing
               geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
                               fill = "#9A8895", alpha = 0.5) +  
               xlab("\nMean Reef Relief") +
               ylab("Predicted Mean Ray MaxN\n") + 
               scale_x_continuous(breaks = seq(min(relief.effect$x), 
                                               max(relief.effect$x), 
                                               by = 0.5)) +
               theme_classic() +                             
               theme(text = element_text(size = 20, color = "black"),
                     axis.text.x = element_text(color = "black"),
                     axis.text.y = element_text(color = "black"))
relief.plot

#check predicted
head(relief.effect)
print(relief.effect, n = Inf)

#rows with min max relief
min_row <- which.min(relief.effect$x) #maxn 0.27 at row 1
max_row <- which.max(relief.effect$x) #maxn 0.14 at row 17

#Calculate relative change (ratio b/t high & low relief) between two points 
#ratio <1 = decrease
#ratio >1 = increase
relative_change <- relief.effect$predicted[max_row] / relief.effect$predicted[min_row]
relative_change #0.53 (maxn at high relief = 53% of maxn at low relief)

#calculating percentage decrease in ray MaxN from low to high relief
(1-relative_change) * 100 #47.3% decrease in maxn from low to high relief

#Variance explained is more complicated. This is basically R-squared, which
#is not easy to calculate from GLMM models in which the residual variance is not constant. 
#There are some proposed methods in the literature. The library MuMin library calcuates
#some of them, that have good help files
?r.squaredGLMM
?r.squaredLR
#Unfortunately these don't seem to work for glmmtmb tweedie models.
#But the performance model has one that works and also does Nakagawa's method
require(performance)
r2(bestmod) #cond = 0.352, marg = 0.068
# Marginal is the variance explained by the fixed effects (relief), 
# and conditional is variance explained by both relief and the random effect together


#####HYPANUS AMERICANUS MAXN MODELS - REEF LEVEL#####

#remove Brazil from data since H. americanus is not there
raysnoBR<-rays[-c(26:32),]

#Tweedie
Ham1<-glmmTMB(Ham_maxn ~ vis + geomorphology + shark_maxn + temp + depth + 
                relief + gravity + ray_fishing + reef_protection_status + 
                rays_landed + bottom_gear,
                family="tweedie", 
                data = raysnoBR)

#check multicollinearity 
Hamc1<-check_collinearity(Ham1)
plot(Hamc1) #remove ray fishing

#model without ray fishing
Ham2<-glmmTMB(Ham_maxn ~ vis + geomorphology + shark_maxn + temp + depth + 
                relief + gravity + reef_protection_status + rays_landed + 
                bottom_gear + (1|country),
                family="tweedie", 
                data = raysnoBR)

#check multicollinearity 
Hamc2<-check_collinearity(Ham2)
plot(Hamc2) #all good

#check residuals with DHARMa
par(mfrow=c(2,2))
simHammod <- simulateResiduals(fittedModel = Ham2, n = 250)
plot(simHammod) #all good

#model averaging
options(na.action = "na.fail")

#model selection table
msHam<-dredge(Ham2, rank = "AIC")

#supported models (delta AIC < or = 2)
subHam<-subset(msHam, delta <= 2, recalc.weights = T) 

#model avg coefficients from supported models only 
HamAvg<-model.avg(subHam)
summary(HamAvg) #islands relative to continental(+), bottom gear, relief, shark maxn (-)
levels(raysnoBR$geomorphology)

#calculate weights
checkw2<-model.avg(subHam, rank = "AIC") #getting weights only for supported models
Weights(checkw2) #weights should sum to 1

#bestmodel
HamBest<- glmmTMB(Ham_maxn ~ geomorphology + shark_maxn + temp + relief +  
                    bottom_gear + (1|country),
                  family="tweedie", 
                  data = raysnoBR)
car::Anova(HamBest) #wald X2

#DATA VIS - EFFECTS
require(ggeffects)
require(ggplot2)
require(patchwork)

#Predicted effects for each covariate with 95% CIs
geo.effect <- ggpredict(HamBest, terms = "geomorphology")
shark.effect <- ggpredict(HamBest, terms = "shark_maxn")
temp.effect <- ggpredict(HamBest, terms = "temp")
relief.effect.Ham <- ggpredict(HamBest, terms = "relief")
gear.effect <- ggpredict(HamBest, terms = "bottom_gear")

#Plots w/95% CIs
#Geomorphology
geo.plot <- ggplot(geo.effect, aes(x = x, y = predicted)) +
                   geom_point(color = "#638E4A", size = 3) +  
                   geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                                     color = "#638E4A", width = 0.1) +  
                   ggtitle("") +
                   xlab("\nGeomorphology") +
                   ylab("Mean Southern Stingray MaxN\n") + 
                   theme_classic() +                             
                   theme(text = element_text(size = 14, color = "black"),
                         axis.text.x = element_text(color = "black"),
                         axis.text.y = element_text(color = "black")) 
geo.plot

#Shark MaxN
shark.plot <- ggplot(shark.effect, aes (x=x, y=predicted)) +
              ggtitle("") +
              geom_line(color = "black", size = 1.1) +  
              geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
                              fill = "#123B6D", alpha = 0.5) +  
              xlab("\nMean Shark MaxN") +
              ylab("") + 
              scale_x_continuous(breaks = seq(min(shark.effect$x), 
                                              max(shark.effect$x), 
                                              by = 0.5)) +
              theme_classic() +                             
              theme(text = element_text(size = 14, color = "black"),
                    axis.text.x = element_text(color = "black"),
                    axis.text.y = element_text(color = "black"))
shark.plot

#SST
temp.plot <- ggplot(temp.effect, aes (x=x, y=predicted)) +
             ggtitle("") +
             geom_line(color = "black", size = 1.1) +  
             geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
                             fill = "#E4C167", alpha = 0.5) +  
             xlab("\nMean Temperature (°C)") +
             ylab("") + 
             scale_x_continuous(breaks = seq(min(temp.effect$x), 
                                  max(temp.effect$x), 
                                  by = 1)) +
             theme_classic() +                             
             theme(text = element_text(size = 14, color = "black"),
                   axis.text.x = element_text(color = "black"),
                   axis.text.y = element_text(color = "black"))
temp.plot

#Reef Relief
reliefHam.plot <- ggplot(relief.effect.Ham, aes (x=x, y=predicted)) +
                  ggtitle("") +
                  geom_line(color = "black", size = 1.1) +  
                  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
                  fill = "#9A8895", alpha = 0.5) +  
                  xlab("\nMean Reef Relief") +
                  ylab("Mean Southern Stingray MaxN\n") + 
                  scale_x_continuous(breaks = seq(min(relief.effect.Ham$x), 
                                                  max(relief.effect.Ham$x), 
                                                  by = 0.5)) +
                  theme_classic() +                             
                  theme(text = element_text(size = 14, color = "black"),
                  axis.text.x = element_text(color = "black"),
                  axis.text.y = element_text(color = "black"))
reliefHam.plot

#Bottom gear
gear.plot <- ggplot(gear.effect, aes(x = x, y = predicted)) +
             geom_point(color = "#6D0712", size = 3) +  
             geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                               color = "#6D0712", width = 0.1) +  
             ggtitle("") +
             xlab("\nBottom Gear") +
             ylab("") + 
             theme_classic() +                             
             theme(text = element_text(size = 14, color = "black"),
                   axis.text.x = element_text(color = "black"),
                   axis.text.y = element_text(color = "black")) 
gear.plot

#Combine plots & add main title
Ham.plot <- (reliefHam.plot | shark.plot | temp.plot) /
             (geo.plot | gear.plot) + 
            plot_annotation(theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))
Ham.plot

#DATA VIS - standardized estimates (Ludecke 2022)
theme_set(theme_classic())

#plot without axis labels first to see order, then add in
#Labels start at bottom
#ref levels for cat variables = continental (geo), no bottom gear
Stan.Est.Ham <- plot_model(HamBest, 
                           vline.color = "black", 
                           sort.est = TRUE, 
                           transform = NULL,
                           title = "",
                           axis.labels = c("Shark MaxN",
                                           "Bottom Gear",
                                           "Reef Relief",
                                           "Temperature (°C)",
                                           "High Island",
                                           "Low Island"))

newtheme = theme(axis.text=element_text(size = 16, color = "black"),
                 axis.text.y=element_text(size = 16, color = "black"),
                 axis.title.x = element_text(size = 20, color = "black"))

Stan.Est.Ham + newtheme + labs(y = "\nStandardized Estimates\n") 

#from global model
Stan.Est <- plot_model(Ham2, 
                       vline.color = "black", 
                       sort.est = TRUE, 
                       transform = NULL,
                       title = "")
Stan.Est

#Covariate effects

#Geomorphology
#check predicted
print(geo.effect, n = Inf)

#Calculate relative change between two points 
relative_change_geo1 <- 0.16 / 0.15 #from high to low islands
relative_change_geo1 #1.067 (greater than one so increase)

relative_change_geo2 <- 0.16 / 0.06 #from continental to low
relative_change_geo2 #2.67 (greater than one so increase)

relative_change_geo3 <- 0.15 / 0.06 #from continental to high
relative_change_geo3 #2.5 (greater than one so increase)

#calculating percentage increase in ray MaxN 
(relative_change_geo1-1) * 100 #6.7% from high to low islands
(relative_change_geo2-1) * 100 #166.7% from continental to low islands
(relative_change_geo3-1) * 100 #150% from continental to high islands

#Shark MaxN
#check predicted
print(shark.effect, n = Inf)

#Calculate relative change between two points 
relative_change_shark <- 0.02 / 0.08 #from low to high shark MaxN
relative_change_shark #0.25 (less than one so decrease)

#calculating percentage decrease in ray MaxN 
(1-relative_change_shark) * 100 #75% from low to high shark MaxN

#SST
#check predicted
print(temp.effect, n = Inf)

#Calculate relative change between two points 
relative_change_temp <- 0.09 / 0.04 #from low to high SST
relative_change_temp #2.25 (>1 = increase)

#calculating percentage increase in ray MaxN 
(relative_change_temp-1) * 100 #125% from low to high SST (24 C to 33 C)

#Reef relief
#check predicted
print(relief.effect.Ham, n = Inf)

#Calculate relative change between two points 
relative_change_HamRelief <- 0.03 / 0.09 #from low to high relief
relative_change_HamRelief #0.33 (<1 = decrease)

#calculating percentage decrease in ray MaxN 
(1-relative_change_HamRelief) * 100 #66.7% from low to high reef complexity

#Bottom Gear
#check predicted
print(gear.effect, n = Inf)

#Calculate relative change between two points 
relative_change_gear <- 0.04 / 0.06 #from no gear to gear
relative_change_gear #0.67 (<1 = decrease)

#calculating percentage decrease in ray MaxN 
(1-relative_change_gear) * 100 #33.3% from no gear to gear 

#R squared for Ham
require(performance)
r2(HamBest) 
# Marginal is the variance explained by the fixed effects (0.25), 
# and conditional is variance explained by fixed & random effects together (0.77)

#To calculate variance explained for a particular variable run the model without the
#variable and compare the marginal R squared:

HamNoShark<- glmmTMB(Ham_maxn ~ geomorphology + temp + relief +  
                                bottom_gear + (1|country),
                                family="tweedie", 
                                data = raysnoBR)
r2(HamNoShark)
# Marginal R squared is reduced from 0.250 to 0.199 when you remove sharks, implying that 
# ~ 5% of the total variance is explained by the shark variable, when
# other variables are included. 
(0.250-0.199)*100 #5.1%

#Interestingly, leaving out other variables, and including only sharks
#only 4% of the variance is explained. That is roughly the same ballpark. Because
#the predictor variables interact, we don't expect to get exactly the same number. 
HamOnlyShark<- glmmTMB(Ham_maxn ~ shark_maxn + (1|country),
                                  family="tweedie", 
                                  data = raysnoBR)
r2(HamOnlyShark) #0.040 marginal (4%)

#RELIEF
HamNoRelief<- glmmTMB(Ham_maxn ~ geomorphology + temp + shark_maxn +  
                                 bottom_gear + (1|country),
                                 family="tweedie", 
                                 data = raysnoBR)
r2(HamNoRelief)
(0.250 - 0.150)*100 #10%

#GEOMORPHOLOGY
HamNoGeo<- glmmTMB(Ham_maxn ~ relief + temp + shark_maxn +  
                              bottom_gear + (1|country),
                              family="tweedie", 
                              data = raysnoBR)
r2(HamNoGeo) 
(0.250-0.188) * 100 #6.2%

#TEMPERATURE
HamNoTemp<- glmmTMB(Ham_maxn ~ geomorphology + relief + shark_maxn +  
                               bottom_gear + (1|country),
                               family="tweedie", 
                               data = raysnoBR)
r2(HamNoTemp) 
(0.250-0.226)*100 #2.4%

#BOTTOM GEAR
HamNoGear<- glmmTMB(Ham_maxn ~ geomorphology + relief + shark_maxn +  
                               temp + (1|country),
                               family="tweedie", 
                               data = raysnoBR)
r2(HamNoGear)
(0.250-0.287)*100 #-3.7% (incl. gear reduces variance explained)

HamOnlyGear<-glmmTMB(Ham_maxn ~ bottom_gear + (1|country), family = "tweedie",
                                data = raysnoBR)
r2(HamOnlyGear) #3.4% - small independent effect but may be less insightful with other variables included in model
