# Playback analysis master script for publication
###
#GAM detection functions and calculation of EDR & EDA
# Nicole Todd
# Email: nicole.t@live.co.uk
# March 2024



######### Detection functions ######
#packages required
install.packages(mgcv)
#GAMs used to estimate the detection function for the detected playback
library(mgcv)

#load in playback data
pb_sher <- read.csv ("playback2.csv")

#fix some variables that may be of interest
pb_sher$recording <- as.factor(pb_sher$recording)
pb_sher$track <- as.factor(pb_sher$track)
pb_sher$transect <- as.factor(pb_sher$transect)
pb_sher$gain <- as.numeric(pb_sher$gain)
pb_sher$id <- as.factor(pb_sher$id)

### CPOD ##
#starting model with only distance
model_c_dist <- gam(detected_cpod_raw ~ s(true_distance, k=5), data= pb_sher,
               family= binomial(link = logit))
summary(model_c_dist)
AIC(model_c_dist)
plot(model_c_dist, rug = TRUE)

#Final CPOD model included distance, SL and transect
model <- gam(detected_cpod_raw ~ s(true_distance, k=5)+ transect +source.level,
             data= pb_sher,family= binomial(link = logit))
summary(model) 
AIC(model)
plot.gam(model, all.terms = T)
anova(model)

## recording subset ###
pb_sher$recording <- as.factor(pb_sher$recording)
pb_sher_rec <- pb_sher %>%
  filter(recording =='rec 3', 'rec 5')

pb_sher_rec <- subset(pb_sher, recording == 'rec 3' | recording == 'rec 5')

model_c_rec <- gam(detected_cpod_raw ~ s(true_distance, k=5) +recording,
                   data= pb_sher_rec,family= binomial(link = logit))
summary(model_c_rec) 

## FPOD ##
#starting model with only distance
model_f_dist <- gam(detected_fpod_raw ~ s(true_distance, k=5), data= pb_sher,
               family= binomial(link = logit))
summary(model_f_dist)
AIC(model_f_dist)
plot(model_f_dist, rug = TRUE)

#Final FPOD model included distance, depth and SL
model <- gam(detected_fpod_raw ~ s(true_distance, k=5) + s(depth, k=5) + s(source.level),
             data= pb_sher,family= binomial(link = logit))
summary(model) 
AIC(model)
## recording subset ###
model_f_rec <- gam(detected_fpod_raw ~ s(true_distance, k=5) +recording,
                   data= pb_sher_rec,family= binomial(link = logit))
summary(model_f_rec) 

## SoundTrap ##
#starting model with only distance
model_st_dist <- gam(detected_st ~ s(true_distance, k=5), data= pb_sher,
               family= binomial(link = logit))
summary(model_st_dist)
AIC(model_st_dist)
plot(model_st_dist, rug = TRUE)

#Final ST model included distance, depth, SL and transect
model <- gam(detected_st ~ s(true_distance, k=5) + s(depth, k=5)+ transect +source.level,
             data= pb_sher,family= binomial(link = logit))
summary(model) 
AIC(model)
## recording subset ###
model_st_rec <- gam(detected_st ~ s(true_distance, k=5) +recording,
                    data= pb_sher_rec,family= binomial(link = logit))
summary(model_st_rec) 
plot.gam(model_st_rec, all.terms = T)


##### make predictions from chosen GAM ####
#no need to set source level as assuming the one level as was one experiment
#predict with more variables
#change variables according to what was included in the model

#e.g. for CPOD
newdata <- data.frame(true_distance = 0:500,
                      depth=median(pb_sher$depth),
                      #transect=median(pb_sher$transect),
                      source.level= median(pb_sher$source.level))

#for recording subset
#newdata <- data.frame(true_distance = 0:500,
                    #  recording= "rec 3")

ob1 <- predict(model, newdata = newdata, se.fit=TRUE) #change model as needed
require(boot)
pred<-inv.logit(ob1$fit)
#confidence intervals
lcl<-inv.logit(ob1$fit-1.96*ob1$se.fit)
#plot detection probability
ucl<-inv.logit(ob1$fit+1.96*ob1$se.fit)
plot(newdata$true_distance,pred,type="l", xlab="Distance (m)", ylab="Detection probability", ylim=c(0,1))+
  lines(newdata$true_distance,lcl,lty=2)+
  lines(newdata$true_distance,ucl,lty=2)+
  rug(pb_sher$true_distance)+
  title('C-POD detection probability')

pred_1<- data.frame(pred) #saving predicted values to data frame for EDR estimation
write.csv(pred_1, "C_pred.csv")#save file 


###### EDR  calc #####

C_pred <- read.csv( "C_pred.csv") #predicted values from predict.gam
#detection probability = g(y) detection fucntion for any given distance y
#For written formula see Kyhn 2012
###### integrate out over predictions and distances#####

fit <- smooth.spline(newdata$true_distance, pred)
smooth <- function(x) x*predict(fit, x)$y
integrate(smooth, 0, 500)
p_est<- (2/250000)*43709.03            #2/w2 * integration result
sqrt(p_est* 250000) # equals EDR

fit <- smooth.spline(newdata$true_distance, lcl)#### compute CI, ucl & lcl
smooth <- function(x) x*predict(fit, x)$y
integrate(smooth, 0, 500)
p_est<- (2/250000)*30845.7                         #2/w2 * integration result
sqrt(p_est* 250000)

fit <- smooth.spline(newdata$true_distance, ucl)#### compute CI, ucl & lcl
smooth <- function(x) x*predict(fit, x)$y
integrate(smooth, 0, 500)
p_est<- (2/250000)*61095.64                        #2/w2 * integration result
sqrt(p_est* 250000)

######looking at EDA######
#compute integrate function first, end step different


eda <-(2*pi)*28828.69  #(2*pi)*integration result
eda/ (100*100) #eda in ha
eda/1000000 #eda in km2

#######
