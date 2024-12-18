# Playback analysis master script for publication
###
#GAM detection functions and calculation of EDR & EDA
# Nicole Todd
# Email: nicole.t@live.co.uk
# created March 2024



######### Detection functions ######
## Estimation of the detection function
# Detection function to be calculated for each device
# based on Kyhn et al., 2012, Nuuttila et al., 2018, and Amundin et al., 2022

#packages required
install.packages(mgcv)
#GAMs used to estimate the detection function for the detected playback
library(mgcv)

#load in playback data
pb_sher <- read.csv ("playback2.2.csv")
#latest dataframe is 2.2 

#fix some variables that may be of interest
pb_sher$recording <- as.factor(pb_sher$recording)
pb_sher$track <- as.factor(pb_sher$track)
pb_sher$transect <- as.factor(pb_sher$transect)
pb_sher$gain <- as.numeric(pb_sher$gain)
pb_sher$id <- as.factor(pb_sher$id)

#### CPOD ####
#starting model with only distance
model_c_dist <- gam(detected_cpod_raw ~ s(true_distance, k=5), data= pb_sher,
               family= binomial(link = logit))
summary(model_c_dist)
AIC(model_c_dist)
plot(model_c_dist, rug = TRUE)


#final CPOD model 
#rerun 23.09
model_c2 <- gam(detected_cpod_raw ~ s(true_distance, k=5)+ transect +SL2,
             data= pb_sher,family= binomial(link = logit))
summary(model_c2) 
AIC(model_c2)#60.11972
#
## recording subset ###
library(dplyr)
pb_sher$recording <- as.factor(pb_sher$recording)
#pb_sher_rec <- pb_sher %>%
 # filter(recording =='rec 3', 'rec 5')

pb_sher_rec <- subset(pb_sher, recording == 'rec 3' | recording == 'rec 5')


model_c_rec <- gam(detected_cpod_raw ~ s(true_distance, k=5) +recording,
                   data= pb_sher_rec,family= binomial(link = logit))
summary(model_c_rec) 
#no significant difference between clicks or buzzes for CPOD

##### FPOD ####
#starting model with only distance
model_f_dist <- gam(detected_fpod_raw ~ s(true_distance, k=5), data= pb_sher,
               family= binomial(link = logit))
summary(model_f_dist)
AIC(model_f_dist)
plot(model_f_dist, rug = TRUE)


#final model with ammended SL
#incl. dist and depth but only dist. significant
model_f <- gam(detected_fpod_raw ~ s(true_distance, k=5) + s(depth, k=5),
             data= pb_sher,family= binomial(link = logit))
summary(model_f) 
AIC(model_f) #90.83254
#23.09 same model as SL not sig

## recording subset ###
model_f_rec <- gam(detected_fpod_raw ~ s(true_distance, k=5) +recording,
                   data= pb_sher_rec,family= binomial(link = logit))
summary(model_f_rec) 
#no significant difference


#### SoundTrap ####
#starting model with only distance
model_st_dist <- gam(detected_st ~ s(true_distance, k=5), data= pb_sher,
               family= binomial(link = logit))
summary(model_st_dist)
AIC(model_st_dist)
plot(model_st_dist, rug = TRUE)

#new final model

model_st2 <- gam(detected_st ~ s(true_distance, k=5) + s(depth, k=5)+ transect +SL2,
             data= pb_sher,family= binomial(link = logit))
summary(model) 
AIC(model) #116.9709


## recording subset ###
model_st_rec <- gam(detected_st ~ s(true_distance, k=5) +recording,
                    data= pb_sher_rec,family= binomial(link = logit))
summary(model_st_rec) 
plot.gam(model_st_rec, all.terms = T)
#significant difference

##### make predictions from chosen GAM ####
#no need to set source level as assuming the one level as was one experiment
#predict with more variables
#change variables according to what was included in the model

#e.g. for CPOD
newdata <- data.frame(true_distance = 0:500,
                      #depth=median(pb_sher$depth),
                      transect=median(pb_sher$transect),
                      SL2= median(pb_sher$SL2))

#for recording subset
#newdata <- data.frame(true_distance = 0:500,
                    #  recording= "rec 3")

# Ensure factor levels match

ob1 <- predict(model_c, newdata = newdata, se.fit=TRUE) #change model as needed
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
pred_1$distance <- 1:nrow(pred_1)

colnames(pred_1)[1] <- "distance"
write.csv(pred_1, "C_pred.csv")#save file 

#same plot using ggplot for customization
library(ggplot2)
CPOD <-  ggplot(pred_1, aes(x = distance, y = pred)) +
  geom_line()+
  geom_line(colour = "black", size = 0.5, alpha = 0.5)+ 
  geom_line(aes(x = distance, y = ucl), colour = "#35b779", size = 1, linetype = 2) +
  geom_line(aes(x = distance, y = lcl), colour = "#35b779", size = 1, linetype = 2)+
  labs(x = "Distance from C-POD (m)",
       y = "Detection probability")+
  theme_classic()+
  ggtitle('C')
CPOD

###### EDR  calc #####

C_pred <- read.csv( "C_pred.csv") #predicted values from predict.gam
#detection probability = g(y) detection function for any given distance y
#For written formula see Kyhn 2012
###### integrate out over predictions and distances#####

fit <- smooth.spline(newdata$true_distance, pred)
smooth <- function(x) x*predict(fit, x)$y
integrate(smooth, 0, 500)
#remember to use integration result for next step
p_est<- (2/250000)*24039.73             #2/w2 * integration result
sqrt(p_est* 250000) # equals EDR

fit <- smooth.spline(newdata$true_distance, lcl)#### compute CI, ucl & lcl
smooth <- function(x) x*predict(fit, x)$y
integrate(smooth, 0, 500)
#remember to use integration result for next step
p_est<- (2/250000)*20058.6               #2/w2 * integration result
sqrt(p_est* 250000)

fit <- smooth.spline(newdata$true_distance, ucl)#### compute CI, ucl & lcl
smooth <- function(x) x*predict(fit, x)$y
integrate(smooth, 0, 500)
#remember to use integration result for next step
p_est<- (2/250000)*29507.63                   #2/w2 * integration result
sqrt(p_est* 250000)

######looking at EDA######
#compute integrate function first, end step different


eda <-(2*pi)*29507.63  #(2*pi)*integration result
eda/ (100*100) #eda in ha
eda/1000000 #eda in km2

#######
