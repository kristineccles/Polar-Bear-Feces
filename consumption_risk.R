##########################################################################
# Human exposure application
# Written By: Kristin Eccles
# Date: September 22nd, 2022
# Updated: May 10th, 2024
# Note: Must run metals_pred_models first
# Note: Figure 4 is made in this script
##########################################################################

library(truncnorm)
library(ggplot2)
# set seeed
set.seed(3432)

#### Estimate Intake ####
# Laird 2009 IHS polar bear meat consumption
# Polar bear meat 9.7 ± 68 Intake, g/wk
# measured Hg 143 ± 113

# make synthetic population of n =1000
#calcualte grams per day of polar bear meat consumption - divide by 7 to estimate daily intake
intake_polar_bear_meat <- rtruncnorm(n = 1000, a=0, b=Inf, mean = (9.7), sd = (68))
intake_polar_bear_meat <- as.data.frame(intake_polar_bear_meat)
# units are g/week

#plot
ggplot(data = intake_polar_bear_meat, aes(x= intake_polar_bear_meat))+
  geom_histogram(fill= " grey")+
  theme_bw()+
  labs(y="Count", x = "Intake Polar Bear Meat g/week")

#### Risk MeHg Muscle ####
# Measured
#get mean and 95% CI of mean
t.test(metals_log$muscle_MeHg, conf.level = 0.95)
#mean of x = -0.4658054 , 95 percent confidence interval: -0.5750938 -0.3565171

#max
10^(max(metals_log$muscle_MeHg, na.rm = TRUE))
# 1.533459 ug/g
# upper 95% of consumption: 9.7 + 68 = 77.7 g/week
# convert to day 11.1 g/day
# 11.1 g/day * 1.533459 ug/g = 17.02139 ug/day
# convert to mg 0.01702139 mg/day

# Predicted muscle MeHg from Feces
#predicted 95% CI
#muscle	(Intercept)	-0.200092237	0.068160702
#muscle	Feces_MeHg	0.389689945	0.085754093

# Use the 95% upper mean of the predicted mean
metals_log$predict_muscle <- (-0.200092237) + (metals_log$Feces_MeHg)* (0.085754093)
t.test(metals_log$predict_muscle, conf.level = 0.95)
#mean of x  = -0.2598257  , 95 percent confidence interval: -0.2701847 -0.2494667
#units are ug/g

#calculate intake of mehg based on consumption, g/w * ug/g, unit = ug/week
# divide by 67 for an average weight person
intake_polar_bear_meat$MeHg_intake_pred <- (intake_polar_bear_meat$intake_polar_bear_meat * ((10^(-0.2494667)/67)))
intake_polar_bear_meat$MeHg_intake_measured <- (intake_polar_bear_meat$intake_polar_bear_meat * ((10^(-0.3565171)/67)))
quantile(intake_polar_bear_meat$MeHg_intake_pred, probs = seq(0.95))

#summarize
max(intake_polar_bear_meat$MeHg_intake_pred)
max(intake_polar_bear_meat$MeHg_intake_measured)

stack_MeHg_muscle <- stack(intake_polar_bear_meat[,2:3])

fig4 <- ggplot(data = stack_MeHg_muscle, aes(x= values, fill = ind))+
  geom_density(alpha=0.5)+
  
  geom_vline(xintercept=(0.1*7), linetype="dashed", color = "black")+
  geom_text(aes(x=0.1*7, y=1.25), label="\nEPA BMDL5 + UF", colour="black", angle=90) +
  
  # IRIS EPA BMDL5 is 0.6 x 10-4 mg/kg/day = 6.02 µg/kg/week 
  geom_vline(xintercept=(6.02), linetype="dashed", color = "black")+
  geom_text(aes(x=6.02, y=1.25), label="\nEPA BMDL5", colour="black", angle=90) +

  geom_vline(xintercept=(0.2*7), linetype="dashed", color = "black")+
  geom_text(aes(x=0.2*7, y=1.25), label="\nHC TDI - Sensitive", colour="black", angle=90) +
  
  geom_vline(xintercept=(0.47*7), linetype="dashed", color = "black")+
  geom_text(aes(x=0.47*7, y=1.25), label="\nHC TDI", colour="black", angle=90) +
  
  theme_bw()+
  
  scale_fill_discrete(labels = c("Predicted Intake", "Measured Intake"))+
  labs(fill = "", x = "MeHg Intake from Polar Bear Meat (ug/kg/week)", y = "Denisty")

fig4
cowplot::save_plot("figure4.jpg", fig4, base_width = 6, base_height = 5.5, dpi = 300)

# population exceeding consumption advisory
#HCTRV
sum(intake_polar_bear_meat$MeHg_intake_measured > (0.2*7), na.rm=TRUE)/1000*100

sum(intake_polar_bear_meat$MeHg_intake_measured > (0.1*7), na.rm=TRUE)/1000*100

#### MeHg Risk Liver ####
# Measured
#get mean and 95% CI of mean
t.test(metals_log$liver_MeHg)
#mean of x = 0.2351122, 95 percent confidence interval:0.1421387 0.3280856

# Predicted
#predicted 95% CI
#liver	(Intercept)	0.572220912	0.072302896
#liver	feaces_MeHg	0.491937108	0.089360316

# Use the 95% upper mean of the predicted mean
metals_log$predict_muscle <- (0.572220912) + (metals_log$feaces_MeHg)* (0.491937108+	0.089360316)
t.test(metals_log$predict_muscle, conf.level = 0.95)
#mean of x  = 0.1673086  , 95 percent confidence interval: 0.23752856 

#calculate intake of mehg based on exposure
intake_polar_bear_meat$MeHg_intake_pred <- intake_polar_bear_meat$intake_polar_bear_meat * 10^(0.23752856) *0.001
intake_polar_bear_meat$MeHg_intake_measured <- intake_polar_bear_meat$intake_polar_bear_meat * 10^(0.23752856) *0.001

ggplot()+
  geom_density(data = intake_polar_bear_meat, aes(x= MeHg_intake_pred, fill = "red"), alpha=0.5)+
  geom_density(data = intake_polar_bear_meat, aes(x= MeHg_intake_measured, fill = "blue"), alpha=0.5)+
  #geom_vline(xintercept=0.0688, linetype="dashed", color = "red")+
  theme_bw()+
  labs(x = "MeHg Intake Polar Bear Meat g/day")

#####################################################################################################
#### Hg ####

#hg mean ug/g  = 0.5911122, sd = 0.4447227
t.test(metals_log$muscle_Hg)
#0.5911122 95% CI: 0.4633730 0.7188515

intake_polar_bear_meat$Hg <- intake_polar_bear_meat$intake_polar_bear_meat * 0.7188515 *0.001
#mg/day
#hg reference value 0.0001 mg/kg-d
#average weight is 80kg
# 0.0001* 80 = 0.008 mg/day
quantile(intake_polar_bear_meat$Hg, probs = c(.05, .5, .95))

ggplot(data = intake_polar_bear_meat, aes(x= Hg))+
  geom_histogram(fill= " grey")+
  theme_bw()+
  geom_vline(xintercept=0.008, linetype="dashed", color = "red")+
  labs(y="Count", x = "Intake Polar Bear Meat g/day")


#intercept -0.369436113	0.037731635
# beta 0.272595203	0.081394562

#Mean feces 1.98 (1.73)  
