##########################################################################
# Predictive models using LOOCV
# Written By: Kristin Eccles
# Date: September 22nd, 2022
# Updated: May 10th, 2024
# Note: Figure 1 and 2 is made from this script
##########################################################################

# Load libraries
library(readxl)
library(rlist)
library(corrplot)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(sjPlot)
library(broom)
library(caret)
library(readr)
library(dplyr)

# Load data
read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

#mysheets <- read_excel_allsheets("PCA_tables_edit2.xlsx")
mysheets <- read_excel_allsheets("orig_data/PCA_tables.xlsx")

##########################################################################
#### Prepare Data ####
Metals_liver <- mysheets$`Metals Liver`
Metals_liver <- as.data.frame(lapply(Metals_liver[-1], as.numeric))
colnames(Metals_liver) <- paste("liver", colnames(Metals_liver), sep = "_")

Metals_fat <- mysheets$`Metals Fat`
Metals_fat <- as.data.frame(lapply(Metals_fat[-1], as.numeric))
colnames(Metals_fat) <- paste("fat", colnames(Metals_fat), sep = "_")

Metals_Feces <- as.data.frame(mysheets$`Metals Feces`)
Metals_Feces <- as.data.frame(lapply(Metals_Feces[-1], as.numeric))
colnames(Metals_Feces) <- paste("Feces", colnames(Metals_Feces), sep = "_")

Metals_muscle <- mysheets$`Metals Muscle`
Metals_muscle <- as.data.frame(lapply(Metals_muscle[-1], as.numeric))
colnames(Metals_muscle) <- paste("muscle", colnames(Metals_muscle), sep = "_")

metals_df <- cbind(Metals_Feces, Metals_liver,Metals_fat, Metals_muscle)
metals_log <- log10(metals_df)

contam_list <- as.data.frame(colnames(metals_log))
colnames(contam_list) <- "metals"
contam_list <- sub(".*_","",contam_list$metals)
contam_list <- unique(contam_list)

write.csv(metals_df, "metals_df.csv")
#### EDA ####
metals <- cor(metals_df, use = "pairwise.complete")
write.csv(metals, "corrmatrix_metals.csv")

testRes = cor.mtest(metals, conf.level = 0.95)
corrplot(metals[1:32,33:(ncol(metals))], tl.col = 'black', 
         p.mat = testRes$p[1:32,33:(ncol(testRes$p))], insig='blank')

#### Regression Model ####
ctrl <- trainControl(method = "boot", number = 1000)

#### Hg ####
# Muscle
Hg_muscle<-lm(muscle_Hg ~ Feces_Hg, data=metals_log)
Hg_muscle_coeff <- as.data.frame(tidy(Hg_muscle))
Hg_muscle_cv <- train(muscle_Hg ~ Feces_Hg, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Hg_muscle_loocv<- as.data.frame(Hg_muscle_cv$results)
sample_n <- summary(Hg_muscle_cv)$df[2]
Hg_muscle_output <- cbind("muscle",Hg_muscle_coeff , Hg_muscle_loocv, sample_n)
colnames(Hg_muscle_output)[1]<-"tissue"

#Liver
Hg_liver<-lm(liver_Hg ~ Feces_Hg, data=metals_log)
Hg_liver_coeff <- as.data.frame(tidy(Hg_liver))
Hg_liver_cv <- train(liver_Hg ~ Feces_Hg, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Hg_liver_loocv<- as.data.frame(Hg_liver_cv$results)
sample_n <- summary(Hg_liver_cv)$df[2]
Hg_liver_output <- cbind("liver",Hg_liver_coeff, Hg_liver_loocv, sample_n)
colnames(Hg_liver_output)[1]<-"tissue"

#fat
Hg_fat<-lm(fat_Hg ~ Feces_Hg, data=metals_log)
Hg_fat_coeff <- as.data.frame(tidy(Hg_fat))
Hg_fat_cv <- train(fat_Hg ~ Feces_Hg, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Hg_fat_loocv<- as.data.frame(Hg_fat_cv$results)
sample_n <- summary(Hg_fat_cv)$df[2]
Hg_fat_output <- cbind("fat", Hg_fat_coeff, Hg_fat_loocv, sample_n)
colnames(Hg_fat_output)[1]<-"tissue"

#write all outputs to csv
Hg_metal_outputs <- rbind (Hg_liver_output, Hg_muscle_output, Hg_fat_output)
write.csv(Hg_metal_outputs, "feces/metals/Hg_metal_outputs.csv")

#make plots 
Hg_fat <- ggplot(data=metals_log, aes(Feces_Hg, fat_Hg))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  #geom_hline(yintercept=log10(1), linetype='dotted')+
  labs(x="Feces Hg (ug/g)", y="Fat Hg (ug/g)")
Hg_fat

muscle_Hg <- ggplot(data=metals_log, aes(Feces_Hg, muscle_Hg))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  #geom_hline(yintercept=log10(1), linetype='dotted')+
  labs(x="Feces Hg (ug/g)", y="Muscle Hg (ug/g)")
muscle_Hg

liver_Hg <- ggplot(data=metals_log, aes(Feces_Hg, liver_Hg))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  #geom_hline(yintercept=log10(1), linetype='dotted')+
  labs(x="Feces Hg (ug/g)", y="Liver Hg (ug/g)")
liver_Hg

Hg_figure=ggarrange(liver_Hg, muscle_Hg,Hg_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
Hg_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_Hg.tif", Hg_figure, width = 20, height = 10, dpi = 300)

#### Pb ####
# Muscle
Pb_muscle<-lm(muscle_Pb ~ Feces_Pb, data=metals_log)
summary(Pb_muscle)
Pb_muscle_coeff <- as.data.frame(tidy(Pb_muscle))
Pb_muscle_cv <- train(muscle_Pb ~ Feces_Pb, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Pb_muscle_loocv<- as.data.frame(Pb_muscle_cv$results)
sample_n <- summary(Pb_muscle_cv)$df[2]
Pb_muscle_output <- cbind("muscle",Pb_muscle_coeff, Pb_muscle_loocv, sample_n)
colnames(Pb_muscle_output)[1]<-"tissue"

#Liver
Pb_liver<-lm(liver_Pb ~ Feces_Pb, data=metals_log)
Pb_liver_coeff <- as.data.frame(tidy(Pb_liver))
Pb_liver_cv <- train(liver_Pb ~ Feces_Pb, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Pb_liver_loocv<- as.data.frame(Pb_liver_cv$results)
sample_n <- summary(Pb_liver_cv)$df[2]
Pb_liver_output <- cbind("liver",Pb_liver_coeff, Pb_liver_loocv, sample_n)
colnames(Pb_liver_output)[1]<-"tissue"

#fat
Pb_fat<-lm(fat_Pb ~ Feces_Pb, data=metals_log)
Pb_fat_coeff <- as.data.frame(tidy(Pb_fat))
Pb_fat_cv <- train(fat_Pb ~ Feces_Pb, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Pb_fat_loocv<- as.data.frame(Pb_fat_cv$results)
sample_n<- summary(Pb_fat_cv)$df[2]
Pb_fat_output <- cbind("fat", Pb_fat_coeff, Pb_fat_loocv, sample_n)
colnames(Pb_fat_output)[1]<-"tissue"

#write all outputs to csv
Pb_metal_outputs <- rbind (Pb_liver_output, Pb_muscle_output,Pb_fat_output)
write.csv(Pb_metal_outputs, "feces/metals/Pb_metal_outputs.csv")

#make plots 
Pb_fat <- ggplot(data=metals_log, aes(Feces_Pb, fat_Pb))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Pb (ug/g)", y="Fat Pb (ug/g)")

muscle_Pb <- ggplot(data=metals_log, aes(Feces_Pb, muscle_Pb))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Pb (ug/g)", y="Muscle Pb (ug/g)")
muscle_Pb

liver_Pb <- ggplot(data=metals_log, aes(Feces_Pb, liver_Pb))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Pb (ug/g)", y="Liver Pb (ug/g)")
liver_Pb

Pb_figure=ggarrange(liver_Pb, muscle_Pb,Pb_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
Pb_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_Pb.tif", Pb_figure, width = 20, height = 10, dpi = 300)

#### Cd ####
# Muscle
Cd_muscle<-lm(muscle_Cd ~ Feces_Cd, data=metals_log)
summary(Cd_muscle)
Cd_muscle_coeff <- as.data.frame(tidy(Cd_muscle))
Cd_muscle_cv <- train(muscle_Cd ~ Feces_Cd, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Cd_muscle_loocv<- as.data.frame(Cd_muscle_cv$results)
sample_n <- summary(Cd_muscle_cv)$df[2]
Cd_muscle_output <- cbind("muscle",Cd_muscle_coeff, Cd_muscle_loocv, sample_n)
colnames(Cd_muscle_output)[1]<-"tissue"

#Liver
Cd_liver<-lm(liver_Cd ~ Feces_Cd, data=metals_log)
Cd_liver_coeff <- as.data.frame(tidy(Cd_liver))
Cd_liver_cv <- train(liver_Cd ~ Feces_Cd, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Cd_liver_loocv<- as.data.frame(Cd_liver_cv$results)
sample_n <- summary(Cd_liver_cv)$df[2]
Cd_liver_output <- cbind("liver",Cd_liver_coeff, Cd_liver_loocv, sample_n)
colnames(Cd_liver_output)[1]<-"tissue"

#fat
Cd_fat<-lm(fat_Cd ~ Feces_Cd, data=metals_log)

Cd_fat_coeff <- as.data.frame(tidy(Cd_fat))
Cd_fat_cv <- train(fat_Cd ~ Feces_Cd, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Cd_fat_loocv<- as.data.frame(Cd_fat_cv$results)
sample_n <- summary(Cd_fat_cv)$df[2]
Cd_fat_output <- cbind("fat", Cd_fat_coeff, Cd_fat_loocv, sample_n)
colnames(Cd_fat_output)[1]<-"tissue"

#write all outputs to csv
Cd_metal_outputs <- rbind (Cd_liver_output, Cd_muscle_output,Cd_fat_output)
write.csv(Cd_metal_outputs, "feces/metals/Cd_metal_outputs.csv")

#make plots 
Cd_fat <- ggplot(data=metals_log, aes(Feces_Cd, fat_Cd))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Cd (ug/g)", y="Fat Cd (ug/g)")

muscle_Cd <- ggplot(data=metals_log, aes(Feces_Cd, muscle_Cd))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Cd (ug/g)", y="Muscle Cd (ug/g)")
muscle_Cd

liver_Cd <- ggplot(data=metals_log, aes(Feces_Cd, liver_Cd))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Cd (ug/g)", y="Liver Cd (ug/g)")
liver_Cd

Cd_figure=ggarrange(liver_Cd, muscle_Cd,Cd_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
Cd_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_Cd.tif", Cd_figure, width = 20, height = 10, dpi = 300)

#### V ####
# Muscle
V_muscle<-lm(muscle_V ~ Feces_V, data=metals_log)
summary(V_muscle)
V_muscle_coeff <- as.data.frame(tidy(V_muscle))
V_muscle_cv <- train(muscle_V ~ Feces_V, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
V_muscle_loocv<- as.data.frame(V_muscle_cv$results)
sample_n <- summary(V_muscle_cv)$df[2]
V_muscle_output <- cbind("muscle",V_muscle_coeff, V_muscle_loocv, sample_n)
colnames(V_muscle_output)[1]<-"tissue"

#Liver
V_liver<-lm(liver_V ~ Feces_V, data=metals_log)
V_liver_coeff <- as.data.frame(tidy(V_liver))
V_liver_cv <- train(liver_V ~ Feces_V, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
V_liver_loocv<- as.data.frame(V_liver_cv$results)
sample_n <- summary(V_liver_cv)$df[2]
V_liver_output <- cbind("liver",V_liver_coeff, V_liver_loocv, sample_n)
colnames(V_liver_output)[1]<-"tissue"

#fat
V_fat<-lm(fat_V ~ Feces_V, data=metals_log)
V_fat_coeff <- as.data.frame(tidy(V_fat))
#LOOCV
V_fat_cv <- train(fat_V ~ Feces_V, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
V_fat_loocv<- as.data.frame(V_fat_cv$results)
sample_n <- summary(V_fat_cv)$df[2]
V_fat_output <- cbind("fat", V_fat_coeff, V_fat_loocv, sample_n)
colnames(V_fat_output)[1]<-"tissue"

#write all outputs to csv
V_metal_outputs <- rbind (V_liver_output, V_muscle_output,V_fat_output)
write.csv(V_metal_outputs, "feces/metals/V_metal_outputs.csv")

#make plots 
V_fat <- ggplot(data=metals_log, aes(Feces_V, fat_V))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces V (ug/g)", y="Fat V")

muscle_V <- ggplot(data=metals_log, aes(Feces_V, muscle_V))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces V (ug/g)", y="Muscle V")
muscle_V

liver_V <- ggplot(data=metals_log, aes(Feces_V, liver_V))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces V (ug/g)", y="Liver V")
liver_V

V_figure=ggarrange(liver_V, muscle_V,V_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
V_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_V.tif", V_figure, width = 20, height = 10, dpi = 300)

#### MeHg ####
# Muscle
MeHg_muscle<-lm(muscle_MeHg ~ Feces_MeHg, data=metals_log)
summary(MeHg_muscle)
MeHg_muscle_coeff <- as.data.frame(tidy(MeHg_muscle))
MeHg_muscle_cv <- train(muscle_MeHg ~ Feces_MeHg, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
MeHg_muscle_loocv<- as.data.frame(MeHg_muscle_cv$results)
sample_n <- summary(MeHg_muscle_cv)$df[2]
MeHg_muscle_output <- cbind("muscle",MeHg_muscle_coeff, MeHg_muscle_loocv, sample_n)
colnames(MeHg_muscle_output)[1]<-"tissue"

#Liver
MeHg_liver<-lm(liver_MeHg ~ Feces_MeHg, data=metals_log)
MeHg_liver_coeff <- as.data.frame(tidy(MeHg_liver))
MeHg_liver_cv <- train(liver_MeHg ~ Feces_MeHg, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
MeHg_liver_loocv<- as.data.frame(MeHg_liver_cv$results)
sample_n <- summary(MeHg_liver_cv)$df[2]
MeHg_liver_output <- cbind("liver",MeHg_liver_coeff, MeHg_liver_loocv, sample_n)
colnames(MeHg_liver_output)[1]<-"tissue"

#fat
MeHg_fat<-lm(fat_MeHg ~ Feces_MeHg, data=metals_log)
MeHg_fat_coeff <- as.data.frame(tidy(MeHg_fat))
MeHg_fat_cv <- train(fat_MeHg ~ Feces_MeHg, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
MeHg_fat_loocv<- as.data.frame(MeHg_fat_cv$results)
sample_n <- summary(MeHg_fat_cv)$df[2]
MeHg_fat_output <- cbind("fat", MeHg_fat_coeff, MeHg_fat_loocv, sample_n)
colnames(MeHg_fat_output)[1]<-"tissue"

#write all outputs to csv
MeHg_metal_outputs <- rbind (MeHg_liver_output, MeHg_muscle_output,MeHg_fat_output)
write.csv(MeHg_metal_outputs, "feces/metals/MeHg_metal_outputs.csv")

#make plots 
MeHg_fat <- ggplot(data=metals_log, aes(Feces_MeHg, fat_MeHg))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Log10 Feces MeHg (ug/g)", y="Log10 Fat MeHg (ug/g)")

muscle_MeHg <- ggplot(data=metals_log, aes(Feces_MeHg, muscle_MeHg))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces MeHg (ug/g)", y="Muscle MeHg (ug/g)")
muscle_MeHg

liver_MeHg <- ggplot(data=metals_log, aes(Feces_MeHg, liver_MeHg))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces MeHg (ug/g)", y="Liver MeHg (ug/g)")
liver_MeHg

MeHg_figure=ggarrange(liver_MeHg, muscle_MeHg,MeHg_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
MeHg_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_MeHg.tif", MeHg_figure, width = 20, height = 10, dpi = 300)


#### Se ####
# Muscle
Se_muscle<-lm(muscle_Se ~ Feces_Se, data=metals_log)
summary(Se_muscle)
Se_muscle_coeff <- as.data.frame(tidy(Se_muscle))
Se_muscle_cv <- train(muscle_Se ~ Feces_Se, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Se_muscle_loocv<- as.data.frame(Se_muscle_cv$results)
sample_n <- summary(Se_muscle_cv)$df[2]
Se_muscle_output <- cbind("muscle",Se_muscle_coeff, Se_muscle_loocv, sample_n)
colnames(Se_muscle_output)[1]<-"tissue"

#Liver
Se_liver<-lm(liver_Se ~ Feces_Se, data=metals_log)
Se_liver_coeff <- as.data.frame(tidy(Se_liver))
Se_liver_cv <- train(liver_Se ~ Feces_Se, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Se_liver_loocv<- as.data.frame(Se_liver_cv$results)
sample_n <- summary(Se_liver_cv)$df[2]
Se_liver_output <- cbind("liver",Se_liver_coeff, Se_liver_loocv, sample_n)
colnames(Se_liver_output)[1]<-"tissue"

#fat
Se_fat<-lm(fat_Se ~ Feces_Se, data=metals_log)
Se_fat_coeff <- as.data.frame(tidy(Se_fat))
Se_fat_cv <- train(fat_Se ~ Feces_Se, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Se_fat_loocv<- as.data.frame(Se_fat_cv$results)
sample_n <- summary(Se_fat_cv)$df[2]
Se_fat_output <- cbind("fat", Se_fat_coeff, Se_fat_loocv, sample_n)
colnames(Se_fat_output)[1]<-"tissue"

#write all outputs to csv
Se_metal_outputs <- rbind (Se_liver_output, Se_muscle_output,Se_fat_output)
write.csv(Se_metal_outputs, "feces/metals/Se_metal_outputs.csv")

#make plots 
Se_fat <- ggplot(data=metals_log, aes(Feces_Se, fat_Se))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Se (ug/g)", y="Fat Se")

muscle_Se <- ggplot(data=metals_log, aes(Feces_Se, muscle_Se))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Se (ug/g)", y="Muscle Se")
muscle_Se

liver_Se <- ggplot(data=metals_log, aes(Feces_Se, liver_Se))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Se (ug/g)", y="Liver Se")
liver_Se

Se_figure=ggarrange(liver_Se, muscle_Se,Se_fat, 
                      labels = c("A", "B", "C"),
                      vjust = 1,
                      hjust = -0.5,
                      ncol = 3, nrow = 1,
                      legend = "right",
                      font.label = list(size = 16))
Se_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_Se.tif", Se_figure, width = 20, height = 10, dpi = 300)


#### As ####
# Muscle
As_muscle<-lm(muscle_As ~ Feces_As, data=metals_log)
summary(As_muscle)
As_muscle_coeff <- as.data.frame(tidy(As_muscle))
As_muscle_cv <- train(muscle_As ~ Feces_As, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
As_muscle_loocv<- as.data.frame(As_muscle_cv$results)
sample_n <- summary(As_muscle_cv)$df[2]
As_muscle_output <- cbind("muscle",As_muscle_coeff, As_muscle_loocv, sample_n)
colnames(As_muscle_output)[1]<-"tissue"

#Liver
As_liver<-lm(liver_As ~ Feces_As, data=metals_log)
As_liver_coeff <- as.data.frame(tidy(As_liver))
As_liver_cv <- train(liver_As ~ Feces_As, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
As_liver_loocv<- as.data.frame(As_liver_cv$results)
sample_n <- summary(As_liver_cv)$df[2]
As_liver_output <- cbind("liver",As_liver_coeff, As_liver_loocv, sample_n)
colnames(As_liver_output)[1]<-"tissue"

#fat
As_fat<-lm(fat_As ~ Feces_As, data=metals_log)
As_fat_coeff <- as.data.frame(tidy(As_fat))
As_fat_cv <- train(fat_As ~ Feces_As, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
As_fat_loocv<- as.data.frame(As_fat_cv$results)
sample_n <- summary(As_fat_cv)$df[2]
As_fat_output <- cbind("fat", As_fat_coeff, As_fat_loocv, sample_n)
colnames(As_fat_output)[1]<-"tissue"

#write all outputs to csv
As_metal_outputs <- rbind (As_liver_output, As_muscle_output,As_fat_output)
write.csv(As_metal_outputs, "feces/metals/As_metal_outputs.csv")

#make plots 
As_fat <- ggplot(data=metals_log, aes(Feces_As, fat_As))+
  geom_point()+
  geom_smooth(method = "lm", As = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces As (ug/g)", y="Fat As (ug/g)")

muscle_As <- ggplot(data=metals_log, aes(Feces_As, muscle_As))+
  geom_point()+
  geom_smooth(method = "lm", As = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces As (ug/g)", y="Muscle As (ug/g)")
muscle_As

liver_As <- ggplot(data=metals_log, aes(Feces_As, liver_As))+
  geom_point()+
  geom_smooth(method = "lm", As = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces As (ug/g)", y="Liver As (ug/g)")
liver_As

As_figure=ggarrange(liver_As, muscle_As,As_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
As_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_As.tif", As_figure, width = 20, height = 10, dpi = 300)


#### Zn ####
# Muscle
Zn_muscle<-lm(muscle_Zn ~ Feces_Zn, data=metals_log)
summary(Zn_muscle)
Zn_muscle_coeff <- as.data.frame(tidy(Zn_muscle))
Zn_muscle_cv <- train(muscle_Zn ~ Feces_Zn, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Zn_muscle_loocv<- as.data.frame(Zn_muscle_cv$results)
sample_n <- summary(Zn_muscle_cv)$df[2]
Zn_muscle_output <- cbind("muscle",Zn_muscle_coeff, Zn_muscle_loocv, sample_n)
colnames(Zn_muscle_output)[1]<-"tissue"

#Liver
Zn_liver<-lm(liver_Zn ~ Feces_Zn, data=metals_log)
Zn_liver_coeff <- as.data.frame(tidy(Zn_liver))
Zn_liver_cv <- train(liver_Zn ~ Feces_Zn, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Zn_liver_loocv<- as.data.frame(Zn_liver_cv$results)
sample_n <- summary(Zn_liver_cv)$df[2]
Zn_liver_output <- cbind("liver",Zn_liver_coeff, Zn_liver_loocv, sample_n)
colnames(Zn_liver_output)[1]<-"tissue"

#fat
Zn_fat<-lm(fat_Zn ~ Feces_Zn, data=metals_log)
Zn_fat_coeff <- as.data.frame(tidy(Zn_fat))
Zn_fat_cv <- train(fat_Zn ~ Feces_Zn, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Zn_fat_loocv<- as.data.frame(Zn_fat_cv$results)
sample_n <- summary(Zn_fat_cv)$df[2]
Zn_fat_output <- cbind("fat", Zn_fat_coeff, Zn_fat_loocv, sample_n)
colnames(Zn_fat_output)[1]<-"tissue"

#write all outputs to csv
Zn_metal_outputs <- rbind (Zn_liver_output, Zn_muscle_output,Zn_fat_output)
write.csv(Zn_metal_outputs, "feces/metals/Zn_metal_outputs.csv")

#make plots 
Zn_fat <- ggplot(data=metals_log, aes(Feces_Zn, fat_Zn))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Zn (ug/g)", y="Fat Zn")

muscle_Zn <- ggplot(data=metals_log, aes(Feces_Zn, muscle_Zn))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Zn (ug/g)", y="Muscle Zn")
muscle_Zn

liver_Zn <- ggplot(data=metals_log, aes(Feces_Zn, liver_Zn))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Zn (ug/g)", y="Liver Zn")
liver_Zn

Zn_figure=ggarrange(liver_Zn, muscle_Zn,Zn_fat, 
                      labels = c("A", "B", "C"),
                      vjust = 1,
                      hjust = -0.5,
                      ncol = 3, nrow = 1,
                      legend = "right",
                      font.label = list(size = 16))
Zn_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_Zn.tif", Zn_figure, width = 20, height = 10, dpi = 300)


#### K ####
# Muscle
K_muscle<-lm(muscle_K ~ Feces_K, data=metals_log)
summary(K_muscle)
K_muscle_coeff <- as.data.frame(tidy(K_muscle))
K_muscle_cv <- train(muscle_K ~ Feces_K, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
K_muscle_loocv<- as.data.frame(K_muscle_cv$results)
sample_n <- summary(K_muscle_cv)$df[2]
K_muscle_output <- cbind("muscle",K_muscle_coeff, K_muscle_loocv, sample_n)
colnames(K_muscle_output)[1]<-"tissue"

#Liver
K_liver<-lm(liver_K ~ Feces_K, data=metals_log)
K_liver_coeff <- as.data.frame(tidy(K_liver))
K_liver_cv <- train(liver_K ~ Feces_K, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
K_liver_loocv<- as.data.frame(K_liver_cv$results)
sample_n <- summary(K_liver_cv)$df[2]
K_liver_output <- cbind("liver",K_liver_coeff, K_liver_loocv, sample_n)
colnames(K_liver_output)[1]<-"tissue"

#fat
K_fat<-lm(fat_K ~ Feces_K, data=metals_log)
K_fat_coeff <- as.data.frame(tidy(K_fat))
K_fat_cv <- train(fat_K ~ Feces_K, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
K_fat_loocv<- as.data.frame(K_fat_cv$results)
sample_n <- summary(K_fat_cv)$df[2]
K_fat_output <- cbind("fat", K_fat_coeff, K_fat_loocv, sample_n)
colnames(K_fat_output)[1]<-"tissue"

#write all outputs to csv
K_metal_outputs <- rbind (K_liver_output, K_muscle_output,K_fat_output)
write.csv(K_metal_outputs, "feces/metals/K_metal_outputs.csv")

#make plots 
K_fat <- ggplot(data=metals_log, aes(Feces_K, fat_K))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces K (ug/g)", y="Fat K")

muscle_K <- ggplot(data=metals_log, aes(Feces_K, muscle_K))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces K (ug/g)", y="Muscle K")
muscle_K

liver_K <- ggplot(data=metals_log, aes(Feces_K, liver_K))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces K (ug/g)", y="Liver K")
liver_K

K_figure=ggarrange(liver_K, muscle_K,K_fat, 
                      labels = c("A", "B", "C"),
                      vjust = 1,
                      hjust = -0.5,
                      ncol = 3, nrow = 1,
                      legend = "right",
                      font.label = list(size = 16))
K_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_K.tif", K_figure, width = 20, height = 10, dpi = 300)


#### Fe ####
# Muscle
Fe_muscle<-lm(muscle_Fe ~ Feces_Fe, data=metals_log)
summary(Fe_muscle)
Fe_muscle_coeff <- as.data.frame(tidy(Fe_muscle))
Fe_muscle_cv <- train(muscle_Fe ~ Feces_Fe, data = metals_log, method = "lm", trControl = ctrl, na.acFeon=na.omit)
Fe_muscle_loocv<- as.data.frame(Fe_muscle_cv$results)
sample_n <- summary(Fe_muscle_cv)$df[2]
Fe_muscle_output <- cbind("muscle",Fe_muscle_coeff, Fe_muscle_loocv, sample_n)
colnames(Fe_muscle_output)[1]<-"tissue"

#Liver
Fe_liver<-lm(liver_Fe ~ Feces_Fe, data=metals_log)
Fe_liver_coeff <- as.data.frame(tidy(Fe_liver))
Fe_liver_cv <- train(liver_Fe ~ Feces_Fe, data = metals_log, method = "lm", trControl = ctrl, na.acFeon=na.omit)
Fe_liver_loocv<- as.data.frame(Fe_liver_cv$results)
sample_n <- summary(Fe_liver_cv)$df[2]
Fe_liver_output <- cbind("liver",Fe_liver_coeff, Fe_liver_loocv, sample_n)
colnames(Fe_liver_output)[1]<-"tissue"

#fat
Fe_fat<-lm(fat_Fe ~ Feces_Fe, data=metals_log)
Fe_fat_coeff <- as.data.frame(tidy(Fe_fat))
Fe_fat_cv <- train(fat_Fe ~ Feces_Fe, data = metals_log, method = "lm", trControl = ctrl, na.acFeon=na.omit)
Fe_fat_loocv<- as.data.frame(Fe_fat_cv$results)
sample_n <- summary(Fe_fat_cv)$df[2]
Fe_fat_output <- cbind("fat", Fe_fat_coeff, Fe_fat_loocv, sample_n)
colnames(Fe_fat_output)[1]<-"tissue"

#write all outputs to csv
Fe_metal_outputs <- rbind (Fe_liver_output, Fe_muscle_output,Fe_fat_output)
write.csv(Fe_metal_outputs, "feces/metals/Fe_metal_outputs.csv")

#make plots 
Fe_fat <- ggplot(data=metals_log, aes(Feces_Fe, fat_Fe))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Fe (ug/g)", y="Fat Fe")

muscle_Fe <- ggplot(data=metals_log, aes(Feces_Fe, muscle_Fe))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Fe (ug/g)", y="Muscle Fe")
muscle_Fe

liver_Fe <- ggplot(data=metals_log, aes(Feces_Fe, liver_Fe))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Fe (ug/g)", y="Liver Fe")
liver_Fe

Fe_figure=ggarrange(liver_Fe, muscle_Fe,Fe_fat, 
                      labels = c("A", "B", "C"),
                      vjust = 1,
                      hjust = -0.5,
                      ncol = 3, nrow = 1,
                      legend = "right",
                      font.label = list(size = 16))
Fe_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_Fe.tif", Fe_figure, width = 20, height = 10, dpi = 300)

#### Al ####
# Muscle
Al_muscle<-lm(muscle_Al ~ Feces_Al, data=metals_log)
summary(Al_muscle)
Al_muscle_coeff <- as.data.frame(tidy(Al_muscle))
Al_muscle_cv <- train(muscle_Al ~ Feces_Al, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Al_muscle_loocv<- as.data.frame(Al_muscle_cv$results)
sample_n <- summary(Al_muscle_cv)$df[2]
Al_muscle_output <- cbind("muscle",Al_muscle_coeff, Al_muscle_loocv, sample_n)
colnames(Al_muscle_output)[1]<-"tissue"

#Liver
Al_liver<-lm(liver_Al ~ Feces_Al, data=metals_log)
Al_liver_coeff <- as.data.frame(tidy(Al_liver))
Al_liver_cv <- train(liver_Al ~ Feces_Al, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Al_liver_loocv<- as.data.frame(Al_liver_cv$results)
sample_n <- summary(Al_liver_cv)$df[2]
Al_liver_output <- cbind("liver",Al_liver_coeff, Al_liver_loocv, sample_n)
colnames(Al_liver_output)[1]<-"tissue"

#fat
Al_fat<-lm(fat_Al ~ Feces_Al, data=metals_log)
Al_fat_coeff <- as.data.frame(tidy(Al_fat))
Al_fat_cv <- train(fat_Al ~ Feces_Al, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Al_fat_loocv<- as.data.frame(Al_fat_cv$results)
sample_n <- summary(Al_fat_cv)$df[2]
Al_fat_output <- cbind("fat", Al_fat_coeff, Al_fat_loocv, sample_n)
colnames(Al_fat_output)[1]<-"tissue"

#write all outputs to csv
Al_metal_outputs <- rbind (Al_liver_output, Al_muscle_output,Al_fat_output)
write.csv(Al_metal_outputs, "feces/metals/Al_metal_outputs.csv")

#make plots 
Al_fat <- ggplot(data=metals_log, aes(Feces_Al, fat_Al))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Al (ug/g)", y="Fat Al")

muscle_Al <- ggplot(data=metals_log, aes(Feces_Al, muscle_Al))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Al (ug/g)", y="Muscle Al")
muscle_Al

liver_Al <- ggplot(data=metals_log, aes(Feces_Al, liver_Al))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Al (ug/g)", y="Liver Al")
liver_Al

Al_figure=ggarrange(liver_Al, muscle_Al,Al_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
Al_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_Al.tif", Al_figure, width = 20, height = 10, dpi = 300)

#### Ba ####
# Muscle
Ba_muscle<-lm(muscle_Ba ~ Feces_Ba, data=metals_log)
summary(Ba_muscle)
Ba_muscle_coeff <- as.data.frame(tidy(Ba_muscle))
Ba_muscle_cv <- train(muscle_Ba ~ Feces_Ba, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Ba_muscle_loocv<- as.data.frame(Ba_muscle_cv$results)
sample_n <- summary(Ba_muscle_cv)$df[2]
Ba_muscle_output <- cbind("muscle",Ba_muscle_coeff, Ba_muscle_loocv, sample_n)
colnames(Ba_muscle_output)[1]<-"tissue"

#Liver
Ba_liver<-lm(liver_Ba ~ Feces_Ba, data=metals_log)
Ba_liver_coeff <- as.data.frame(tidy(Ba_liver))
Ba_liver_cv <- train(liver_Ba ~ Feces_Ba, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Ba_liver_loocv<- as.data.frame(Ba_liver_cv$results)
sample_n <- summary(Ba_liver_cv)$df[2]
Ba_liver_output <- cbind("liver",Ba_liver_coeff, Ba_liver_loocv, sample_n)
colnames(Ba_liver_output)[1]<-"tissue"

#fat
Ba_fat<-lm(fat_Ba ~ Feces_Ba, data=metals_log)
Ba_fat_coeff <- as.data.frame(tidy(Ba_fat))
Ba_fat_cv <- train(fat_Ba ~ Feces_Ba, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Ba_fat_loocv<- as.data.frame(Ba_fat_cv$results)
sample_n <- summary(Ba_fat_cv)$df[2]
Ba_fat_output <- cbind("fat", Ba_fat_coeff, Ba_fat_loocv, sample_n)
colnames(Ba_fat_output)[1]<-"tissue"

#write all outputs to csv
Ba_metal_outputs <- rbind (Ba_liver_output, Ba_muscle_output,Ba_fat_output)
write.csv(Ba_metal_outputs, "feces/metals/Ba_metal_outputs.csv")

#make plots 
Ba_fat <- ggplot(data=metals_log, aes(Feces_Ba, fat_Ba))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Ba (ug/g)", y="Fat Ba")

muscle_Ba <- ggplot(data=metals_log, aes(Feces_Ba, muscle_Ba))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Ba (ug/g)", y="Muscle Ba")
muscle_Ba

liver_Ba <- ggplot(data=metals_log, aes(Feces_Ba, liver_Ba))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Ba (ug/g)", y="Liver Ba")
liver_Ba

Ba_figure=ggarrange(liver_Ba, muscle_Ba,Ba_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
Ba_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_Ba.tif", Ba_figure, width = 20, height = 10, dpi = 300)

#### Be ####
# Muscle
Be_muscle<-lm(muscle_Be ~ Feces_Be, data=metals_log)
summary(Be_muscle)
Be_muscle_coeff <- as.data.frame(tidy(Be_muscle))
Be_muscle_cv <- train(muscle_Be ~ Feces_Be, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Be_muscle_loocv<- as.data.frame(Be_muscle_cv$results)
sample_n <- summary(Be_muscle_cv)$df[2]
Be_muscle_output <- cbind("muscle",Be_muscle_coeff, Be_muscle_loocv, sample_n)
colnames(Be_muscle_output)[1]<-"tissue"

#Liver
Be_liver<-lm(liver_Be ~ Feces_Be, data=metals_log)
Be_liver_coeff <- as.data.frame(tidy(Be_liver))
Be_liver_cv <- train(liver_Be ~ Feces_Be, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Be_liver_loocv<- as.data.frame(Be_liver_cv$results)
sample_n <- summary( Be_liver_cv)$df[2]
Be_liver_output <- cbind("liver",Be_liver_coeff, Be_liver_loocv, sample_n)
colnames(Be_liver_output)[1]<-"tissue"

#fat
Be_fat<-lm(fat_Be ~ Feces_Be, data=metals_log)
Be_fat_coeff <- as.data.frame(tidy(Be_fat))
Be_fat_cv <- train(fat_Be ~ Feces_Be, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Be_fat_loocv<- as.data.frame(Be_fat_cv$results)
sample_n <- summary( Be_fat_cv)$df[2]
Be_fat_output <- cbind("fat", Be_fat_coeff, Be_fat_loocv, sample_n)
colnames(Be_fat_output)[1]<-"tissue"

#write all outputs to csv
Be_metal_outputs <- rbind (Be_liver_output, Be_muscle_output,Be_fat_output)
write.csv(Be_metal_outputs, "feces/metals/Be_metal_outputs.csv")

#make plots 
Be_fat <- ggplot(data=metals_log, aes(Feces_Be, fat_Be))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Be (ug/g)", y="Fat Be")

muscle_Be <- ggplot(data=metals_log, aes(Feces_Be, muscle_Be))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Be (ug/g)", y="Muscle Be")
muscle_Be

liver_Be <- ggplot(data=metals_log, aes(Feces_Be, liver_Be))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Be (ug/g)", y="Liver Be")
liver_Be

Be_figure=ggarrange(liver_Be, muscle_Be,Be_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
Be_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_Be.tif", Be_figure, width = 20, height = 10, dpi = 300)

# #### B ####
# # Muscle
# B_muscle<-lm(muscle_B ~ Feces_B, data=metals_log)
# B_muscle_coeff <- as.data.frame(tidy(B_muscle))
# B_muscle_cv <- train(muscle_B ~ Feces_B, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
# B_muscle_loocv<- as.data.frame(B_muscle_cv$results)
# sample_n <- summary( B_liver_cv)$df[2]
# B_muscle_output <- cbind("muscle",B_muscle_coeff, B_muscle_loocv, sample_n)
# colnames(B_muscle_output)[1]<-"tissue"
# 
# #Liver
# B_liver<-lm(liver_B ~ Feces_B, data=metals_log)
# B_liver_coeff <- as.data.frame(tidy(B_liver))
# B_liver_cv <- train(liver_B ~ Feces_B, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
# B_liver_loocv<- as.data.frame(B_liver_cv$results)
# sample_n <- summary( B_liver_cv)$df[2]
# B_liver_output <- cbind("liver",B_liver_coeff, B_liver_loocv, sample_n)
# colnames(B_liver_output)[1]<-"tissue"
# 
# #fat
# B_fat<-lm(fat_B ~ Feces_B, data=metals_log)
# B_fat_coeff <- as.data.frame(tidy(B_fat))
# B_fat_cv <- train(fat_B ~ Feces_B, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
# B_fat_loocv<- as.data.frame(B_fat_cv$results)
# sample_n <- summary( B_fat_cv)$df[2]
# B_fat_output <- cbind("fat", B_fat_coeff, B_fat_loocv, sample_n)
# colnames(B_fat_output)[1]<-"tissue"
# 
# #write all outputs to csv
# B_metal_outputs <- rbind (B_liver_output, B_muscle_output,B_fat_output)
# write.csv(B_metal_outputs, "feces/metals/B_metal_outputs.csv")
# 
# #make plots 
# B_fat <- ggplot(data=metals_log, aes(Feces_B, fat_B))+
#   geom_point()+
#   geom_smooth(method = "lm", se = TRUE, color="red")+
#   theme_bw()+
#   labs(x="Feces B (ug/g)", y="Fat B")
# 
# muscle_B <- ggplot(data=metals_log, aes(Feces_B, muscle_B))+
#   geom_point()+
#   geom_smooth(method = "lm", se = TRUE, color="red")+
#   theme_bw()+
#   labs(x="Feces B (ug/g)", y="Muscle B")
# muscle_B
# 
# liver_B <- ggplot(data=metals_log, aes(Feces_B, liver_B))+
#   geom_point()+
#   geom_smooth(method = "lm", se = TRUE, color="red")+
#   theme_bw()+
#   labs(x="Feces B (ug/g)", y="Liver B")
# liver_B
# 
# B_figure=ggarrange(liver_B, muscle_B,B_fat, 
#                     labels = c("A", "B", "C"),
#                     vjust = 1,
#                     hjust = -0.5,
#                     ncol = 3, nrow = 1,
#                     legend = "right",
#                     font.label = list(size = 16))
# B_figure
# #Plot figures with dpi=300
# save_plot("feces/metals/metals_B.tif", B_figure, width = 20, height = 10, dpi = 300)

#### Sb ####
# Muscle
Sb_muscle<-lm(muscle_Sb ~ Feces_Sb, data=metals_log)
summary(Sb_muscle)
Sb_muscle_coeff <- as.data.frame(tidy(Sb_muscle))
Sb_muscle_cv <- train(muscle_Sb ~ Feces_Sb, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Sb_muscle_loocv<- as.data.frame(Sb_muscle_cv$results)
sample_n <- summary( Sb_muscle_cv)$df[2]
Sb_muscle_output <- cbind("muscle",Sb_muscle_coeff, Sb_muscle_loocv, sample_n)
colnames(Sb_muscle_output)[1]<-"tissue"

#Liver
Sb_liver<-lm(liver_Sb ~ Feces_Sb, data=metals_log)
Sb_liver_coeff <- as.data.frame(tidy(Sb_liver))
Sb_liver_cv <- train(liver_Sb ~ Feces_Sb, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Sb_liver_loocv<- as.data.frame(Sb_liver_cv$results)
sample_n <- summary( Sb_liver_cv)$df[2]
Sb_liver_output <- cbind("liver",Sb_liver_coeff, Sb_liver_loocv, sample_n)
colnames(Sb_liver_output)[1]<-"tissue"

#fat
Sb_fat<-lm(fat_Sb ~ Feces_Sb, data=metals_log)
Sb_fat_coeff <- as.data.frame(tidy(Sb_fat))
Sb_fat_cv <- train(fat_Sb ~ Feces_Sb, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Sb_fat_loocv<- as.data.frame(Sb_fat_cv$results)
sample_n <- summary( Sb_fat_cv)$df[2]
Sb_fat_output <- cbind("fat", Sb_fat_coeff, Sb_fat_loocv, sample_n)
colnames(Sb_fat_output)[1]<-"tissue"

#write all outputs to csv
Sb_metal_outputs <- rbind (Sb_liver_output, Sb_muscle_output,Sb_fat_output)
write.csv(Sb_metal_outputs, "feces/metals/Sb_metal_outputs.csv")

#make plots 
Sb_fat <- ggplot(data=metals_log, aes(Feces_Sb, fat_Sb))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Sb (ug/g)", y="Fat Sb")

muscle_Sb <- ggplot(data=metals_log, aes(Feces_Sb, muscle_Sb))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Sb (ug/g)", y="Muscle Sb")
muscle_Sb

liver_Sb <- ggplot(data=metals_log, aes(Feces_Sb, liver_Sb))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Sb (ug/g)", y="Liver Sb")
liver_Sb

Sb_figure=ggarrange(liver_Sb, muscle_Sb,Sb_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
Sb_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_Sb.tif", Sb_figure, width = 20, height = 10, dpi = 300)

#### Cr ####
# Muscle
Cr_muscle<-lm(muscle_Cr ~ Feces_Cr, data=metals_log)
summary(Cr_muscle)
Cr_muscle_coeff <- as.data.frame(tidy(Cr_muscle))
Cr_muscle_cv <- train(muscle_Cr ~ Feces_Cr, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Cr_muscle_loocv<- as.data.frame(Cr_muscle_cv$results)
sample_n <- summary( Cr_muscle_cv)$df[2]
Cr_muscle_output <- cbind("muscle",Cr_muscle_coeff, Cr_muscle_loocv, sample_n)
colnames(Cr_muscle_output)[1]<-"tissue"

#Liver
Cr_liver<-lm(liver_Cr ~ Feces_Cr, data=metals_log)
Cr_liver_coeff <- as.data.frame(tidy(Cr_liver))
Cr_liver_cv <- train(liver_Cr ~ Feces_Cr, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Cr_liver_loocv<- as.data.frame(Cr_liver_cv$results)
sample_n <- summary( Cr_liver_cv)$df[2]
Cr_liver_output <- cbind("liver",Cr_liver_coeff, Cr_liver_loocv, sample_n)
colnames(Cr_liver_output)[1]<-"tissue"

#fat
Cr_fat<-lm(fat_Cr ~ Feces_Cr, data=metals_log)
Cr_fat_coeff <- as.data.frame(tidy(Cr_fat))
Cr_fat_cv <- train(fat_Cr ~ Feces_Cr, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Cr_fat_loocv<- as.data.frame(Cr_fat_cv$results)
sample_n <- summary( Cr_fat_cv)$df[2]
Cr_fat_output <- cbind("fat", Cr_fat_coeff, Cr_fat_loocv, sample_n)
colnames(Cr_fat_output)[1]<-"tissue"

#write all outputs to csv
Cr_metal_outputs <- rbind (Cr_liver_output, Cr_muscle_output,Cr_fat_output)
write.csv(Cr_metal_outputs, "feces/metals/Cr_metal_outputs.csv")

#make plots 
Cr_fat <- ggplot(data=metals_log, aes(Feces_Cr, fat_Cr))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Cr (ug/g)", y="Fat Cr")

muscle_Cr <- ggplot(data=metals_log, aes(Feces_Cr, muscle_Cr))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Cr (ug/g)", y="Muscle Cr")
muscle_Cr

liver_Cr <- ggplot(data=metals_log, aes(Feces_Cr, liver_Cr))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Cr (ug/g)", y="Liver Cr")
liver_Cr

Cr_figure=ggarrange(liver_Cr, muscle_Cr,Cr_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
Cr_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_Cr.tif", Cr_figure, width = 20, height = 10, dpi = 300)


#### Co ####
# Muscle
Co_muscle<-lm(muscle_Co ~ Feces_Co, data=metals_log)
summary(Co_muscle)
Co_muscle_coeff <- as.data.frame(tidy(Co_muscle))
Co_muscle_cv <- train(muscle_Co ~ Feces_Co, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Co_muscle_loocv<- as.data.frame(Co_muscle_cv$results)
sample_n <- summary( Co_muscle_cv)$df[2]
Co_muscle_output <- cbind("muscle",Co_muscle_coeff, Co_muscle_loocv, sample_n)
colnames(Co_muscle_output)[1]<-"tissue"

#Liver
Co_liver<-lm(liver_Co ~ Feces_Co, data=metals_log)
Co_liver_coeff <- as.data.frame(tidy(Co_liver))
Co_liver_cv <- train(liver_Co ~ Feces_Co, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Co_liver_loocv<- as.data.frame(Co_liver_cv$results)
sample_n <- summary( Co_liver_cv)$df[2]
Co_liver_output <- cbind("liver",Co_liver_coeff, Co_liver_loocv, sample_n)
colnames(Co_liver_output)[1]<-"tissue"

#fat
Co_fat<-lm(fat_Co ~ Feces_Co, data=metals_log)
Co_fat_coeff <- as.data.frame(tidy(Co_fat))
Co_fat_cv <- train(fat_Co ~ Feces_Co, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Co_fat_loocv<- as.data.frame(Co_fat_cv$results)
sample_n <- summary( Co_fat_cv)$df[2]
Co_fat_output <- cbind("fat", Co_fat_coeff, Co_fat_loocv, sample_n)
colnames(Co_fat_output)[1]<-"tissue"

#write all outputs to csv
Co_metal_outputs <- rbind (Co_liver_output, Co_muscle_output,Co_fat_output)
write.csv(Co_metal_outputs, "feces/metals/Co_metal_outputs.csv")

#make plots 
Co_fat <- ggplot(data=metals_log, aes(Feces_Co, fat_Co))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Co (ug/g)", y="Fat Co")

muscle_Co <- ggplot(data=metals_log, aes(Feces_Co, muscle_Co))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Co (ug/g)", y="Muscle Co")
muscle_Co

liver_Co <- ggplot(data=metals_log, aes(Feces_Co, liver_Co))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Co (ug/g)", y="Liver Co")
liver_Co

Co_figure=ggarrange(liver_Co, muscle_Co,Co_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
Co_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_Co.tif", Co_figure, width = 20, height = 10, dpi = 300)


#### Cu ####
# Muscle
Cu_muscle<-lm(muscle_Cu ~ Feces_Cu, data=metals_log)
summary(Cu_muscle)
Cu_muscle_coeff <- as.data.frame(tidy(Cu_muscle))
Cu_muscle_cv <- train(muscle_Cu ~ Feces_Cu, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Cu_muscle_loocv<- as.data.frame(Cu_muscle_cv$results)
sample_n <- summary( Cu_muscle_cv)$df[2]
Cu_muscle_output <- cbind("muscle",Cu_muscle_coeff, Cu_muscle_loocv, sample_n)
colnames(Cu_muscle_output)[1]<-"tissue"

#Liver
Cu_liver<-lm(liver_Cu ~ Feces_Cu, data=metals_log)
Cu_liver_coeff <- as.data.frame(tidy(Cu_liver))
Cu_liver_cv <- train(liver_Cu ~ Feces_Cu, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Cu_liver_loocv<- as.data.frame(Cu_liver_cv$results)
sample_n <- summary( Cu_liver_cv)$df[2]
Cu_liver_output <- cbind("liver",Cu_liver_coeff, Cu_liver_loocv, sample_n)
colnames(Cu_liver_output)[1]<-"tissue"

#fat
Cu_fat<-lm(fat_Cu ~ Feces_Cu, data=metals_log)
Cu_fat_coeff <- as.data.frame(tidy(Cu_fat))
Cu_fat_cv <- train(fat_Cu ~ Feces_Cu, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Cu_fat_loocv<- as.data.frame(Cu_fat_cv$results)
sample_n <- summary( Cu_fat_cv)$df[2]
Cu_fat_output <- cbind("fat", Cu_fat_coeff, Cu_fat_loocv, sample_n)
colnames(Cu_fat_output)[1]<-"tissue"

#write all outputs to csv
Cu_metal_outputs <- rbind (Cu_liver_output, Cu_muscle_output,Cu_fat_output)
write.csv(Cu_metal_outputs, "feces/metals/Cu_metal_outputs.csv")

#make plots 
Cu_fat <- ggplot(data=metals_log, aes(Feces_Cu, fat_Cu))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Cu (ug/g)", y="Fat Cu")

muscle_Cu <- ggplot(data=metals_log, aes(Feces_Cu, muscle_Cu))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Cu (ug/g)", y="Muscle Cu")
muscle_Cu

liver_Cu <- ggplot(data=metals_log, aes(Feces_Cu, liver_Cu))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Cu (ug/g)", y="Liver Cu")
liver_Cu

Cu_figure=ggarrange(liver_Cu, muscle_Cu,Cu_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
Cu_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_Cu.tif", Cu_figure, width = 20, height = 10, dpi = 300)


#### Mg ####
# Muscle
Mg_muscle<-lm(muscle_Mg ~ Feces_Mg, data=metals_log)
summary(Mg_muscle)
Mg_muscle_coeff <- as.data.frame(tidy(Mg_muscle))
Mg_muscle_cv <- train(muscle_Mg ~ Feces_Mg, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Mg_muscle_loocv<- as.data.frame(Mg_muscle_cv$results)
sample_n <- summary( Mg_muscle_cv)$df[2]
Mg_muscle_output <- cbind("muscle",Mg_muscle_coeff, Mg_muscle_loocv, sample_n)
colnames(Mg_muscle_output)[1]<-"tissue"

#Liver
Mg_liver<-lm(liver_Mg ~ Feces_Mg, data=metals_log)
Mg_liver_coeff <- as.data.frame(tidy(Mg_liver))
Mg_liver_cv <- train(liver_Mg ~ Feces_Mg, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Mg_liver_loocv<- as.data.frame(Mg_liver_cv$results)
sample_n <- summary( Mg_liver_cv)$df[2]
Mg_liver_output <- cbind("liver",Mg_liver_coeff, Mg_liver_loocv, sample_n)
colnames(Mg_liver_output)[1]<-"tissue"

#fat
Mg_fat<-lm(fat_Mg ~ Feces_Mg, data=metals_log)
Mg_fat_coeff <- as.data.frame(tidy(Mg_fat))
Mg_fat_cv <- train(fat_Mg ~ Feces_Mg, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Mg_fat_loocv<- as.data.frame(Mg_fat_cv$results)
sample_n <- summary( Mg_fat_cv)$df[2]
Mg_fat_output <- cbind("fat", Mg_fat_coeff, Mg_fat_loocv, sample_n)
colnames(Mg_fat_output)[1]<-"tissue"

#write all outputs to csv
Mg_metal_outputs <- rbind (Mg_liver_output, Mg_muscle_output,Mg_fat_output)
write.csv(Mg_metal_outputs, "feces/metals/Mg_metal_outputs.csv")

#make plots 
Mg_fat <- ggplot(data=metals_log, aes(Feces_Mg, fat_Mg))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Mg (ug/g)", y="Fat Mg")

muscle_Mg <- ggplot(data=metals_log, aes(Feces_Mg, muscle_Mg))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Mg (ug/g)", y="Muscle Mg")
muscle_Mg

liver_Mg <- ggplot(data=metals_log, aes(Feces_Mg, liver_Mg))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Mg (ug/g)", y="Liver Mg")
liver_Mg

Mg_figure=ggarrange(liver_Mg, muscle_Mg,Mg_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
Mg_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_Mg.tif", Mg_figure, width = 20, height = 10, dpi = 300)


#### Mn ####
# Muscle
Mn_muscle<-lm(muscle_Mn ~ Feces_Mn, data=metals_log)
summary(Mn_muscle)
Mn_muscle_coeff <- as.data.frame(tidy(Mn_muscle))
Mn_muscle_cv <- train(muscle_Mn ~ Feces_Mn, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Mn_muscle_loocv<- as.data.frame(Mn_muscle_cv$results)
sample_n <- summary( Mn_muscle_cv)$df[2]
Mn_muscle_output <- cbind("muscle",Mn_muscle_coeff, Mn_muscle_loocv, sample_n)
colnames(Mn_muscle_output)[1]<-"tissue"

#Liver
Mn_liver<-lm(liver_Mn ~ Feces_Mn, data=metals_log)
Mn_liver_coeff <- as.data.frame(tidy(Mn_liver))
Mn_liver_cv <- train(liver_Mn ~ Feces_Mn, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Mn_liver_loocv<- as.data.frame(Mn_liver_cv$results)
sample_n <- summary( Mn_liver_cv)$df[2]
Mn_liver_output <- cbind("liver",Mn_liver_coeff, Mn_liver_loocv, sample_n)
colnames(Mn_liver_output)[1]<-"tissue"

#fat
Mn_fat<-lm(fat_Mn ~ Feces_Mn, data=metals_log)
Mn_fat_coeff <- as.data.frame(tidy(Mn_fat))
Mn_fat_cv <- train(fat_Mn ~ Feces_Mn, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Mn_fat_loocv<- as.data.frame(Mn_fat_cv$results)
sample_n <- summary( Mn_fat_cv)$df[2]
Mn_fat_output <- cbind("fat", Mn_fat_coeff, Mn_fat_loocv, sample_n)
colnames(Mn_fat_output)[1]<-"tissue"

#write all outputs to csv
Mn_metal_outputs <- rbind (Mn_liver_output, Mn_muscle_output,Mn_fat_output)
write.csv(Mn_metal_outputs, "feces/metals/Mn_metal_outputs.csv")

#make plots 
Mn_fat <- ggplot(data=metals_log, aes(Feces_Mn, fat_Mn))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Mn (ug/g)", y="Fat Mn")

muscle_Mn <- ggplot(data=metals_log, aes(Feces_Mn, muscle_Mn))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Mn (ug/g)", y="Muscle Mn")
muscle_Mn

liver_Mn <- ggplot(data=metals_log, aes(Feces_Mn, liver_Mn))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Mn (ug/g)", y="Liver Mn")
liver_Mn

Mn_figure=ggarrange(liver_Mn, muscle_Mn,Mn_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
Mn_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_Mn.tif", Mn_figure, width = 20, height = 10, dpi = 300)


#### Mo ####
# Muscle
Mo_muscle<-lm(muscle_Mo ~ Feces_Mo, data=metals_log)
summary(Mo_muscle)
Mo_muscle_coeff <- as.data.frame(tidy(Mo_muscle))
Mo_muscle_cv <- train(muscle_Mo ~ Feces_Mo, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Mo_muscle_loocv<- as.data.frame(Mo_muscle_cv$results)
sample_n <- summary( Mo_muscle_cv)$df[2]
Mo_muscle_output <- cbind("muscle",Mo_muscle_coeff, Mo_muscle_loocv, sample_n)
colnames(Mo_muscle_output)[1]<-"tissue"

#Liver
Mo_liver<-lm(liver_Mo ~ Feces_Mo, data=metals_log)
Mo_liver_coeff <- as.data.frame(tidy(Mo_liver))
Mo_liver_cv <- train(liver_Mo ~ Feces_Mo, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Mo_liver_loocv<- as.data.frame(Mo_liver_cv$results)
sample_n <- summary( Mo_liver_cv)$df[2]
Mo_liver_output <- cbind("liver",Mo_liver_coeff, Mo_liver_loocv, sample_n)
colnames(Mo_liver_output)[1]<-"tissue"

#fat
Mo_fat<-lm(fat_Mo ~ Feces_Mo, data=metals_log)
Mo_fat_coeff <- as.data.frame(tidy(Mo_fat))
Mo_fat_cv <- train(fat_Mo ~ Feces_Mo, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Mo_fat_loocv<- as.data.frame(Mo_fat_cv$results)
sample_n <- summary( Mo_fat_cv)$df[2]
Mo_fat_output <- cbind("fat", Mo_fat_coeff, Mo_fat_loocv, sample_n)
colnames(Mo_fat_output)[1]<-"tissue"

#write all outputs to csv
Mo_metal_outputs <- rbind (Mo_liver_output, Mo_muscle_output,Mo_fat_output)
write.csv(Mo_metal_outputs, "feces/metals/Mo_metal_outputs.csv")

#make plots 
Mo_fat <- ggplot(data=metals_log, aes(Feces_Mo, fat_Mo))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Mo (ug/g)", y="Fat Mo")

muscle_Mo <- ggplot(data=metals_log, aes(Feces_Mo, muscle_Mo))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Mo (ug/g)", y="Muscle Mo")
muscle_Mo

liver_Mo <- ggplot(data=metals_log, aes(Feces_Mo, liver_Mo))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Mo (ug/g)", y="Liver Mo")
liver_Mo

Mo_figure=ggarrange(liver_Mo, muscle_Mo,Mo_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
Mo_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_Mo.tif", Mo_figure, width = 20, height = 10, dpi = 300)


#### Ni ####
# Muscle
Ni_muscle<-lm(muscle_Ni ~ Feces_Ni, data=metals_log)
summary(Ni_muscle)
Ni_muscle_coeff <- as.data.frame(tidy(Ni_muscle))
Ni_muscle_cv <- train(muscle_Ni ~ Feces_Ni, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Ni_muscle_loocv<- as.data.frame(Ni_muscle_cv$results)
sample_n <- summary( Ni_muscle_cv)$df[2]
Ni_muscle_output <- cbind("muscle",Ni_muscle_coeff, Ni_muscle_loocv, sample_n)
colnames(Ni_muscle_output)[1]<-"tissue"

#Liver
Ni_liver<-lm(liver_Ni ~ Feces_Ni, data=metals_log)
Ni_liver_coeff <- as.data.frame(tidy(Ni_liver))
Ni_liver_cv <- train(liver_Ni ~ Feces_Ni, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Ni_liver_loocv<- as.data.frame(Ni_liver_cv$results)
sample_n <- summary( Ni_liver_cv)$df[2]
Ni_liver_output <- cbind("liver",Ni_liver_coeff, Ni_liver_loocv, sample_n)
colnames(Ni_liver_output)[1]<-"tissue"

#fat
Ni_fat<-lm(fat_Ni ~ Feces_Ni, data=metals_log)
Ni_fat_coeff <- as.data.frame(tidy(Ni_fat))
Ni_fat_cv <- train(fat_Ni ~ Feces_Ni, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Ni_fat_loocv<- as.data.frame(Ni_fat_cv$results)
sample_n <- summary( Ni_fat_cv)$df[2]
Ni_fat_output <- cbind("fat", Ni_fat_coeff, Ni_fat_loocv, sample_n)
colnames(Ni_fat_output)[1]<-"tissue"

#write all outputs to csv
Ni_metal_outputs <- rbind (Ni_liver_output, Ni_muscle_output,Ni_fat_output)
write.csv(Ni_metal_outputs, "feces/metals/Ni_metal_outputs.csv")

#make plots 
Ni_fat <- ggplot(data=metals_log, aes(Feces_Ni, fat_Ni))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Ni (ug/g)", y="Fat Ni")

muscle_Ni <- ggplot(data=metals_log, aes(Feces_Ni, muscle_Ni))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Ni (ug/g)", y="Muscle Ni")
muscle_Ni

liver_Ni <- ggplot(data=metals_log, aes(Feces_Ni, liver_Ni))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Ni (ug/g)", y="Liver Ni")
liver_Ni

Ni_figure=ggarrange(liver_Ni, muscle_Ni,Ni_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
Ni_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_Ni.tif", Ni_figure, width = 20, height = 10, dpi = 300)


#### P ####
# Muscle
P_muscle<-lm(muscle_P ~ Feces_P, data=metals_log)
summary(P_muscle)
P_muscle_coeff <- as.data.frame(tidy(P_muscle))
P_muscle_cv <- train(muscle_P ~ Feces_P, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
P_muscle_loocv<- as.data.frame(P_muscle_cv$results)
sample_n <- summary( P_muscle_cv)$df[2]
P_muscle_output <- cbind("muscle",P_muscle_coeff, P_muscle_loocv, sample_n)
colnames(P_muscle_output)[1]<-"tissue"

#Liver
P_liver<-lm(liver_P ~ Feces_P, data=metals_log)
P_liver_coeff <- as.data.frame(tidy(P_liver))
P_liver_cv <- train(liver_P ~ Feces_P, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
P_liver_loocv<- as.data.frame(P_liver_cv$results)
sample_n <- summary( P_liver_cv)$df[2]
P_liver_output <- cbind("liver",P_liver_coeff, P_liver_loocv, sample_n)
colnames(P_liver_output)[1]<-"tissue"

#fat
P_fat<-lm(fat_P ~ Feces_P, data=metals_log)
P_fat_coeff <- as.data.frame(tidy(P_fat))
P_fat_cv <- train(fat_P ~ Feces_P, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
P_fat_loocv<- as.data.frame(P_fat_cv$results)
sample_n <- summary( P_fat_cv)$df[2]
P_fat_output <- cbind("fat", P_fat_coeff, P_fat_loocv, sample_n)
colnames(P_fat_output)[1]<-"tissue"

#write all outputs to csv
P_metal_outputs <- rbind (P_liver_output, P_muscle_output,P_fat_output)
write.csv(P_metal_outputs, "feces/metals/P_metal_outputs.csv")

#make plots 
P_fat <- ggplot(data=metals_log, aes(Feces_P, fat_P))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces P (ug/g)", y="Fat P")

muscle_P <- ggplot(data=metals_log, aes(Feces_P, muscle_P))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces P (ug/g)", y="Muscle P")
muscle_P

liver_P <- ggplot(data=metals_log, aes(Feces_P, liver_P))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces P (ug/g)", y="Liver P")
liver_P

P_figure=ggarrange(liver_P, muscle_P,P_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
P_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_P.tif", P_figure, width = 20, height = 10, dpi = 300)


#### Ag ####
# Muscle
Ag_muscle<-lm(muscle_Ag ~ Feces_Ag, data=metals_log)
summary(Ag_muscle)
Ag_muscle_coeff <- as.data.frame(tidy(Ag_muscle))
Ag_muscle_cv <- train(muscle_Ag ~ Feces_Ag, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Ag_muscle_loocv<- as.data.frame(Ag_muscle_cv$results)
sample_n <- summary( Ag_muscle_cv)$df[2]
Ag_muscle_output <- cbind("muscle",Ag_muscle_coeff, Ag_muscle_loocv, sample_n)
colnames(Ag_muscle_output)[1]<-"tissue"

#Liver
Ag_liver<-lm(liver_Ag ~ Feces_Ag, data=metals_log)
Ag_liver_coeff <- as.data.frame(tidy(Ag_liver))
Ag_liver_cv <- train(liver_Ag ~ Feces_Ag, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Ag_liver_loocv<- as.data.frame(Ag_liver_cv$results)
sample_n <- summary( Ag_liver_cv)$df[2]
Ag_liver_output <- cbind("liver",Ag_liver_coeff, Ag_liver_loocv, sample_n)
colnames(Ag_liver_output)[1]<-"tissue"

#fat
Ag_fat<-lm(fat_Ag ~ Feces_Ag, data=metals_log)
Ag_fat_coeff <- as.data.frame(tidy(Ag_fat))
Ag_fat_cv <- train(fat_Ag ~ Feces_Ag, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Ag_fat_loocv<- as.data.frame(Ag_fat_cv$results)
sample_n <- summary( Ag_fat_cv)$df[2]
Ag_fat_output <- cbind("fat", Ag_fat_coeff, Ag_fat_loocv, sample_n)
colnames(Ag_fat_output)[1]<-"tissue"

#write all outputs to csv
Ag_metal_outputs <- rbind (Ag_liver_output, Ag_muscle_output,Ag_fat_output)
write.csv(Ag_metal_outputs, "feces/metals/Ag_metal_outputs.csv")

#make plots 
Ag_fat <- ggplot(data=metals_log, aes(Feces_Ag, fat_Ag))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Ag (ug/g)", y="Fat Ag")

muscle_Ag <- ggplot(data=metals_log, aes(Feces_Ag, muscle_Ag))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Ag (ug/g)", y="Muscle Ag")
muscle_Ag

liver_Ag <- ggplot(data=metals_log, aes(Feces_Ag, liver_Ag))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Ag (ug/g)", y="Liver Ag")
liver_Ag

Ag_figure=ggarrange(liver_Ag, muscle_Ag,Ag_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
Ag_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_Ag.tif", Ag_figure, width = 20, height = 10, dpi = 300)


#### Na ####
# Muscle
Na_muscle<-lm(muscle_Na ~ Feces_Na, data=metals_log)
summary(Na_muscle)
Na_muscle_coeff <- as.data.frame(tidy(Na_muscle))
Na_muscle_cv <- train(muscle_Na ~ Feces_Na, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Na_muscle_loocv<- as.data.frame(Na_muscle_cv$results)
sample_n <- summary( Na_muscle_cv)$df[2]
Na_muscle_output <- cbind("muscle",Na_muscle_coeff, Na_muscle_loocv, sample_n)
colnames(Na_muscle_output)[1]<-"tissue"

#Liver
Na_liver<-lm(liver_Na ~ Feces_Na, data=metals_log)
Na_liver_coeff <- as.data.frame(tidy(Na_liver))
Na_liver_cv <- train(liver_Na ~ Feces_Na, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Na_liver_loocv<- as.data.frame(Na_liver_cv$results)
sample_n <- summary( Na_liver_cv)$df[2]
Na_liver_output <- cbind("liver",Na_liver_coeff, Na_liver_loocv, sample_n)
colnames(Na_liver_output)[1]<-"tissue"

#fat
Na_fat<-lm(fat_Na ~ Feces_Na, data=metals_log)
Na_fat_coeff <- as.data.frame(tidy(Na_fat))
Na_fat_cv <- train(fat_Na ~ Feces_Na, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Na_fat_loocv<- as.data.frame(Na_fat_cv$results)
sample_n <- summary( Na_fat_cv)$df[2]
Na_fat_output <- cbind("fat", Na_fat_coeff, Na_fat_loocv, sample_n)
colnames(Na_fat_output)[1]<-"tissue"

#write all outputs to csv
Na_metal_outputs <- rbind (Na_liver_output, Na_muscle_output,Na_fat_output)
write.csv(Na_metal_outputs, "feces/metals/Na_metal_outputs.csv")

#make plots 
Na_fat <- ggplot(data=metals_log, aes(Feces_Na, fat_Na))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Na (ug/g)", y="Fat Na")

muscle_Na <- ggplot(data=metals_log, aes(Feces_Na, muscle_Na))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Na (ug/g)", y="Muscle Na")
muscle_Na

liver_Na <- ggplot(data=metals_log, aes(Feces_Na, liver_Na))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Na (ug/g)", y="Liver Na")
liver_Na

Na_figure=ggarrange(liver_Na, muscle_Na,Na_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
Na_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_Na.tif", Na_figure, width = 20, height = 10, dpi = 300)


#### Sr ####
# Muscle
Sr_muscle<-lm(muscle_Sr ~ Feces_Sr, data=metals_log)
summary(Sr_muscle)
Sr_muscle_coeff <- as.data.frame(tidy(Sr_muscle))
Sr_muscle_cv <- train(muscle_Sr ~ Feces_Sr, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Sr_muscle_loocv<- as.data.frame(Sr_muscle_cv$results)
sample_n <- summary( Sr_muscle_cv)$df[2]
Sr_muscle_output <- cbind("muscle",Sr_muscle_coeff, Sr_muscle_loocv, sample_n)
colnames(Sr_muscle_output)[1]<-"tissue"

#Liver
Sr_liver<-lm(liver_Sr ~ Feces_Sr, data=metals_log)
Sr_liver_coeff <- as.data.frame(tidy(Sr_liver))
Sr_liver_cv <- train(liver_Sr ~ Feces_Sr, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Sr_liver_loocv<- as.data.frame(Sr_liver_cv$results)
sample_n <- summary( Sr_liver_cv)$df[2]
Sr_liver_output <- cbind("liver",Sr_liver_coeff, Sr_liver_loocv, sample_n)
colnames(Sr_liver_output)[1]<-"tissue"

#fat
Sr_fat<-lm(fat_Sr ~ Feces_Sr, data=metals_log)
Sr_fat_coeff <- as.data.frame(tidy(Sr_fat))
Sr_fat_cv <- train(fat_Sr ~ Feces_Sr, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Sr_fat_loocv<- as.data.frame(Sr_fat_cv$results)
sample_n <- summary( Sr_fat_cv)$df[2]
Sr_fat_output <- cbind("fat", Sr_fat_coeff, Sr_fat_loocv, sample_n)
colnames(Sr_fat_output)[1]<-"tissue"

#write all outputs to csv
Sr_metal_outputs <- rbind (Sr_liver_output, Sr_muscle_output,Sr_fat_output)
write.csv(Sr_metal_outputs, "feces/metals/Sr_metal_outputs.csv")

#make plots 
Sr_fat <- ggplot(data=metals_log, aes(Feces_Sr, fat_Sr))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Sr (ug/g)", y="Fat Sr")

muscle_Sr <- ggplot(data=metals_log, aes(Feces_Sr, muscle_Sr))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Sr (ug/g)", y="Muscle Sr")
muscle_Sr

liver_Sr <- ggplot(data=metals_log, aes(Feces_Sr, liver_Sr))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Sr (ug/g)", y="Liver Sr")
liver_Sr

Sr_figure=ggarrange(liver_Sr, muscle_Sr,Sr_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
Sr_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_Sr.tif", Sr_figure, width = 20, height = 10, dpi = 300)


#### S ####
# Muscle
S_muscle<-lm(muscle_S ~ Feces_S, data=metals_log)
summary(S_muscle)
S_muscle_coeff <- as.data.frame(tidy(S_muscle))
S_muscle_cv <- train(muscle_S ~ Feces_S, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
S_muscle_loocv<- as.data.frame(S_muscle_cv$results)
sample_n <- summary( S_muscle_cv)$df[2]
S_muscle_output <- cbind("muscle",S_muscle_coeff, S_muscle_loocv, sample_n)
colnames(S_muscle_output)[1]<-"tissue"

#Liver
S_liver<-lm(liver_S ~ Feces_S, data=metals_log)
S_liver_coeff <- as.data.frame(tidy(S_liver))
S_liver_cv <- train(liver_S ~ Feces_S, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
S_liver_loocv<- as.data.frame(S_liver_cv$results)
sample_n <- summary( S_liver_cv)$df[2]
S_liver_output <- cbind("liver",S_liver_coeff, S_liver_loocv, sample_n)
colnames(S_liver_output)[1]<-"tissue"

#fat
S_fat<-lm(fat_S ~ Feces_S, data=metals_log)
S_fat_coeff <- as.data.frame(tidy(S_fat))
S_fat_cv <- train(fat_S ~ Feces_S, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
S_fat_loocv<- as.data.frame(S_fat_cv$results)
sample_n <- summary( S_fat_cv)$df[2]
S_fat_output <- cbind("fat", S_fat_coeff, S_fat_loocv, sample_n)
colnames(S_fat_output)[1]<-"tissue"

#write all outputs to csv
S_metal_outputs <- rbind (S_liver_output, S_muscle_output,S_fat_output)
write.csv(S_metal_outputs, "feces/metals/S_metal_outputs.csv")

#make plots 
S_fat <- ggplot(data=metals_log, aes(Feces_S, fat_S))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces S (ug/g)", y="Fat S")

muscle_S <- ggplot(data=metals_log, aes(Feces_S, muscle_S))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces S (ug/g)", y="Muscle S")
muscle_S

liver_S <- ggplot(data=metals_log, aes(Feces_S, liver_S))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces S (ug/g)", y="Liver S")
liver_S

S_figure=ggarrange(liver_S, muscle_S,S_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
S_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_S.tif", S_figure, width = 20, height = 10, dpi = 300)


#### Sn ####
# Muscle
Sn_muscle<-lm(muscle_Sn ~ Feces_Sn, data=metals_log)
summary(Sn_muscle)
Sn_muscle_coeff <- as.data.frame(tidy(Sn_muscle))
Sn_muscle_cv <- train(muscle_Sn ~ Feces_Sn, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Sn_muscle_loocv<- as.data.frame(Sn_muscle_cv$results)
sample_n <- summary( Sn_muscle_cv)$df[2]
Sn_muscle_output <- cbind("muscle",Sn_muscle_coeff, Sn_muscle_loocv, sample_n)
colnames(Sn_muscle_output)[1]<-"tissue"

#Liver
Sn_liver<-lm(liver_Sn ~ Feces_Sn, data=metals_log)
Sn_liver_coeff <- as.data.frame(tidy(Sn_liver))
Sn_liver_cv <- train(liver_Sn ~ Feces_Sn, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Sn_liver_loocv<- as.data.frame(Sn_liver_cv$results)
sample_n <- summary( Sn_liver_cv)$df[2]
Sn_liver_output <- cbind("liver",Sn_liver_coeff, Sn_liver_loocv, sample_n)
colnames(Sn_liver_output)[1]<-"tissue"

#fat
Sn_fat<-lm(fat_Sn ~ Feces_Sn, data=metals_log)
Sn_fat_coeff <- as.data.frame(tidy(Sn_fat))
Sn_fat_cv <- train(fat_Sn ~ Feces_Sn, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Sn_fat_loocv<- as.data.frame(Sn_fat_cv$results)
sample_n <- summary( Sn_fat_cv)$df[2]
Sn_fat_output <- cbind("fat", Sn_fat_coeff, Sn_fat_loocv, sample_n)
colnames(Sn_fat_output)[1]<-"tissue"

#write all outputs to csv
Sn_metal_outputs <- rbind (Sn_liver_output, Sn_muscle_output,Sn_fat_output)
write.csv(Sn_metal_outputs, "feces/metals/Sn_metal_outputs.csv")

#make plots 
Sn_fat <- ggplot(data=metals_log, aes(Feces_Sn, fat_Sn))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Sn (ug/g)", y="Fat Sn")

muscle_Sn <- ggplot(data=metals_log, aes(Feces_Sn, muscle_Sn))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Sn (ug/g)", y="Muscle Sn")
muscle_Sn

liver_Sn <- ggplot(data=metals_log, aes(Feces_Sn, liver_Sn))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Sn (ug/g)", y="Liver Sn")
liver_Sn

Sn_figure=ggarrange(liver_Sn, muscle_Sn,Sn_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
Sn_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_Sn.tif", Sn_figure, width = 20, height = 10, dpi = 300)


#### U ####
# Muscle
U_muscle<-lm(muscle_U ~ Feces_U, data=metals_log)
summary(U_muscle)
U_muscle_coeff <- as.data.frame(tidy(U_muscle))
U_muscle_cv <- train(muscle_U ~ Feces_U, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
U_muscle_loocv<- as.data.frame(U_muscle_cv$results)
sample_n <- summary( U_muscle_cv)$df[2]
U_muscle_output <- cbind("muscle",U_muscle_coeff, U_muscle_loocv, sample_n)
colnames(U_muscle_output)[1]<-"tissue"

#Liver
U_liver<-lm(liver_U ~ Feces_U, data=metals_log)
U_liver_coeff <- as.data.frame(tidy(U_liver))
U_liver_cv <- train(liver_U ~ Feces_U, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
U_liver_loocv<- as.data.frame(U_liver_cv$results)
sample_n <- summary( U_liver_cv)$df[2]
U_liver_output <- cbind("liver",U_liver_coeff, U_liver_loocv, sample_n)
colnames(U_liver_output)[1]<-"tissue"

#fat
U_fat<-lm(fat_U ~ Feces_U, data=metals_log)
U_fat_coeff <- as.data.frame(tidy(U_fat))
U_fat_cv <- train(fat_U ~ Feces_U, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
U_fat_loocv<- as.data.frame(U_fat_cv$results)
sample_n <- summary( U_fat_cv)$df[2]
U_fat_output <- cbind("fat", U_fat_coeff, U_fat_loocv, sample_n)
colnames(U_fat_output)[1]<-"tissue"

#write all outputs to csv
U_metal_outputs <- rbind (U_liver_output, U_muscle_output,U_fat_output)
write.csv(U_metal_outputs, "feces/metals/U_metal_outputs.csv")

#make plots 
U_fat <- ggplot(data=metals_log, aes(Feces_U, fat_U))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces U (ug/g)", y="Fat U")

muscle_U <- ggplot(data=metals_log, aes(Feces_U, muscle_U))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces U (ug/g)", y="Muscle U")
muscle_U

liver_U <- ggplot(data=metals_log, aes(Feces_U, liver_U))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces U (ug/g)", y="Liver U")
liver_U

U_figure=ggarrange(liver_U, muscle_U,U_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
U_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_U.tif", U_figure, width = 20, height = 10, dpi = 300)


#### Ti ####
# Muscle
Ti_muscle<-lm(muscle_Ti ~ Feces_Ti, data=metals_log)
summary(Ti_muscle)
Ti_muscle_coeff <- as.data.frame(tidy(Ti_muscle))
Ti_muscle_cv <- train(muscle_Ti ~ Feces_Ti, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Ti_muscle_loocv<- as.data.frame(Ti_muscle_cv$results)
sample_n <- summary( Ti_muscle_cv)$df[2]
Ti_muscle_output <- cbind("muscle",Ti_muscle_coeff, Ti_muscle_loocv, sample_n)
colnames(Ti_muscle_output)[1]<-"tissue"

#Liver
Ti_liver<-lm(liver_Ti ~ Feces_Ti, data=metals_log)
Ti_liver_coeff <- as.data.frame(tidy(Ti_liver))
Ti_liver_cv <- train(liver_Ti ~ Feces_Ti, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Ti_liver_loocv<- as.data.frame(Ti_liver_cv$results)
sample_n <- summary( Ti_liver_cv)$df[2]
Ti_liver_output <- cbind("liver",Ti_liver_coeff, Ti_liver_loocv, sample_n)
colnames(Ti_liver_output)[1]<-"tissue"

#fat
Ti_fat<-lm(fat_Ti ~ Feces_Ti, data=metals_log)
Ti_fat_coeff <- as.data.frame(tidy(Ti_fat))
Ti_fat_cv <- train(fat_Ti ~ Feces_Ti, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Ti_fat_loocv<- as.data.frame(Ti_fat_cv$results)
sample_n <- summary( Ti_fat_cv)$df[2]
Ti_fat_output <- cbind("fat", Ti_fat_coeff, Ti_fat_loocv, sample_n)
colnames(Ti_fat_output)[1]<-"tissue"

#write all outputs to csv
Ti_metal_outputs <- rbind (Ti_liver_output, Ti_muscle_output,Ti_fat_output)
write.csv(Ti_metal_outputs, "feces/metals/Ti_metal_outputs.csv")

#make plots 
Ti_fat <- ggplot(data=metals_log, aes(Feces_Ti, fat_Ti))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Ti (ug/g)", y="Fat Ti")

muscle_Ti <- ggplot(data=metals_log, aes(Feces_Ti, muscle_Ti))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Ti (ug/g)", y="Muscle Ti")
muscle_Ti

liver_Ti <- ggplot(data=metals_log, aes(Feces_Ti, liver_Ti))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Ti (ug/g)", y="Liver Ti")
liver_Ti

Ti_figure=ggarrange(liver_Ti, muscle_Ti,Ti_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
Ti_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_Ti.tif", Ti_figure, width = 20, height = 10, dpi = 300)

#### TI ####
# Muscle
TI_muscle<-lm(muscle_TI ~ Feces_TI, data=metals_log)
summary(TI_muscle)
TI_muscle_coeff <- as.data.frame(tidy(TI_muscle))
TI_muscle_cv <- train(muscle_TI ~ Feces_TI, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
TI_muscle_loocv<- as.data.frame(TI_muscle_cv$results)
sample_n <- summary( TI_muscle_cv)$df[2]
TI_muscle_output <- cbind("muscle",TI_muscle_coeff, TI_muscle_loocv, sample_n)
colnames(TI_muscle_output)[1]<-"tissue"

#Liver
TI_liver<-lm(liver_TI ~ Feces_TI, data=metals_log)
TI_liver_coeff <- as.data.frame(tidy(TI_liver))
TI_liver_cv <- train(liver_TI ~ Feces_TI, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
TI_liver_loocv<- as.data.frame(TI_liver_cv$results)
sample_n <- summary( TI_liver_cv)$df[2]
TI_liver_output <- cbind("liver",TI_liver_coeff, TI_liver_loocv, sample_n)
colnames(TI_liver_output)[1]<-"tissue"

#fat
TI_fat<-lm(fat_TI ~ Feces_TI, data=metals_log)
TI_fat_coeff <- as.data.frame(tidy(TI_fat))
TI_fat_cv <- train(fat_TI ~ Feces_TI, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
TI_fat_loocv<- as.data.frame(TI_fat_cv$results)
sample_n <- summary( TI_fat_cv)$df[2]
TI_fat_output <- cbind("fat", TI_fat_coeff, TI_fat_loocv, sample_n)
colnames(TI_fat_output)[1]<-"tissue"

#write all outputs to csv
TI_metal_outputs <- rbind (TI_liver_output, TI_muscle_output,TI_fat_output)
write.csv(TI_metal_outputs, "feces/metals/Tl_metal_outputs.csv")

#make plots 
TI_fat <- ggplot(data=metals_log, aes(Feces_TI, fat_TI))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces TI (ug/g)", y="Fat TI")

muscle_TI <- ggplot(data=metals_log, aes(Feces_TI, muscle_TI))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces TI (ug/g)", y="Muscle TI")
muscle_TI

liver_TI <- ggplot(data=metals_log, aes(Feces_TI, liver_TI))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces TI (ug/g)", y="Liver TI")
liver_TI

TI_figure=ggarrange(liver_TI, muscle_TI,TI_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
TI_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_TI.tif", TI_figure, width = 20, height = 10, dpi = 300)

#### Ca ####
# Muscle
Ca_muscle<-lm(muscle_Ca ~ Feces_Ca, data=metals_log)
summary(Ca_muscle)
Ca_muscle_coeff <- as.data.frame(tidy(Ca_muscle))
Ca_muscle_cv <- train(muscle_Ca ~ Feces_Ca, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Ca_muscle_loocv<- as.data.frame(Ca_muscle_cv$results)
sample_n <- summary( Ca_muscle_cv)$df[2]
Ca_muscle_output <- cbind("muscle",Ca_muscle_coeff, Ca_muscle_loocv, sample_n)
colnames(Ca_muscle_output)[1]<-"tissue"

#Liver
Ca_liver<-lm(liver_Ca ~ Feces_Ca, data=metals_log)
Ca_liver_coeff <- as.data.frame(tidy(Ca_liver))
Ca_liver_cv <- train(liver_Ca ~ Feces_Ca, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Ca_liver_loocv<- as.data.frame(Ca_liver_cv$results)
sample_n <- summary( Ca_liver_cv)$df[2]
Ca_liver_output <- cbind("liver",Ca_liver_coeff, Ca_liver_loocv, sample_n)
colnames(Ca_liver_output)[1]<-"tissue"

#fat
Ca_fat<-lm(fat_Ca ~ Feces_Ca, data=metals_log)
Ca_fat_coeff <- as.data.frame(tidy(Ca_fat))
Ca_fat_cv <- train(fat_Ca ~ Feces_Ca, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
Ca_fat_loocv<- as.data.frame(Ca_fat_cv$results)
sample_n <- summary( Ca_fat_cv)$df[2]
Ca_fat_output <- cbind("fat", Ca_fat_coeff, Ca_fat_loocv, sample_n)
colnames(Ca_fat_output)[1]<-"tissue"

#write all outputs to csv
Ca_metal_outputs <- rbind (Ca_liver_output, Ca_muscle_output,Ca_fat_output)
write.csv(Ca_metal_outputs, "feces/metals/Ca_metal_outputs.csv")

#make plots 
Ca_fat <- ggplot(data=metals_log, aes(Feces_Ca, fat_Ca))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Ca (ug/g)", y="Fat Ca")

muscle_Ca <- ggplot(data=metals_log, aes(Feces_Ca, muscle_Ca))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Ca (ug/g)", y="Muscle Ca")
muscle_Ca

liver_Ca <- ggplot(data=metals_log, aes(Feces_Ca, liver_Ca))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color="red")+
  theme_bw()+
  labs(x="Feces Ca (ug/g)", y="Liver Ca")
liver_Ca

Ca_figure=ggarrange(liver_Ca, muscle_Ca,Ca_fat, 
                    labels = c("A", "B", "C"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 1,
                    legend = "right",
                    font.label = list(size = 16))
Ca_figure
#Plot figures with dpi=300
save_plot("feces/metals/metals_Ca.tif", Ca_figure, width = 20, height = 10, dpi = 300)


#### Add Together ####
data_all <- list.files(path = "feces/metals/",    
                       pattern = "*.csv", full.names = TRUE) %>% 
  lapply(read_csv) %>%                                           
  bind_rows %>%
  as.data.frame()
data_all
write.csv(data_all, "feces/metals/metals_data_all_loocv.csv", row.names = FALSE)

# Plot 
r2_plots <- ggplot(data = data_all, aes(x=Rsquared))+
  geom_histogram(bins = 10, aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")+ 
  theme_bw()+
  facet_grid(vars(tissue))+
  labs(x = "R-squared", y = "Density")
r2_plots




#### Composite Figure ####

#Plot figures with dpi=300
metals_composite= ggarrange(liver_As, muscle_As,As_fat, 
                            liver_Cd, muscle_Cd,Cd_fat, 
                            liver_Hg, muscle_Hg,Hg_fat, 
                            liver_MeHg, muscle_MeHg,MeHg_fat,  
                    labels = c("A", "B", "C",
                               "D","E","F",
                               "G","H","I",
                               "J","K","L"),
                    vjust = 1,
                    hjust = -0.5,
                    ncol = 3, nrow = 4,
                    legend = "right",
                    font.label = list(size = 14))
metals_composite

ggsave("feces/metals/metals_composite.jpg", metals_composite, width = 10, height = 10, dpi = 300)

#############################################################################################

#### Feces Relationships ####

#### Regression Model ####
ctrl <- trainControl(method = "boot", number = 1000)

#### Hg ####
# Muscle
MeHg_Hg_muscle<-lm(Feces_Hg ~ muscle_MeHg, data=metals_log)
MeHg_Hg_muscle_coeff <- as.data.frame(tidy(MeHg_Hg_muscle))
MeHg_Hg_muscle_cv <- train(Feces_Hg ~ muscle_MeHg, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
MeHg_Hg_muscle_loocv<- as.data.frame(MeHg_Hg_muscle_cv$results)
MeHg_sample_n <- summary(MeHg_Hg_muscle_cv)$df[2]
MeHg_Hg_muscle_output <- cbind("muscle",MeHg_Hg_muscle_coeff , MeHg_Hg_muscle_loocv, MeHg_sample_n)
colnames(MeHg_Hg_muscle_output)[1]<-"tissue"

# liver
MeHg_Hg_liver<-lm(Feces_Hg ~ liver_MeHg, data=metals_log)
MeHg_Hg_liver_coeff <- as.data.frame(tidy(MeHg_Hg_liver))
MeHg_Hg_liver_cv <- train(Feces_Hg ~ liver_MeHg, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
MeHg_Hg_liver_loocv<- as.data.frame(MeHg_Hg_liver_cv$results)
MeHg_sample_n <- summary(MeHg_Hg_liver_cv)$df[2]
MeHg_Hg_liver_output <- cbind("liver",MeHg_Hg_liver_coeff , MeHg_Hg_liver_loocv, MeHg_sample_n)
colnames(MeHg_Hg_liver_output)[1]<-"tissue"

# fat
MeHg_Hg_fat<-lm(Feces_Hg ~ fat_MeHg, data=metals_log)
MeHg_Hg_fat_coeff <- as.data.frame(tidy(MeHg_Hg_fat))
MeHg_Hg_fat_cv <- train(Feces_Hg ~ fat_MeHg, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
MeHg_Hg_fat_loocv<- as.data.frame(MeHg_Hg_fat_cv$results)
MeHg_sample_n <- summary(MeHg_Hg_fat_cv)$df[2]
MeHg_Hg_fat_output <- cbind("fat",MeHg_Hg_fat_coeff , MeHg_Hg_fat_loocv, MeHg_sample_n)
colnames(MeHg_Hg_fat_output)[1]<-"tissue"

# feces
MeHg_Hg_feces<-lm(Feces_Hg ~ Feces_MeHg, data=metals_log)
MeHg_Hg_feces_coeff <- as.data.frame(tidy(MeHg_Hg_feces))
MeHg_Hg_feces_cv <- train(Feces_Hg ~ Feces_MeHg, data = metals_log, method = "lm", trControl = ctrl, na.action=na.omit)
MeHg_Hg_feces_loocv<- as.data.frame(MeHg_Hg_feces_cv$results)
MeHg_sample_n <- summary(MeHg_Hg_feces_cv)$df[2]
MeHg_Hg_feces_output <- cbind("feces",MeHg_Hg_feces_coeff , MeHg_Hg_feces_loocv, MeHg_sample_n)
colnames(MeHg_Hg_feces_output)[1]<-"tissue"

mehg_hg_df <- rbind(MeHg_Hg_feces_output,MeHg_Hg_fat_output, MeHg_Hg_liver_output, MeHg_Hg_muscle_output )
write.csv(mehg_hg_df, "mehg_hg_df.csv")

