##########################################################################
# Summary statisticals for polar bear feces
# Written By: Kristin Eccles
# Date: September 22nd, 2022
# Updated: May 10th, 2024
# Note: Figure 3 is made from this script
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
library(vtable)
library(reshape2)
library(tidyverse)
library(broom)
library(ggpmisc)

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
mysheets <- read_excel_allsheets("orig_data/reg_tables.xlsx")

#### Prepare Data ####
Metals_liver <- mysheets$`Metals Liver`
liver_melt <- melt(Metals_liver, id="id")
liver_melt$tissue <-"Liver"

Metals_fat <- mysheets$`Metals Fat`
fat_melt <- melt(Metals_fat, id="id")
fat_melt$tissue <-"Fat"

Metals_feaces <- as.data.frame(mysheets$`Metals Feaces`)
feces_melt <- melt(Metals_feaces, id="id")
feces_melt$tissue <-"Feces"

Metals_muscle <- mysheets$`Metals Muscle`
muscle_melt <- melt(Metals_muscle, id="id")
muscle_melt$tissue <-"Muscle"

metal_df_stack <- rbind(liver_melt, muscle_melt,fat_melt, feces_melt)
metal_df_stack$value <-as.numeric(metal_df_stack$value)
metal_df_stack$log10value <- log10(metal_df_stack$value)

# read in sample info 
info_df <- read.csv("orig_data/polarbear_summary.csv", na.strings=c("","NA"))
metal_df_stack <- left_join(metal_df_stack, info_df, by = c("id" = "Identification"), keep=FALSE)

#### EDA ####
# summary table 
summary_stat <-metal_df_stack %>%                              
  group_by(tissue,variable) %>% 
  summarize(
            min = min(value, na.rm = TRUE),
            median = median(value,na.rm = TRUE),
            mean = mean(value,na.rm = TRUE),
            sd = sd(value,na.rm = TRUE),
            max = max(value,na.rm = TRUE))%>%
  as.data.frame()
write.csv(summary_stat, "summary_stat.csv")

# corrplot

#### T-test ####
# test if there is a difference between sex and age

age_test <- metal_df_stack %>%
  group_by(tissue, variable) %>%
  summarise(dif = t.test(log10value ~ Age, paired = FALSE, var.equal=FALSE, na.action = "na.omit")$statistic) %>%
  ungroup()%>%
  as.data.frame

age_p <- metal_df_stack %>%
  group_by(tissue, variable) %>%
  summarise(p.value = t.test(log10(value) ~ Age, paired = FALSE, var.equal=FALSE)$p.value,
            p.adj = p.adjust(p.value, method = "bonferroni", n = 25)) %>%
  ungroup()%>%
  as.data.frame

#age_test_sig <- subset(cbind(age_test, age_p), p.value < 0.05/length(age_p))

write.csv(cbind(age_test, age_p), "age_test.csv")

# Sex
sex_test <- metal_df_stack %>%
  group_by(tissue, variable) %>%
  summarise(dif = t.test(log10value ~ Sex, paired = FALSE, var.equal=FALSE)$statistic) %>%
  ungroup()%>%
  as.data.frame

sex_p <- metal_df_stack %>%
  group_by(tissue, variable) %>%
  summarise(p.value = t.test(log10(value) ~ Sex, paired = FALSE, var.equal=FALSE)$p.value,
            p.adj = p.adjust(p.value, method = "bonferroni", n = 25)) %>%
  ungroup()%>%
  as.data.frame

sex_test_sig <- subset(cbind(sex_test, sex_p), p.value < 0.05/25)

write.csv(cbind(sex_test, sex_p), "Sex_test.csv")

#### Difference Plots ####
hg_subset <- subset(metal_df_stack, variable == "Hg")
hg_subset2 <- hg_subset[!is.na(hg_subset$Age),]

A_boxplot= ggplot(hg_subset2, aes(x=tissue, y=log10value, fill=Age)) + 
  stat_boxplot(geom ='errorbar')+
  geom_boxplot()+
  scale_fill_brewer(palette = "YlOrRd")+
  theme_minimal(base_size = 15)+
  xlab("Tissue")+
  ylab("Log10 Hg Concentration (mg/kg)")+
  stat_compare_means(method="t.test", size=4, label.y = 2.5)+
  theme()
A_boxplot

B_boxplot= ggplot(hg_subset, aes(x=tissue, y=log10value, fill=Sex)) + 
  stat_boxplot(geom ='errorbar')+
  geom_boxplot()+
  scale_fill_brewer(palette = "YlOrRd")+
  theme_minimal(base_size = 15)+
  xlab("Tissue")+
  ylab("Log10 Hg Concentration (mg/kg)")+
  stat_compare_means(method="t.test", size=4, label.y = 2.5)+
  theme()
B_boxplot

#Compile
figure1=ggarrange(A_boxplot, B_boxplot,
                  labels = c("A", "B"),
                  vjust = 1,
                  hjust = -0.5,
                  ncol = 1, nrow = 2,
                  common.legend = FALSE,
                  legend = "right",
                  font.label = list(size = 16))
figure1
#Plot figures with dpi=300
cowplot::save_plot("Boxplot_Hg.jpg", figure1, base_width = 10, base_height = 10, dpi = 300)


#### feces Hg ratio ####
#unstack the data
mehg_subset <- subset(metal_df_stack, variable == "MeHg")

mehg_subset <- subset(metal_df_stack, variable == "MeHg")
mehg_subset2 <- mehg_subset[!is.na(mehg_subset$Age),]

C_boxplot= ggplot(mehg_subset2, aes(x=tissue, y=log10value, fill=Age)) + 
  stat_boxplot(geom ='errorbar')+
  geom_boxplot()+
  scale_fill_brewer(palette = "YlOrRd")+
  theme_minimal(base_size = 15)+
  xlab("Tissue")+
  ylab("Log10 MeHg Concentration (mg/kg)")+
  stat_compare_means(method="t.test", size=4, label.y = 2.5)+
  theme()
C_boxplot

D_boxplot= ggplot(mehg_subset, aes(x=tissue, y=log10value, fill=Sex)) + 
  stat_boxplot(geom ='errorbar')+
  geom_boxplot()+
  scale_fill_brewer(palette = "YlOrRd")+
  theme_minimal(base_size = 15)+
  xlab("Tissue")+
  ylab("Log10 MeHg Concentration (mg/kg)")+
  stat_compare_means(method="t.test", size=4, label.y = 2.5)+
  theme()
D_boxplot

#Compile
figure1=ggarrange(A_boxplot, B_boxplot,C_boxplot, D_boxplot,
                  labels = c("A", "B", "C", "D"),
                  vjust = 1,
                  hjust = -0.5,
                  ncol = 2, nrow = 2,
                  common.legend = FALSE,
                  legend = "right",
                  font.label = list(size = 16))
figure1
#Plot figures with dpi=300
cowplot::save_plot("Boxplot_Hg.jpg", figure1, base_width = 15, base_height = 10, dpi = 300)


unstack_hg <- hg_subset %>%
  pivot_wider(names_from = tissue, values_from = log10value, id_cols = id)
colnames(unstack_hg)<-paste(colnames(unstack_hg),"hg",sep="_")

unstack_mehg <- mehg_subset %>%
  pivot_wider(names_from = tissue, values_from = log10value, id_cols = id)
colnames(unstack_mehg)<-paste(colnames(unstack_mehg),"mehg",sep="_")

unstack_all_hg <- left_join(unstack_hg, unstack_mehg, by = c("id_hg" = "id_mehg"), keep = FALSE)
unstack_all_hg <- left_join(unstack_all_hg, info_df, by = c("id_hg" = "Identification"), keep = FALSE)
unstack_all_hg <- unstack_all_hg[!is.na(unstack_all_hg$Age),]

# feces: test for interaction ####
lm_feces_ratio=lm(data=unstack_all_hg, Feces_mehg~ Feces_hg  * Sex * Age)
anova(lm_feces_ratio)
summary(lm_feces_ratio)
AICc(lm_feces_ratio) #10.40307

lm_feces_ratio2=lm(data=unstack_all_hg, feces_mehg~ feces_hg  + Sex + Age)
anova(lm_feces_ratio2)
summary(lm_feces_ratio2)
AICc(lm_feces_ratio2) #6.220311


#### Feces Relationships ####
# Plot 
p1 <- ggplot(data=unstack_all_hg, aes(x=Feces_hg, y=Liver_mehg, colour=Age, shape=Sex))+
  geom_point(size=5)+
  scale_colour_brewer(palette = "YlOrRd")+
  theme_minimal(base_size = 15)+
  geom_smooth(data=unstack_all_hg, aes(x=Feces_hg, y=Liver_mehg, shape = NULL, color= NULL), 
              method='lm', formula= y~x, color = "black")+
  #stat_cor(data=unstack_all_hg, aes(x=Feces_hg, y=Liver_mehg , shape = NULL, color= NULL, 
  #                                  label = after_stat(rr.label)), color = "black", geom = "label")+
  xlab("Log10 Feces Hg (ug/g)")+
  ylab("Log10 Liver MeHg (ug/g)")
p1

p2 <- ggplot(data=unstack_all_hg, aes(x=Feces_hg, y=Muscle_mehg, colour=Age, shape=Sex))+
  geom_point(size=5)+
  scale_colour_brewer(palette = "YlOrRd")+
  geom_smooth(data=unstack_all_hg, aes(x=Feces_hg, y=Muscle_mehg, shape = NULL, color= NULL), 
              method='lm', formula= y~x, color = "black")+
  #stat_cor(data=unstack_all_hg, aes(x=Feces_hg, y=Muscle_mehg , shape = NULL, color= NULL, 
  #                                  label = after_stat(rr.label)), color = "black", geom = "label")+
  theme_minimal(base_size = 15)+
  xlab("Log10 Feces Hg (ug/g)")+
  ylab("Log10 Muscle MeHg (ug/g)")
p2

p3 <- ggplot(data=unstack_all_hg, aes(x=Feces_hg, y=Feces_mehg, colour=Age, shape=Sex))+
  geom_point(size=5)+
  scale_colour_brewer(palette = "YlOrRd")+
  theme_minimal(base_size = 15)+
  geom_smooth(data=unstack_all_hg, aes(x=Feces_hg, y=Feces_mehg, shape = NULL, color= NULL), 
              method='lm', formula= y~x, color = "black")+
 # stat_cor(data=unstack_all_hg, aes(x=Feces_hg, y=Feces_mehg , shape = NULL, color= NULL, 
  #                                  label = after_stat(rr.label)), color = "black", geom = "label")+
  xlab("Log10 Feces Hg (ug/g)")+
  ylab("Log10 Feces MeHg (ug/g)")
p3

p4 <- ggplot(data=unstack_all_hg, aes(x=Feces_hg, y=Fat_mehg, colour=Age, shape=Sex))+
  geom_point(size=5)+
  scale_colour_brewer(palette = "YlOrRd")+
  geom_smooth(data=unstack_all_hg, aes(x=Feces_hg, y=Fat_mehg, shape = NULL, color= NULL), 
              method='lm', formula= y~x, color = "black")+
  #stat_cor(data=unstack_all_hg, aes(x=Feces_hg, y=Fat_mehg , shape = NULL, color= NULL, 
  #                                  label = after_stat(rr.label)), color = "black", geom = "label")+
  theme_minimal(base_size = 15)+
  xlab("Log10 Feces Hg (ug/g)")+
  ylab("Log10 Fat MeHg (ug/g)")
p4


#Compile
figure2=ggarrange(p1, p2, p3, p4,
                  labels = c("A", "B", "C", "D"),
                  vjust = 1,
                  hjust = -0.5,
                  ncol = 2, nrow = 2,
                  common.legend = TRUE,
                  legend = "bottom",
                  font.label = list(size = 16))
figure2
#Plot figures with dpi=300
cowplot::save_plot("mehg_hg_tissue.jpg", figure2, base_width = 10, base_height = 10, dpi = 300)

