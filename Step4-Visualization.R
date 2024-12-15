rm(list=ls())
library(tidyverse)
library(arrow)
source('functions.R')
library(ggplot2)
library(patchwork)

DataStep1 = read_feather('D:/UKBiobank/ML_Stroke/Data/memb_map_df_all.feather')
colnames(DataStep1)

covarites = colnames(DataStep1)[c(55,56,57,58,64,65,66,
                                  80:120,123,124:153,28:54,
                                  14:27,
                                  28:54,2:9, 156, 159,162,165,168,
                                  171,174,177,180,183,186,189,192,195)]
DataStep2 = DataStep1 %>% select(covarites)
colnames(DataStep2)
DataStep3 = DataStep2 %>%
  mutate_at(c(2,3,5,6,7,24:34,
              35:75,78), as.factor)


DataBox = DataStep3 %>% 
  select(Cluster, `lv_stroke_volume_f24102_2_0`:`la_ejection_fraction_f24113_2_0`)

theme_my = theme(panel.grid = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              axis.title = element_text(size = 14/2*1.25),
              axis.line = element_line(size = 1/2*1.25),
              axis.text = element_text(size = 14/2*1.25),
              axis.ticks =  element_line(size = 1/2*1.25),
              axis.text.y.right = element_text(color="grey70"), 
              axis.title.y.right = element_text(color="grey70"),
              plot.title = element_text(face = 'bold',
                                        size = 14/2*1.25))

A = ggplot(DataBox, aes(x = Cluster, y = `lv_stroke_volume_f24102_2_0`, fill = Cluster)) +
  geom_boxplot()+
  labs(title = "A", 
       x = "Phenotype", 
       y = "LV Stroke Volume (mL)") +
  theme_my +
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("1", "2"))) +
  scale_fill_manual(values = c("1" = "blue", "2" = "red"),
                    name = "Phenotype")

B = ggplot(DataBox, aes(x = Cluster, y = `lv_myocardial_mass_f24105_2_0`, fill = Cluster)) +
  geom_boxplot()+
  labs(title = "B", 
       x = "Phenotype", 
       y = "LV Myocardial Mass (g)") +
  theme_my +
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("1", "2"))) +
  scale_fill_manual(values = c("1" = "blue", "2" = "red"),
                    name = "Phenotype")

C = ggplot(DataBox, aes(x = Cluster, y = `lv_end_diastolic_volume_f24100_2_0`, fill = Cluster)) +
  geom_boxplot()+
  labs(title = "C", 
       x = "Phenotype", 
       y = "LV End Diastolic Volume (mL)") +
  theme_my +
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("1", "2"))) +
  scale_fill_manual(values = c("1" = "blue", "2" = "red"),
                    name = "Phenotype")


D = ggplot(DataBox, aes(x = Cluster, y = `lv_end_systolic_volume_f24101_2_0`, fill = Cluster)) +
  geom_boxplot()+
  labs(title = "D", 
       x = "Phenotype", 
       y = "LV End Systolic Volume (mL)") +
  theme_my +
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("1", "2"))) +
  scale_fill_manual(values = c("1" = "blue", "2" = "red"),
                    name = "Phenotype")

E = ggplot(DataBox, aes(x = Cluster, y = `lv_ejection_fraction_f24103_2_0`, fill = Cluster)) +
  geom_boxplot()+
  labs(title = "E", 
       x = "Phenotype", 
       y = "LV Ejection Fraction (%)") +
  theme_my +
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("1", "2"))) +
  scale_fill_manual(values = c("1" = "blue", "2" = "red"),
                    name = "Phenotype")

`F` = ggplot(DataBox, aes(x = Cluster, y = `lv_longitudinal_strain_global_f24181_2_0`, fill = Cluster)) +
  geom_boxplot()+
  labs(title = "D", 
       x = "Phenotype", 
       y = "LV Longitudinal Strain Global (%)") +
  theme_my +
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("1", "2"))) +
  scale_fill_manual(values = c("1" = "blue", "2" = "red"),
                    name = "Phenotype")

G = ggplot(DataBox, aes(x = Cluster, y = `la_maximum_volume_f24110_2_0`, fill = Cluster)) +
  geom_boxplot()+
  labs(title = "G", 
       x = "Phenotype", 
       y = "LA Maximum Volume (mL)") +
  theme_my +
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("1", "2"))) +
  scale_fill_manual(values = c("1" = "blue", "2" = "red"),
                    name = "Phenotype")

H = ggplot(DataBox, aes(x = Cluster, y = `la_ejection_fraction_f24113_2_0`, fill = Cluster)) +
  geom_boxplot()+
  labs(title = "H", 
       x = "Phenotype", 
       y = "LA Ejection Fraction (%)") +
  theme_my +
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("1", "2"))) +
  scale_fill_manual(values = c("1" = "blue", "2" = "red"),
                    name = "Phenotype")

tiff("figures/V2-Boxplot.tiff", units = "cm",width = 24, height = 18,
     res = 300)
A + B + `C` +D + E+`F`+ G+H + plot_layout(nrow = 2,guides = 'collect')
dev.off()
