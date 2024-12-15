rm(list=ls())
library(tidyverse)
library(arrow)
source('functions.R')

# DataStep1 = read_feather('D:/UKBiobank/ML_Stroke/Data/memb_map_df_demographics.feather')
DataStep1 = read_feather('D:/UKBiobank/ML_Stroke/Data/memb_map_df_all.feather')

temp_cmr = read_feather('D:/UKBiobank/ML_Stroke/Data/temp_cmr.feather')

DataStep2 = DataStep1 %>% 
  left_join(temp_cmr, by = 'eid')

colnames(DataStep2)
# colnames(DataStep2)

# covarites = colnames(DataStep1)[c(87,88,89,90,96,97,98,
#                                   112:152,155,156:185,42:68,
#                                   73:86)]

covarites = colnames(DataStep2)[c(55,56,57,58,64,65,66,
                                  80:120,123,124:153,28:54,
                                  14:27,
                                  28:54,2:9, 156, 159,162,165,168,
                                  171,174,177,180,183,186,189,192,195,
                                  197:204)]
DataStep3 = DataStep2 %>% select(covarites)
colnames(DataStep3)
DataStep4 = DataStep3 %>%
  mutate_at(c(2,3,5,6,7,24:34,
              35:78,107:120), as.factor)

continous_variables <- DataStep4 %>% 
  dplyr::select(where(is.numeric)) %>% 
  colnames()

catVars <- DataStep4 %>% 
  dplyr::select(where(is.factor)) %>% 
  colnames()

myvars <- colnames(DataStep4)

#基线表1
tab1 <- CreateTableOne(vars = myvars, 
                       strata = "Cluster", 
                       data = DataStep4, 
                       factorVars = catVars,
                       addOverall = TRUE) # 增加overall

tabMat1 <- print(tab1,
                 catDigits = 1,contDigits = 1,
                 quote = FALSE, # 不显示引号
                 noSpaces = TRUE, # 删除用于在R控制台中对齐文本的空格
                 printToggle = FALSE)

write.csv(tabMat1, 'tables/V2-Table1.csv')


