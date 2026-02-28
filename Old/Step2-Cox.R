rm(list=ls())
library(tidyverse)
library(arrow)
source('functions.R')

# DataStep1 = read_feather('D:/UKBiobank/ML_Stroke/Data/memb_map_df_demographics.feather')
# colnames(DataStep1)
# 
# covarites = colnames(DataStep1)[c(87,88,89,90,96,97,98,
#                                   112:152,155,156:185,42:68,
#                                   73:86)]
# DataStep2 = DataStep1 %>% select(covarites)
# colnames(DataStep2)
# DataStep3 = DataStep2 %>%
#   mutate_at(c(2,3,5,6,7,24:34,
#               35:76,107:120), as.factor)

# DataStep1 = read_feather('D:/UKBiobank/ML_Stroke/Data/memb_map_df_demographics.feather')
DataStep1 = read_feather('D:/UKBiobank/ML_Stroke/Data/memb_map_df_all.feather')
colnames(DataStep1)
# colnames(DataStep2)

# covarites = colnames(DataStep1)[c(87,88,89,90,96,97,98,
#                                   112:152,155,156:185,42:68,
#                                   73:86)]

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


# Cox regression
fit_Stroke_Model1 = coxph(Surv(`ALL Stroke_event_time`, `ALL Stroke`) ~ 
              Cluster, data = DataStep3)
Stroke_Model1 = cox_result_out(fit = fit_Stroke_Model1) %>% 
  filter(term == 'Cluster2') %>% 
  cbind(Var = 'StrokeModel1')

fit_Stroke_Model2 = coxph(Surv(`ALL Stroke_event_time`, `ALL Stroke`) ~ 
                            Cluster+`Age at recruitment`+ `Sex`+ `Townsend deprivation index at recruitment`+`Smoking status`+
                            `Alcohol intake frequency` + `Systolic blood pressure, automated reading` +
                            `Diastolic blood pressure, automated reading` + `Body mass index (BMI)`, data = DataStep3)
Stroke_Model2 = cox_result_out(fit = fit_Stroke_Model2) %>% 
  filter(term == 'Cluster2') %>% 
  cbind(Var = 'StrokeModel2')



fit_Stroke_Model3 = coxph(Surv(`ALL Stroke_event_time`, `ALL Stroke`) ~ 
              Cluster+`Age at recruitment`+ `Sex`+ `Townsend deprivation index at recruitment`+`Smoking status`+
              `Alcohol intake frequency` + `Systolic blood pressure, automated reading` +
              `Diastolic blood pressure, automated reading` + `Body mass index (BMI)`+
              `Atrial fibrillation`+`Type 2 diabetes`+`Coronary heart disease`+`Heart failure`, data = DataStep3)
Stroke_Model3 = cox_result_out(fit = fit_Stroke_Model3) %>% 
  filter(term == 'Cluster2') %>% 
  cbind(Var = 'StrokeModel3')


fit_Death_Model1 = coxph(Surv(`death_allcause_event_time`, `death_allcause_event`) ~ 
              Cluster, data = DataStep3)
Death_Model1 = cox_result_out(fit = fit_Death_Model1) %>% 
  filter(term == 'Cluster2') %>% 
  cbind(Var = 'DeathModel1')

fit_Death_Model2 = coxph(Surv(`death_allcause_event_time`, `death_allcause_event`) ~ 
                           Cluster+`Age at recruitment`+ `Sex`+ `Townsend deprivation index at recruitment`+`Smoking status`+
                           `Alcohol intake frequency` + `Systolic blood pressure, automated reading` +
                           `Diastolic blood pressure, automated reading` + `Body mass index (BMI)`, data = DataStep3)
Death_Model2 = cox_result_out(fit = fit_Death_Model2) %>% 
  filter(term == 'Cluster2') %>% 
  cbind(Var = 'DeathModel2')

fit_Death_Model3 = coxph(Surv(`death_allcause_event_time`, `death_allcause_event`) ~ 
                           Cluster+`Age at recruitment`+ `Sex`+ `Townsend deprivation index at recruitment`+`Smoking status`+
                           `Alcohol intake frequency` + `Systolic blood pressure, automated reading` +
                           `Diastolic blood pressure, automated reading` + `Body mass index (BMI)`+
                           `Atrial fibrillation`+`Type 2 diabetes`+`Coronary heart disease`+`Heart failure`, data = DataStep3)
Death_Model3 = cox_result_out(fit = fit_Death_Model3) %>% 
  filter(term == 'Cluster2') %>% 
  cbind(Var = 'DeathModel3')


# competing risk

DataStep4 = DataStep3 %>% 
  mutate(comp_event = case_when(`death_allcause_event` == 0 & `ALL Stroke` == 0 ~ 0,
                                `death_allcause_event` == 0 & `ALL Stroke` == 1 ~ 1,
                                `death_allcause_event` == 1 & `ALL Stroke` == 0 ~ 2,
                                `death_allcause_event` == 1 & `ALL Stroke` == 1 & `ALL Stroke_event_time` <= `death_allcause_event_time` ~ 1,
                                `death_allcause_event` == 1 & `ALL Stroke` == 1 & `ALL Stroke_event_time` > `death_allcause_event_time` ~ 2),
         comp_event_time = case_when(`ALL Stroke_event_time` <= `death_allcause_event_time` ~ `ALL Stroke_event_time`,
                                     `ALL Stroke_event_time` > `death_allcause_event_time` ~ `death_allcause_event_time`))

library(cmprsk)
Comp_Stroke_Model1 = crr(ftime = DataStep4$comp_event_time, fstatus = DataStep4$comp_event, cov1 = DataStep4[, c("Cluster")])
Competing_Model = cox_result_out(fit = Comp_Stroke_Model1) %>% 
  filter(term == 'Cluster') %>% 
  cbind(Var = 'CompModel2')

Stroke_Model1 %>% 
  rbind(Stroke_Model2) %>% 
  rbind(Stroke_Model3) %>% 
  rbind(Death_Model1) %>% 
  rbind(Death_Model2) %>%
  rbind(Death_Model3) %>% 
  rbind(Competing_Model) %>% 
  write.csv('./tables/V2-Cox.csv')

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
DataStep4$comp_event2 <- factor(DataStep4$comp_event, 
                               levels = c(0, 1, 2), 
                               labels = c("None", "Stroke", "Death"))


fit2 <- survfit(Surv(comp_event_time, comp_event2, type="mstate") ~ Phenotype, data=DataStep4 %>% 
                  mutate(Phenotype = Cluster))
p = ggcompetingrisks(fit2,
                 ylim = c(0,0.05),
                 xlab = "Follow up, years",            # 自定义X轴标签
                 ylab = "Cumulative incidence",     # 自定义Y轴标签
                 title = "",  # 设置图表标题
                 palette = "Dark2",                # 使用配色方案 'Dark2'
                 font.title = c(16, "bold"),       # 设置标题字体大小和样式
                 font.x = c(14, "plain"),          # 设置X轴标签的字体大小和样式
                 font.y = c(14, "plain"),          # 设置Y轴标签的字体大小和样式
                 font.legend = c(12, "italic"),    # 设置图例字体大小和样式
                 risk.table = TRUE,                # 显示风险表
                 risk.table.col = "strata",        # 风险表按分层变量着色
                 linetype = "strata",              # 线型按分层变量区分
                 conf.int = TRUE)+
  theme_my

tiff('figures/V2-CompKM.tiff', units = "cm",width = 20, height = 12,
     res = 300)
p
dev.off()

fit_death <- survfit(Surv(death_allcause_event_time, death_allcause_event) ~ Cluster, data = DataStep3)
fit_stroke <- survfit(Surv(`ALL Stroke_event_time`, `ALL Stroke`) ~ Cluster, data = DataStep3)

km_Cluster <-
  km_plot(
    fit_death,
    DataStep3,
    title = 'B',
    p_x = 0.2,
    p_y = 0.975,
    p_x_add = 0.9,
    censor.size = 1 / 1.5*1.5,
    surv.plot.height = 0.15*1.5,
    size = 0.5*1.5,
    legend = c(0.25,0.525),
    palette = c("#98714D", "#0C8D66"),
    legend.labs = c('Phenotype 1','Phenotype 2'),
    tables.height = 0.3,
    fontsize = 2*1.5,
    font.legend = 6*1.5,
    text_size = 6*1.5,
    line_size = 0.1*1.5,
    p_size = 2*1.5,
    break.x.by = 2,
    ylim = c(0.95,1),
    ylab = 'Mortality free probability'
  )

km_Cluster$table <- customize_labels(
  km_Cluster$table,
  font.title    = c(6*1.5, "bold.italic", "black")
)

customize_labels(km_Cluster, font.caption = 2)

km_Cluster2 <-
  km_plot(
    fit_stroke,
    DataStep3,
    title = '',
    p_x = 0.2,
    p_y = 0.975,
    p_x_add = 0.9,
    censor.size = 1 / 1.5*1.5,
    surv.plot.height = 0.15*1.5,
    size = 0.5*1.5,
    legend = c(0.25,0.525),
    palette = c("#98714D", "#0C8D66"),
    legend.labs = c('Phenotype 1','Phenotype 2'),
    tables.height = 0.3,
    fontsize = 2*1.5,
    font.legend = 6*1.5,
    text_size = 6*1.5,
    line_size = 0.1*1.5,
    p_size = 2*1.5,
    break.x.by = 2,
    ylim = c(0.95,1),
    ylab = 'Stroke free probability'
  )

km_Cluster2$table <- customize_labels(
  km_Cluster2$table,
  font.title    = c(6*1.5, "bold.italic", "black")
)

customize_labels(km_Cluster2, font.caption = 2)

library(patchwork)
tiff('figures/V2-KM.tiff', units = "cm",width = 16, height = 9,
     res = 300)
km_Cluster2$plot
dev.off()




# Stroke subtypes SAH
fit_Stroke_Model1 = coxph(Surv(`ICD 10 - Subarachnoid haemorrhage_event_time`, `ICD 10 - Subarachnoid haemorrhage`) ~ 
                            Cluster, data = DataStep3)
SAH_Model1 = cox_result_out(fit = fit_Stroke_Model1) %>% 
  filter(term == 'Cluster2') %>% 
  cbind(Var = 'SAHModel1')

fit_Stroke_Model2 = coxph(Surv(`ICD 10 - Subarachnoid haemorrhage_event_time`, `ICD 10 - Subarachnoid haemorrhage`) ~ 
                            Cluster+`Age at recruitment`+ `Sex`+ `Townsend deprivation index at recruitment`+`Smoking status`+
                            `Alcohol intake frequency` + `Systolic blood pressure, automated reading` +
                            `Diastolic blood pressure, automated reading` + `Body mass index (BMI)`, data = DataStep3)
SAH_Model2 = cox_result_out(fit = fit_Stroke_Model2) %>% 
  filter(term == 'Cluster2') %>% 
  cbind(Var = 'SAHModel2')

fit_Stroke_Model3 = coxph(Surv(`ICD 10 - Subarachnoid haemorrhage_event_time`, `ICD 10 - Subarachnoid haemorrhage`) ~ 
                            Cluster+`Age at recruitment`+ `Sex`+ `Townsend deprivation index at recruitment`+`Smoking status`+
                            `Alcohol intake frequency` + `Systolic blood pressure, automated reading` +
                            `Diastolic blood pressure, automated reading` + `Body mass index (BMI)`+
                            `Atrial fibrillation`+`Type 2 diabetes`+`Coronary heart disease`+`Heart failure`, data = DataStep3)
SAH_Model3 = cox_result_out(fit = fit_Stroke_Model3) %>% 
  filter(term == 'Cluster2') %>% 
  cbind(Var = 'SAHModel3')

# Intracerebral haemorrhage
fit_Stroke_Model1 = coxph(Surv(`ICD 10 - Intracerebral haemorrhage_event_time`, `ICD 10 - Intracerebral haemorrhage`) ~ 
                            Cluster, data = DataStep3)
ICH_Model1 = cox_result_out(fit = fit_Stroke_Model1) %>% 
  filter(term == 'Cluster2') %>% 
  cbind(Var = 'ICHModel1')

fit_Stroke_Model2 = coxph(Surv(`ICD 10 - Intracerebral haemorrhage_event_time`, `ICD 10 - Intracerebral haemorrhage`) ~ 
                            Cluster+`Age at recruitment`+ `Sex`+ `Townsend deprivation index at recruitment`+`Smoking status`+
                            `Alcohol intake frequency` + `Systolic blood pressure, automated reading` +
                            `Diastolic blood pressure, automated reading` + `Body mass index (BMI)`, data = DataStep3)
ICH_Model2 = cox_result_out(fit = fit_Stroke_Model2) %>% 
  filter(term == 'Cluster2') %>% 
  cbind(Var = 'ICHModel2')

fit_Stroke_Model3 = coxph(Surv(`ICD 10 - Intracerebral haemorrhage_event_time`, `ICD 10 - Intracerebral haemorrhage`) ~ 
                            Cluster+`Age at recruitment`+ `Sex`+ `Townsend deprivation index at recruitment`+`Smoking status`+
                            `Alcohol intake frequency` + `Systolic blood pressure, automated reading` +
                            `Diastolic blood pressure, automated reading` + `Body mass index (BMI)`+
                            `Atrial fibrillation`+`Type 2 diabetes`+`Coronary heart disease`+`Heart failure`, data = DataStep3)
ICH_Model3 = cox_result_out(fit = fit_Stroke_Model3) %>% 
  filter(term == 'Cluster2') %>% 
  cbind(Var = 'ICHModel3')


# Cerebral infarction
fit_Stroke_Model1 = coxph(Surv(`ICD 10 - Cerebral infarction_event_time`, `ICD 10 - Cerebral infarction`) ~ 
                            Cluster, data = DataStep3)
CI_Model1 = cox_result_out(fit = fit_Stroke_Model1) %>% 
  filter(term == 'Cluster2') %>% 
  cbind(Var = 'CIModel1')

fit_Stroke_Model2 = coxph(Surv(`ICD 10 - Cerebral infarction_event_time`, `ICD 10 - Cerebral infarction`) ~ 
                            Cluster+`Age at recruitment`+ `Sex`+ `Townsend deprivation index at recruitment`+`Smoking status`+
                            `Alcohol intake frequency` + `Systolic blood pressure, automated reading` +
                            `Diastolic blood pressure, automated reading` + `Body mass index (BMI)`, data = DataStep3)
CI_Model2 = cox_result_out(fit = fit_Stroke_Model2) %>% 
  filter(term == 'Cluster2') %>% 
  cbind(Var = 'CIModel2')

fit_Stroke_Model3 = coxph(Surv(`ICD 10 - Cerebral infarction_event_time`, `ICD 10 - Cerebral infarction`) ~ 
                            Cluster+`Age at recruitment`+ `Sex`+ `Townsend deprivation index at recruitment`+`Smoking status`+
                            `Alcohol intake frequency` + `Systolic blood pressure, automated reading` +
                            `Diastolic blood pressure, automated reading` + `Body mass index (BMI)`+
                            `Atrial fibrillation`+`Type 2 diabetes`+`Coronary heart disease`+`Heart failure`, data = DataStep3)
CI_Model3 = cox_result_out(fit = fit_Stroke_Model3) %>% 
  filter(term == 'Cluster2') %>% 
  cbind(Var = 'CIModel3')


# Other nontraumatic intracranial haemorrhage
fit_Stroke_Model1 = coxph(Surv(`ICD 10 - Other nontraumatic intracranial haemorrhage_event_time`, `ICD 10 - Other nontraumatic intracranial haemorrhage`) ~ 
                            Cluster, data = DataStep3)
OICH_Model1 = cox_result_out(fit = fit_Stroke_Model1) %>% 
  filter(term == 'Cluster2') %>% 
  cbind(Var = 'OICHModel1')

fit_Stroke_Model2 = coxph(Surv(`ICD 10 - Other nontraumatic intracranial haemorrhage_event_time`, `ICD 10 - Other nontraumatic intracranial haemorrhage`) ~ 
                            Cluster+`Age at recruitment`+ `Sex`+ `Townsend deprivation index at recruitment`+`Smoking status`+
                            `Alcohol intake frequency` + `Systolic blood pressure, automated reading` +
                            `Diastolic blood pressure, automated reading` + `Body mass index (BMI)`, data = DataStep3)
OICH_Model2 = cox_result_out(fit = fit_Stroke_Model2) %>% 
  filter(term == 'Cluster2') %>% 
  cbind(Var = 'OICHModel2')

fit_Stroke_Model3 = coxph(Surv(`ICD 10 - Other nontraumatic intracranial haemorrhage_event_time`, `ICD 10 - Other nontraumatic intracranial haemorrhage`) ~ 
                            Cluster+`Age at recruitment`+ `Sex`+ `Townsend deprivation index at recruitment`+`Smoking status`+
                            `Alcohol intake frequency` + `Systolic blood pressure, automated reading` +
                            `Diastolic blood pressure, automated reading` + `Body mass index (BMI)`+
                            `Atrial fibrillation`+`Type 2 diabetes`+`Coronary heart disease`+`Heart failure`, data = DataStep3)
OICH_Model3 = cox_result_out(fit = fit_Stroke_Model3) %>% 
  filter(term == 'Cluster2') %>% 
  cbind(Var = 'OICHModel3')


# Stroke, not specified as haemorrhage or infarction
fit_Stroke_Model1 = coxph(Surv(`ICD 10 - Stroke, not specified as haemorrhage or infarction_event_time`, 
                               `ICD 10 - Stroke, not specified as haemorrhage or infarction`) ~ 
                            Cluster, data = DataStep3)
UnStroke_Model1 = cox_result_out(fit = fit_Stroke_Model1) %>% 
  filter(term == 'Cluster2') %>% 
  cbind(Var = 'UnStrokeModel1')

fit_Stroke_Model2 = coxph(Surv(`ICD 10 - Stroke, not specified as haemorrhage or infarction_event_time`, 
                               `ICD 10 - Stroke, not specified as haemorrhage or infarction`) ~ 
                            Cluster+`Age at recruitment`+ `Sex`+ `Townsend deprivation index at recruitment`+`Smoking status`+
                            `Alcohol intake frequency` + `Systolic blood pressure, automated reading` +
                            `Diastolic blood pressure, automated reading` + `Body mass index (BMI)`, data = DataStep3)
UnStroke_Model2 = cox_result_out(fit = fit_Stroke_Model2) %>% 
  filter(term == 'Cluster2') %>% 
  cbind(Var = 'UnStrokeModel2')

fit_Stroke_Model3 = coxph(Surv(`ICD 10 - Stroke, not specified as haemorrhage or infarction_event_time`, 
                               `ICD 10 - Stroke, not specified as haemorrhage or infarction`) ~ 
                            Cluster+`Age at recruitment`+ `Sex`+ `Townsend deprivation index at recruitment`+`Smoking status`+
                            `Alcohol intake frequency` + `Systolic blood pressure, automated reading` +
                            `Diastolic blood pressure, automated reading` + `Body mass index (BMI)`+
                            `Atrial fibrillation`+`Type 2 diabetes`+`Coronary heart disease`+`Heart failure`, data = DataStep3)
UnStroke_Model3 = cox_result_out(fit = fit_Stroke_Model3) %>% 
  filter(term == 'Cluster2') %>% 
  cbind(Var = 'UnStrokeModel3')


SAH_Model1 %>% 
  rbind(SAH_Model2) %>% 
  rbind(SAH_Model3) %>% 
  rbind(ICH_Model1) %>% 
  rbind(ICH_Model2) %>% 
  rbind(ICH_Model3) %>% 
  rbind(CI_Model1) %>% 
  rbind(CI_Model2) %>% 
  rbind(CI_Model3) %>% 
  rbind(OICH_Model1) %>% 
  rbind(OICH_Model2) %>% 
  rbind(OICH_Model3) %>% 
  rbind(UnStroke_Model1) %>% 
  rbind(UnStroke_Model2) %>% 
  rbind(UnStroke_Model3) %>% 
  write.csv('./tables/V2-Sup_Cox.csv')
